#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <tuple>

#include <ruckig/block.hpp>
#include <ruckig/brake.hpp>
#include <ruckig/input_parameter.hpp>
#include <ruckig/profile.hpp>
#include <ruckig/position.hpp>
#include <ruckig/velocity.hpp>


namespace ruckig {

// Forward declare alternative OTG algorithms for friend class
template <size_t> class Reflexxes;


//! Interface for the generated trajectory.
template<size_t DOFs>
class Trajectory {
    // Allow alternative OTG algorithms to directly access members (i.e. duration)
    friend class Reflexxes<DOFs>;

    constexpr static double eps {std::numeric_limits<double>::epsilon()};

    //! Set of current profiles for each DoF
    std::array<Profile, DOFs> profiles;

    double duration;
    std::array<double, DOFs> independent_min_durations;

    //! Is the trajectory phase synchronizable?
    static bool is_phase_synchronizable(
        const InputParameter<DOFs>& inp,
        const std::array<double, DOFs>& vMax,
        const std::array<double, DOFs>& vMin,
        const std::array<double, DOFs>& aMax,
        const std::array<double, DOFs>& aMin,
        const std::array<double, DOFs>& jMax,
        Profile::Direction limiting_direction,
        size_t limiting_dof,
        std::array<double, DOFs>& new_max_velocity,
        std::array<double, DOFs>& new_max_acceleration,
        std::array<double, DOFs>& new_min_acceleration,
        std::array<double, DOFs>& new_max_jerk
    ) {
        using Direction = Profile::Direction;

        // Get scaling factor of first DoF
        std::array<double, DOFs> pd;

        bool pd_found_nonzero {false};
        double v0_scale, a0_scale, vf_scale, af_scale;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            pd[dof] = inp.target_position[dof] - inp.current_position[dof];

            if (!pd_found_nonzero && std::abs(pd[dof]) > eps) {
                v0_scale = inp.current_velocity[dof] / pd[dof];
                a0_scale = inp.current_acceleration[dof] / pd[dof];
                vf_scale = inp.target_velocity[dof] / pd[dof];
                af_scale = inp.target_acceleration[dof] / pd[dof];
                pd_found_nonzero = true;
            }
        }

        if (!pd_found_nonzero) { // position difference is zero everywhere...
            return false;
        }

        const double max_jerk_limiting = (limiting_direction == Direction::UP) ? jMax[limiting_dof] : -jMax[limiting_dof];
        const double max_vel_limiting = (limiting_direction == Direction::UP) ? vMax[limiting_dof] : vMin[limiting_dof];
        const double max_acc_limiting = (limiting_direction == Direction::UP) ? aMax[limiting_dof] : aMin[limiting_dof];
        const double min_acc_limiting = (limiting_direction == Direction::UP) ? aMin[limiting_dof] : aMax[limiting_dof];

        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (dof == limiting_dof) {
                continue;
            }

            const double scale = pd[dof] / pd[limiting_dof];
            const double eps_colinear {10 * eps};

            // Are the vectors colinear?
            if (
                std::abs(inp.current_velocity[dof] - v0_scale * pd[dof]) > eps_colinear
                || std::abs(inp.current_acceleration[dof] - a0_scale * pd[dof]) > eps_colinear
                || std::abs(inp.target_velocity[dof] - vf_scale * pd[dof]) > eps_colinear
                || std::abs(inp.target_acceleration[dof] - af_scale * pd[dof]) > eps_colinear
                || std::abs(scale) > 1.0
            ) {
                return false;
            }

            // Are the old kinematic limits met?
            const Direction new_direction = ((limiting_direction == Direction::UP && scale >= 0.0) || (limiting_direction == Direction::DOWN && scale <= 0.0)) ? Direction::UP : Direction::DOWN;
            const double old_max_jerk = (new_direction == Direction::UP) ? jMax[dof] : -jMax[dof];
            const double old_max_vel = (new_direction == Direction::UP) ? vMax[dof] : vMin[dof];
            const double old_max_acc = (new_direction == Direction::UP) ? aMax[dof] : aMin[dof];
            const double old_min_acc = (new_direction == Direction::UP) ? aMin[dof] : aMax[dof];

            new_max_velocity[dof] = scale * max_vel_limiting;
            new_max_acceleration[dof] = scale * max_acc_limiting;
            new_min_acceleration[dof] = scale * min_acc_limiting;
            new_max_jerk[dof] = scale * max_jerk_limiting;

            if (
                std::abs(old_max_vel) < std::abs(new_max_velocity[dof])
                || std::abs(old_max_acc) < std::abs(new_max_acceleration[dof])
                || std::abs(old_min_acc) < std::abs(new_max_acceleration[dof])
                || std::abs(old_max_jerk) < std::abs(new_max_jerk[dof])
            ) {
                return false;
            }
        }

        return true;
    }

public:
    //! Calculate the time-optimal waypoint-based trajectory
    template<bool throw_error, bool return_error_at_maximal_duration>
    Result calculate(const InputParameter<DOFs>& inp, double delta_time) {
        std::array<Block, DOFs> blocks;
        std::array<double, DOFs> p0s, v0s, a0s; // Starting point of profiles without brake trajectory
        std::array<double, DOFs> inp_min_velocity, inp_min_acceleration;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            auto& p = profiles[dof];
            p.pf = inp.current_position[dof];
            p.vf = inp.current_velocity[dof];
            p.af = inp.current_acceleration[dof];

            if (!inp.enabled[dof]) {
                p.t_sum[6] = 0.0;
                continue;
            }

            inp_min_velocity[dof] = inp.min_velocity ? inp.min_velocity.value()[dof] : -inp.max_velocity[dof];
            inp_min_acceleration[dof] = inp.min_acceleration ? inp.min_acceleration.value()[dof] : -inp.max_acceleration[dof];

            // Calculate brake (if input exceeds or will exceed limits)
            switch (inp.interface) {
                case Interface::Position: {
                    Brake::get_position_brake_trajectory(inp.current_velocity[dof], inp.current_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof], p.t_brakes, p.j_brakes);
                } break;
                case Interface::Velocity: {
                    Brake::get_velocity_brake_trajectory(inp.current_acceleration[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof], p.t_brakes, p.j_brakes);
                } break;
            }

            p.t_brake = p.t_brakes[0] + p.t_brakes[1];
            p0s[dof] = inp.current_position[dof];
            v0s[dof] = inp.current_velocity[dof];
            a0s[dof] = inp.current_acceleration[dof];

            // Integrate brake pre-trajectory
            for (size_t i = 0; i < 2 && p.t_brakes[i] > 0; ++i) {
                p.p_brakes[i] = p0s[dof];
                p.v_brakes[i] = v0s[dof];
                p.a_brakes[i] = a0s[dof];
                std::tie(p0s[dof], v0s[dof], a0s[dof]) = Profile::integrate(p.t_brakes[i], p0s[dof], v0s[dof], a0s[dof], p.j_brakes[i]);
            }

            bool found_profile;
            switch (inp.interface) {
                case Interface::Position: {
                    PositionStep1 step1 {p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                    found_profile = step1.get_profile(p, blocks[dof]);
                } break;
                case Interface::Velocity: {
                    VelocityStep1 step1 {p0s[dof], v0s[dof], a0s[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                    found_profile = step1.get_profile(p, blocks[dof]);
                } break;
            }

            if (!found_profile) {
                if constexpr (throw_error) {
                    throw std::runtime_error("[ruckig] error in step 1, dof: " + std::to_string(dof) + " input: " + inp.to_string());
                }
                return Result::ErrorExecutionTimeCalculation;
            }

            independent_min_durations[dof] = blocks[dof].t_min;
        }

        int limiting_dof; // The DoF that doesn't need step 2
        const bool discrete_duration = (inp.duration_discretization == DurationDiscretization::Discrete);
        const bool found_synchronization = Block::synchronize<DOFs>(blocks, inp.minimum_duration, duration, limiting_dof, profiles, discrete_duration, delta_time);
        if (!found_synchronization) {
            if constexpr (throw_error) {
                throw std::runtime_error("[ruckig] error in time synchronization: " + std::to_string(duration));
            }
            return Result::ErrorSynchronizationCalculation;
        }

        if constexpr (return_error_at_maximal_duration) {
            if (duration > 7.6e3) {
                return Result::ErrorTrajectoryDuration;
            }
        }

        if (duration == 0.0) {
            return Result::Working;
        }

        if (inp.synchronization == Synchronization::None) {
            for (size_t dof = 0; dof < DOFs; ++dof) {
                if (!inp.enabled[dof] || dof == limiting_dof) {
                    continue;
                }

                profiles[dof] = blocks[dof].p_min;
            }
            return Result::Working;
        }

        if (inp.synchronization == Synchronization::Phase && inp.interface == Interface::Position) {
            std::array<double, DOFs> new_max_velocity, new_max_acceleration, new_min_acceleration, new_max_jerk;
            if (is_phase_synchronizable(inp, inp.max_velocity, inp_min_velocity, inp.max_acceleration, inp_min_acceleration, inp.max_jerk, profiles[limiting_dof].direction, limiting_dof, new_max_velocity, new_max_acceleration, new_min_acceleration, new_max_jerk)) {
                bool found_time_synchronization {true};
                for (size_t dof = 0; dof < DOFs; ++dof) {
                    if (!inp.enabled[dof] || dof == limiting_dof) {
                        continue;
                    }

                    Profile& p = profiles[dof];
                    const double t_profile = duration - p.t_brake.value_or(0.0);

                    p.t = profiles[limiting_dof].t; // Copy timing information from limiting DoF
                    p.set_boundary(inp.current_position[dof], inp.current_velocity[dof], inp.current_acceleration[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof]);

                    // Profile::Limits::NONE is a small hack, as there is no specialization for that in the check function
                    switch (p.jerk_signs) {
                        case Profile::JerkSigns::UDDU: {
                            if (!p.check<Profile::JerkSigns::UDDU, Profile::Limits::NONE>(t_profile, new_max_jerk[dof], new_max_velocity[dof], new_max_acceleration[dof], new_min_acceleration[dof])) {
                                found_time_synchronization = false;
                            }
                        } break;
                        case Profile::JerkSigns::UDUD: {
                            if (!p.check<Profile::JerkSigns::UDUD, Profile::Limits::NONE>(t_profile, new_max_jerk[dof], new_max_velocity[dof], new_max_acceleration[dof], new_min_acceleration[dof])) {
                                found_time_synchronization = false;
                            }
                        } break;
                    }

                    p.limits = profiles[limiting_dof].limits; // After check method call to set correct limits
                }

                if (found_time_synchronization) {
                    return Result::Working;
                }
            }
        }

        // The general case
        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (!inp.enabled[dof] || dof == limiting_dof) {
                continue;
            }

            Profile& p = profiles[dof];
            const double t_profile = duration - p.t_brake.value_or(0.0);

            if (inp.synchronization == Synchronization::TimeIfNecessary && std::abs(inp.target_velocity[dof]) < eps && std::abs(inp.target_acceleration[dof]) < eps) {
                p = blocks[dof].p_min;
                continue;
            }

            // Check if the final time corresponds to an extremal profile calculated in step 1
            if (std::abs(t_profile - blocks[dof].t_min) < eps) {
                p = blocks[dof].p_min;
                continue;
            } else if (blocks[dof].a && std::abs(t_profile - blocks[dof].a->right) < eps) {
                p = blocks[dof].a->profile;
                continue;
            } else if (blocks[dof].b && std::abs(t_profile - blocks[dof].b->right) < eps) {
                p = blocks[dof].b->profile;
                continue;
            }

            bool found_time_synchronization;
            switch (inp.interface) {
                case Interface::Position: {
                    PositionStep2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                    found_time_synchronization = step2.get_profile(p);
                } break;
                case Interface::Velocity: {
                    VelocityStep2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                    found_time_synchronization = step2.get_profile(p);
                } break;
            }
            if (!found_time_synchronization) {
                if constexpr (throw_error) {
                    throw std::runtime_error("[ruckig] error in step 2 in dof: " + std::to_string(dof) + " for t sync: " + std::to_string(duration) + " input: " + inp.to_string());
                }
                return Result::ErrorSynchronizationCalculation;
            }
            // std::cout << dof << " profile step2: " << p.to_string() << std::endl;
        }

        return Result::Working;
    }

    //! Get the kinematic state at a given time
    void at_time(double time, std::array<double, DOFs>& new_position, std::array<double, DOFs>& new_velocity, std::array<double, DOFs>& new_acceleration) const {
        if (time >= duration) {
            // Keep constant acceleration
            for (size_t dof = 0; dof < DOFs; ++dof) {
                std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(time - duration, profiles[dof].pf, profiles[dof].vf, profiles[dof].af, 0);
            }
            return;
        }

        for (size_t dof = 0; dof < DOFs; ++dof) {
            const Profile& p = profiles[dof];

            double t_diff = time;
            if (p.t_brake) {
                if (t_diff < p.t_brake.value()) {
                    const size_t index = (t_diff < p.t_brakes[0]) ? 0 : 1;
                    if (index > 0) {
                        t_diff -= p.t_brakes[index - 1];
                    }

                    std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff, p.p_brakes[index], p.v_brakes[index], p.a_brakes[index], p.j_brakes[index]);
                    continue;
                } else {
                    t_diff -= p.t_brake.value();
                }
            }

            // Non-time synchronization
            if (t_diff >= p.t_sum[6]) {
                // Keep constant acceleration
                std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff - p.t_sum[6], p.pf, p.vf, p.af, 0);
                continue;
            }

            const auto index_ptr = std::upper_bound(p.t_sum.begin(), p.t_sum.end(), t_diff);
            const size_t index = std::distance(p.t_sum.begin(), index_ptr);

            if (index > 0) {
                t_diff -= p.t_sum[index - 1];
            }

            std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff, p.p[index], p.v[index], p.a[index], p.j[index]);
        }
    }

    //! Get the duration of the (synchronized) trajectory
    double get_duration() const {
        return duration;
    }

    //! Get the minimum duration of each independent DoF
    std::array<double, DOFs> get_independent_min_durations() const {
        return independent_min_durations;
    }

    //! Get the min/max values of the position for each DoF
    std::array<PositionExtrema, DOFs> get_position_extrema() const {
        std::array<PositionExtrema, DOFs> result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = profiles[dof].get_position_extrema();
        }
        return result;
    }
};

} // namespace ruckig
