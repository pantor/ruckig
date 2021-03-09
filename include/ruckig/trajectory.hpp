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

//! Interface for the generated trajectory.
template<size_t DOFs>
class Trajectory {
    constexpr static double eps {std::numeric_limits<double>::epsilon()};

public:
    //! Duration of the synchronized trajectory
    double duration;

    //! Set of current profiles for each DoF
    std::array<Profile, DOFs> profiles;

    //! Minimum duration of each independent DoF
    std::array<double, DOFs> independent_min_durations;

    //! Calculate time-optimal waypoint-based trajectories
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
            for (size_t i = 0; p.t_brakes[i] > 0 && i < 2; ++i) {
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

    //! Get the output parameter for the given time
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

    //! Get the min/max values of the position for each DoF and the current trajectory
    std::array<PositionExtrema, DOFs> get_position_extrema() {
        std::array<PositionExtrema, DOFs> result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = profiles[dof].get_position_extrema();
        }
        return result;
    }
};

} // namespace ruckig
