#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <tuple>
#include <type_traits>

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
    template<class T> using Vector = typename std::conditional<DOFs >= 1, std::array<T, DOFs>, std::vector<T>>::type;
    template<class T> using VectorIntervals = typename std::conditional<DOFs >= 1, std::array<T, 3*DOFs+1>, std::vector<T>>::type;

    // Allow alternative OTG algorithms to directly access members (i.e. duration)
    friend class Reflexxes<DOFs>;

    constexpr static double eps {std::numeric_limits<double>::epsilon()};

    // Set of current profiles for each DoF
    Vector<Profile> profiles;

    double duration {0.0};
    Vector<double> independent_min_durations;

    Vector<double> pd;

    VectorIntervals<double> possible_t_syncs;
    VectorIntervals<int> idx;

    Vector<Block> blocks;
    Vector<double> p0s, v0s, a0s; // Starting point of profiles without brake trajectory
    Vector<double> inp_min_velocity, inp_min_acceleration;

    Vector<double> new_max_velocity, new_min_velocity, new_max_acceleration, new_min_acceleration, new_max_jerk; // For phase synchronization

    Vector<PositionExtrema> position_extrema;

    //! Is the trajectory phase synchronizable?
    bool is_phase_synchronizable(
        const InputParameter<DOFs>& inp,
        const Vector<double>& vMax, const Vector<double>& vMin,
        const Vector<double>& aMax, const Vector<double>& aMin,
        const Vector<double>& jMax,
        Profile::Direction limiting_direction,
        size_t limiting_dof,
        Vector<double>& new_max_velocity, Vector<double>& new_min_velocity,
        Vector<double>& new_max_acceleration, Vector<double>& new_min_acceleration,
        Vector<double>& new_max_jerk
    ) {
        using Direction = Profile::Direction;

        // Get scaling factor of first DoF
        bool pd_found_nonzero {false};
        double v0_scale, a0_scale, vf_scale, af_scale;
        for (size_t dof = 0; dof < pd.size(); ++dof) {
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
        const double min_vel_limiting = (limiting_direction == Direction::UP) ? vMin[limiting_dof] : vMax[limiting_dof];
        const double max_acc_limiting = (limiting_direction == Direction::UP) ? aMax[limiting_dof] : aMin[limiting_dof];
        const double min_acc_limiting = (limiting_direction == Direction::UP) ? aMin[limiting_dof] : aMax[limiting_dof];

        for (size_t dof = 0; dof < pd.size(); ++dof) {
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
            const double old_min_vel = (new_direction == Direction::UP) ? vMin[dof] : vMax[dof];
            const double old_max_acc = (new_direction == Direction::UP) ? aMax[dof] : aMin[dof];
            const double old_min_acc = (new_direction == Direction::UP) ? aMin[dof] : aMax[dof];

            new_max_velocity[dof] = scale * max_vel_limiting;
            new_min_velocity[dof] = scale * min_vel_limiting;
            new_max_acceleration[dof] = scale * max_acc_limiting;
            new_min_acceleration[dof] = scale * min_acc_limiting;
            new_max_jerk[dof] = scale * max_jerk_limiting;

            if (
                std::abs(old_max_vel) < std::abs(new_max_velocity[dof])
                || std::abs(old_min_vel) < std::abs(new_min_velocity[dof])
                || std::abs(old_max_acc) < std::abs(new_max_acceleration[dof])
                || std::abs(old_min_acc) < std::abs(new_min_acceleration[dof])
                || std::abs(old_max_jerk) < std::abs(new_max_jerk[dof])
            ) {
                return false;
            }
        }

        return true;
    }

    bool synchronize(const Vector<Block>& blocks, std::optional<double> t_min, double& t_sync, int& limiting_dof, Vector<Profile>& profiles, bool discrete_duration, double delta_time) {
        if (degrees_of_freedom == 1 && !t_min && !discrete_duration) {
            limiting_dof = 0;
            t_sync = blocks[0].t_min;
            profiles[0] = blocks[0].p_min;
            return true;
        }

        // Possible t_syncs are the start times of the intervals and optional t_min
        bool any_interval {false};
        for (size_t dof = 0; dof < degrees_of_freedom; ++dof) {
            possible_t_syncs[dof] = blocks[dof].t_min;
            possible_t_syncs[degrees_of_freedom + dof] = blocks[dof].a ? blocks[dof].a->right : std::numeric_limits<double>::infinity();
            possible_t_syncs[2 * degrees_of_freedom + dof] = blocks[dof].b ? blocks[dof].b->right : std::numeric_limits<double>::infinity();
            any_interval |= blocks[dof].a || blocks[dof].b;
        }
        possible_t_syncs[3 * degrees_of_freedom] = t_min.value_or(std::numeric_limits<double>::infinity());
        any_interval |= t_min.has_value();

        if (discrete_duration) {
            for (size_t i = 0; i < possible_t_syncs.size(); ++i) {
                possible_t_syncs[i] = std::ceil(possible_t_syncs[i] / delta_time) * delta_time;
            }
        }

        // Test them in sorted order
        auto idx_end = any_interval ? idx.end() : idx.begin() + degrees_of_freedom;
        std::iota(idx.begin(), idx_end, 0);
        std::sort(idx.begin(), idx_end, [&possible_t_syncs=possible_t_syncs](size_t i, size_t j) { return possible_t_syncs[i] < possible_t_syncs[j]; });

        // Start at last tmin (or worse)
        for (auto i = idx.begin() + degrees_of_freedom - 1; i != idx_end; ++i) {
            const double possible_t_sync = possible_t_syncs[*i];
            if (std::any_of(blocks.begin(), blocks.end(), [possible_t_sync](const Block& block){ return block.is_blocked(possible_t_sync); }) || possible_t_sync < t_min.value_or(0.0)) {
                continue;
            }

            t_sync = possible_t_sync;
            if (*i == 3*degrees_of_freedom) { // Optional t_min
                limiting_dof = -1;
                return true;
            }

            const auto div = std::div(*i, degrees_of_freedom);
            limiting_dof = div.rem;
            switch (div.quot) {
                case 0: {
                    profiles[limiting_dof] = blocks[limiting_dof].p_min;
                } break;
                case 1: {
                    profiles[limiting_dof] = blocks[limiting_dof].a->profile;
                } break;
                case 2: {
                    profiles[limiting_dof] = blocks[limiting_dof].b->profile;
                } break;
            }
            return true;
        }

        return false;
    }

public:
    size_t degrees_of_freedom;

    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    Trajectory(): degrees_of_freedom(DOFs) { }

    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    Trajectory(size_t dofs): degrees_of_freedom(dofs) {
        blocks.resize(dofs);
        p0s.resize(dofs);
        v0s.resize(dofs);
        a0s.resize(dofs);
        inp_min_velocity.resize(dofs);
        inp_min_acceleration.resize(dofs);
        profiles.resize(dofs);
        independent_min_durations.resize(dofs);
        pd.resize(dofs);

        new_max_velocity.resize(dofs);
        new_min_velocity.resize(dofs);
        new_max_acceleration.resize(dofs);
        new_min_acceleration.resize(dofs);
        new_max_jerk.resize(dofs);

        position_extrema.resize(dofs);

        possible_t_syncs.resize(3*dofs+1);
        idx.resize(3*dofs+1);
    }

    //! Calculate the time-optimal waypoint-based trajectory
    template<bool throw_error, bool return_error_at_maximal_duration>
    Result calculate(const InputParameter<DOFs>& inp, double delta_time) {
        for (size_t dof = 0; dof < profiles.size(); ++dof) {
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
            switch (inp.control_interface) {
                case ControlInterface::Position: {
                    Brake::get_position_brake_trajectory(inp.current_velocity[dof], inp.current_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof], p.t_brakes, p.j_brakes);
                } break;
                case ControlInterface::Velocity: {
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
            switch (inp.control_interface) {
                case ControlInterface::Position: {
                    PositionStep1 step1 {p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                    found_profile = step1.get_profile(p, blocks[dof]);
                } break;
                case ControlInterface::Velocity: {
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
            // std::cout << dof << " profile step1: " << blocks[dof].to_string() << std::endl;
        }

        int limiting_dof; // The DoF that doesn't need step 2
        const bool discrete_duration = (inp.duration_discretization == DurationDiscretization::Discrete);
        const bool found_synchronization = synchronize(blocks, inp.minimum_duration, duration, limiting_dof, profiles, discrete_duration, delta_time);
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
            for (size_t dof = 0; dof < blocks.size(); ++dof) {
                if (!inp.enabled[dof] || dof == limiting_dof) {
                    continue;
                }

                profiles[dof] = blocks[dof].p_min;
            }
            return Result::Working;
        }

        if (inp.synchronization == Synchronization::Phase && inp.control_interface == ControlInterface::Position) {
            if (is_phase_synchronizable(inp, inp.max_velocity, inp_min_velocity, inp.max_acceleration, inp_min_acceleration, inp.max_jerk, profiles[limiting_dof].direction, limiting_dof, new_max_velocity, new_min_velocity, new_max_acceleration, new_min_acceleration, new_max_jerk)) {
                bool found_time_synchronization {true};
                for (size_t dof = 0; dof < profiles.size(); ++dof) {
                    if (!inp.enabled[dof] || dof == limiting_dof) {
                        continue;
                    }

                    Profile& p = profiles[dof];
                    const double t_profile = duration - p.t_brake.value_or(0.0);

                    p.t = profiles[limiting_dof].t; // Copy timing information from limiting DoF
                    p.jerk_signs = profiles[limiting_dof].jerk_signs;
                    p.set_boundary(inp.current_position[dof], inp.current_velocity[dof], inp.current_acceleration[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof]);

                    // Profile::Limits::NONE is a small hack, as there is no specialization for that in the check function
                    switch (p.jerk_signs) {
                        case Profile::JerkSigns::UDDU: {
                            if (!p.check_with_timing<Profile::JerkSigns::UDDU, Profile::Limits::NONE>(t_profile, new_max_jerk[dof], new_max_velocity[dof], new_min_velocity[dof], new_max_acceleration[dof], new_min_acceleration[dof])) {
                                found_time_synchronization = false;
                            }
                        } break;
                        case Profile::JerkSigns::UDUD: {
                            if (!p.check_with_timing<Profile::JerkSigns::UDUD, Profile::Limits::NONE>(t_profile, new_max_jerk[dof], new_max_velocity[dof], new_min_velocity[dof], new_max_acceleration[dof], new_min_acceleration[dof])) {
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
        for (size_t dof = 0; dof < profiles.size(); ++dof) {
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
            switch (inp.control_interface) {
                case ControlInterface::Position: {
                    PositionStep2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                    found_time_synchronization = step2.get_profile(p);
                } break;
                case ControlInterface::Velocity: {
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

    //! Continue the trajectory calculation
    template<bool throw_error, bool return_error_at_maximal_duration>
    Result continue_calculation(const InputParameter<DOFs>&, double) {
        return Result::Error;
    }

    //! Get the kinematic state at a given time

    //! The Python wrapper takes `time` as an argument, and returns `new_position`, `new_velocity`, and `new_acceleration` instead.
    void at_time(double time, Vector<double>& new_position, Vector<double>& new_velocity, Vector<double>& new_acceleration) const {
        if constexpr (DOFs == 0) {
            if (degrees_of_freedom != new_position.size() || degrees_of_freedom != new_velocity.size() || degrees_of_freedom != new_acceleration.size()) {
                throw std::runtime_error("[ruckig] mismatch in degrees of freedom (vector size).");
            }
        }

        if (time >= duration) {
            // Keep constant acceleration
            for (size_t dof = 0; dof < profiles.size(); ++dof) {
                std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(time - duration, profiles[dof].pf, profiles[dof].vf, profiles[dof].af, 0);
            }
            return;
        }

        for (size_t dof = 0; dof < profiles.size(); ++dof) {
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

    //! Get the durations when the intermediate waypoints are reached
    std::vector<double> get_intermediate_durations() const {
        return {duration};
    }

    //! Get the minimum duration of each independent DoF
    Vector<double> get_independent_min_durations() const {
        return independent_min_durations;
    }

    //! Get the min/max values of the position for each DoF
    Vector<PositionExtrema> get_position_extrema() {
        for (size_t dof = 0; dof < profiles.size(); ++dof) {
            position_extrema[dof] = profiles[dof].get_position_extrema();
        }
        return position_extrema;
    }

    //! Get the time that this trajectory passes a specific position of a given DoF the first time

    //! If the position is passed, this method returns true, otherwise false
    //! The Python wrapper takes `dof` and `position` as arguments and returns `time` (or `None`) instead
    bool get_first_time_at_position(size_t dof, double position, double& time) const {
        if (dof >= degrees_of_freedom) {
            return false;
        }

        double v, a;
        return profiles[dof].get_first_state_at_position(position, time, v, a);
    }
};

} // namespace ruckig
