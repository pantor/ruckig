#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <limits>
#include <math.h>
#include <numeric>
#include <optional>
#include <tuple>

#include <ruckig/parameter.hpp>
#include <ruckig/trajectory.hpp>
#include <ruckig/steps.hpp>


namespace ruckig {

//! Main class for the Ruckig algorithm.
template<size_t DOFs, bool throw_error = false, bool return_error_at_maximal_duration = true>
class Ruckig {
    //! Current input, only for comparison for recalculation
    InputParameter<DOFs> current_input;

    static bool synchronize(const std::array<Block, DOFs>& blocks, std::optional<double> t_min, double& t_sync, int& limiting_dof, std::array<Profile, DOFs>& profiles) {
        if (DOFs == 1 && !t_min) {
            limiting_dof = 0;
            t_sync = blocks[0].t_min;
            profiles[0] = blocks[0].p_min;
            return true;
        }

        // Possible t_syncs are the start times of the intervals and optional t_min
        std::array<double, 3*DOFs+1> possible_t_syncs;
        std::array<int, 3*DOFs+1> idx;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            possible_t_syncs[3 * dof] = blocks[dof].t_min;
            possible_t_syncs[3 * dof + 1] = blocks[dof].a ? blocks[dof].a->right : std::numeric_limits<double>::infinity();
            possible_t_syncs[3 * dof + 2] = blocks[dof].b ? blocks[dof].b->right : std::numeric_limits<double>::infinity();
        }
        possible_t_syncs[3 * DOFs] = t_min.value_or(std::numeric_limits<double>::infinity());

        // Test them in sorted order
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&possible_t_syncs](size_t i, size_t j) { return possible_t_syncs[i] < possible_t_syncs[j]; });

        // Start at last tmin (or worse)
        for (auto i = idx.begin() + DOFs - 1; i != idx.end(); ++i) {
            const double possible_t_sync = possible_t_syncs[*i];
            if (std::any_of(blocks.begin(), blocks.end(), [possible_t_sync](auto block){ return block.is_blocked(possible_t_sync); }) || possible_t_sync < t_min.value_or(0.0)) {
                continue;
            }

            t_sync = possible_t_sync;
            if (*i == 3*DOFs) { // Optional t_min
                limiting_dof = -1;
                return true;
            }

            const auto div = std::div(*i, 3);
            limiting_dof = div.quot;
            switch (div.rem) {
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

    Result calculate(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        current_input = input;

        if (!validate_input(input)) {
            return Result::ErrorInvalidInput;
        }

        InputParameter<DOFs>& inp = current_input;
        Trajectory<DOFs>& trajectory = output.trajectory;

        std::array<Block, DOFs> blocks;
        std::array<double, DOFs> p0s, v0s, a0s; // Starting point of profiles without brake trajectory
        std::array<double, DOFs> inp_min_velocity, inp_min_acceleration;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            Profile& p = trajectory.profiles[dof];

            if (!inp.enabled[dof]) {
                p.pf = inp.current_position[dof];
                p.vf = inp.current_velocity[dof];
                p.af = inp.current_acceleration[dof];
                p.t_sum[6] = 0.0;
                continue;
            }

            inp_min_velocity[dof] = inp.min_velocity ? inp.min_velocity.value()[dof] : -inp.max_velocity[dof];
            inp_min_acceleration[dof] = inp.min_acceleration ? inp.min_acceleration.value()[dof] : -inp.max_acceleration[dof];

            // Calculate brake (if input exceeds or will exceed limits)
            Brake::get_brake_trajectory(inp.current_velocity[dof], inp.current_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof], p.t_brakes, p.j_brakes);

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

            Step1 step1 {p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
            bool found_profile = step1.get_profile(p, blocks[dof]);
            if (!found_profile) {
                if constexpr (throw_error) {
                    throw std::runtime_error("[ruckig] error in step 1, dof: " + std::to_string(dof) + " input: " + input.to_string());
                }
                return Result::ErrorExecutionTimeCalculation;
            }

            trajectory.independent_min_durations[dof] = blocks[dof].t_min;
        }

        int limiting_dof; // The DoF that doesn't need step 2
        bool found_synchronization = synchronize(blocks, inp.minimum_duration, trajectory.duration, limiting_dof, trajectory.profiles);
        if (!found_synchronization) {
            if constexpr (throw_error) {
                throw std::runtime_error("[ruckig] error in time synchronization: " + std::to_string(trajectory.duration));
            }
            return Result::ErrorSynchronizationCalculation;
        }

        if constexpr (return_error_at_maximal_duration) {
            if (trajectory.duration > 7.6e3) {
                return Result::ErrorTrajectoryDuration;
            }
        }

        if (trajectory.duration > 0.0) {
            for (size_t dof = 0; dof < DOFs; ++dof) {
                if (!inp.enabled[dof] || dof == limiting_dof) {
                    continue;
                }

                Profile& p = trajectory.profiles[dof];
                const double t_profile = trajectory.duration - p.t_brake.value_or(0.0);

                // Check if the final time corresponds to an extremal profile calculated in step 1
                if (std::abs(t_profile - blocks[dof].t_min) < std::numeric_limits<double>::epsilon()) {
                    p = blocks[dof].p_min;
                    continue;
                } else if (blocks[dof].a && std::abs(t_profile - blocks[dof].a->right) < std::numeric_limits<double>::epsilon()) {
                    p = blocks[dof].a->profile;
                    continue;
                } else if (blocks[dof].b && std::abs(t_profile - blocks[dof].b->right) < std::numeric_limits<double>::epsilon()) {
                    p = blocks[dof].b->profile;
                    continue;
                }

                Step2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp_min_acceleration[dof], inp.max_jerk[dof]};
                bool found_time_synchronization = step2.get_profile(p);
                if (!found_time_synchronization) {
                    if constexpr (throw_error) {
                        throw std::runtime_error("[ruckig] error in step 2 in dof: " + std::to_string(dof) + " for t sync: " + std::to_string(trajectory.duration) + " input: " + input.to_string());
                    }
                    return Result::ErrorSynchronizationCalculation;
                }
                // std::cout << dof << " profile step2: " << p.to_string() << std::endl;
            }
        }

        output.time = 0.0;
        output.new_calculation = true;
        return Result::Working;
    }

public:
    // Just a shorter notation
    using Input = InputParameter<DOFs>;
    using Output = OutputParameter<DOFs>;

    //! Time step between updates (cycle time) in [s]
    double delta_time;

    explicit Ruckig() { }
    explicit Ruckig(double delta_time): delta_time(delta_time) { }

    bool validate_input(const InputParameter<DOFs>& input) {
        if (input.type == InputParameter<DOFs>::Type::Velocity) {
            return false;
        }

        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (input.max_velocity[dof] <= std::numeric_limits<double>::min()) {
                return false;
            }

            if (input.min_velocity && input.min_velocity.value()[dof] >= -std::numeric_limits<double>::min()) {
                return false;
            }

            if (input.max_acceleration[dof] <= std::numeric_limits<double>::min()) {
                return false;
            }

            if (input.min_acceleration && input.min_acceleration.value()[dof] >= -std::numeric_limits<double>::min()) {
                return false;
            }

            if (input.max_jerk[dof] <= std::numeric_limits<double>::min()) {
                return false;
            }

            if (std::isnan(input.target_position[dof])) {
                return false;
            }

            if (input.min_velocity) {
                if (input.target_velocity[dof] > input.max_velocity[dof] || input.target_velocity[dof] < input.min_velocity.value()[dof]) {
                    return false;
                }

            } else {
                if (std::abs(input.target_velocity[dof]) > input.max_velocity[dof]) {
                    return false;
                }
            }

            if (input.min_acceleration) {
                if (input.target_acceleration[dof] > input.max_acceleration[dof] || input.target_acceleration[dof] < input.min_acceleration.value()[dof]) {
                    return false;
                }

            } else {
                if (std::abs(input.target_acceleration[dof]) > input.max_acceleration[dof]) {
                    return false;
                }
            }

            double max_target_acceleration;
            if (input.min_velocity && input.target_velocity[dof] < 0) {
                max_target_acceleration = std::sqrt(-2 * input.max_jerk[dof] * (input.min_velocity.value()[dof] - input.target_velocity[dof]));
            } else {
                max_target_acceleration = std::sqrt(2 * input.max_jerk[dof] * (input.max_velocity[dof] - std::abs(input.target_velocity[dof])));
            }
            if (std::abs(input.target_acceleration[dof]) > max_target_acceleration) {
                return false;
            }
        }

        return true;
    }

    Result update(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        auto start = std::chrono::high_resolution_clock::now();

        output.new_calculation = false;

        if (input != current_input) {
            auto result = calculate(input, output);
            if (result != Result::Working) {
                return result;
            }
        }

        output.time += delta_time;
        output.trajectory.at_time(output.time, output.new_position, output.new_velocity, output.new_acceleration);

        auto stop = std::chrono::high_resolution_clock::now();
        output.calculation_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000.0;

        if (output.time > output.trajectory.duration) {
            return Result::Finished;
        }

        current_input.current_position = output.new_position;
        current_input.current_velocity = output.new_velocity;
        current_input.current_acceleration = output.new_acceleration;
        return Result::Working;
    }
};

} // namespace ruckig
