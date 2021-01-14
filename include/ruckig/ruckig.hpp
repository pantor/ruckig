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
#include <ruckig/profile.hpp>
#include <ruckig/steps.hpp>


namespace ruckig {

template<size_t DOFs, bool throw_error = false>
class Ruckig {
    //! Current input, only for comparison for recalculation
    InputParameter<DOFs> current_input;

    //! Normalized input for calculating the trajectory
    InputParameter<DOFs> normalized_input;

    //! Scale that normalizes the input
    double scale;

    //! Current time in [s]
    double t;
    
    //! Duration of the current trajectory in [s]
    double tf;

    //! Set of current profiles for each DoF
    std::array<Profile, DOFs> profiles;

    static bool abs_compare(double a, double b) {
        return (std::abs(a) < std::abs(b));
    }

    void normalize_input(const InputParameter<DOFs>& input) {
        const auto [vMax_min, vMax_max] = std::minmax_element(input.max_velocity.cbegin(), input.max_velocity.cend(), abs_compare);
        const auto [aMax_min, aMax_max] = std::minmax_element(input.max_acceleration.cbegin(), input.max_acceleration.cend(), abs_compare);
        const auto [jMax_min, jMax_max] = std::minmax_element(input.max_jerk.cbegin(), input.max_jerk.cend(), abs_compare);

        const double min_value = std::min({*vMax_min, *aMax_min, *jMax_min});
        const double max_value = std::max({*vMax_max, *aMax_max, *jMax_max});
        scale = 1.0; // / min_value;

        normalized_input = input;
        // normalized_input.scale(scale);
    }

    bool synchronize(const std::array<Block, DOFs>& blocks, std::optional<double> t_min, double& t_sync, size_t& limiting_dof, std::array<Profile, DOFs>& profiles) {
        if (DOFs == 1 && !t_min) {
            limiting_dof = 0;
            t_sync = blocks[0].t_min;
            profiles[limiting_dof] = blocks[0].p_min;
            return true;
        }

        // Possible t_syncs are the start times of the intervals
        std::array<double, 3*DOFs> possible_t_syncs;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            const auto& block = blocks[dof];
            possible_t_syncs[3 * dof] = block.t_min;
            possible_t_syncs[3 * dof + 1] = block.a ? block.a->right : std::numeric_limits<double>::infinity();
            possible_t_syncs[3 * dof + 2] = block.b ? block.b->right : std::numeric_limits<double>::infinity();
        }

        // Test them in sorted order
        std::array<size_t, 3*DOFs> idx;
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&possible_t_syncs](size_t i, size_t j) { return possible_t_syncs[i] < possible_t_syncs[j]; });

        for (size_t i: idx) {
            const double possible_t_sync = possible_t_syncs[i];
            if (std::any_of(blocks.begin(), blocks.end(), [possible_t_sync](auto block){ return block.is_blocked(possible_t_sync); }) || possible_t_sync < t_min.value_or(0.0)) {
                continue;
            }

            t_sync = possible_t_sync;
            limiting_dof = std::ceil((i + 1.0) / 3) - 1;
            switch (i % 3) {
                case 0: {
                    profiles[limiting_dof] = blocks[limiting_dof].p_min;
                } break;
                case 1: {
                    profiles[limiting_dof] = blocks[limiting_dof].p_a.value();
                } break;
                case 2: {
                    profiles[limiting_dof] = blocks[limiting_dof].p_b.value();
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

        normalize_input(input);
        InputParameter<DOFs>& inp = normalized_input;

        std::array<Block, DOFs> blocks;
        std::array<double, DOFs> p0s, v0s, a0s; // Starting point of profiles without brake trajectory
        std::array<double, DOFs> inp_min_velocity;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (!inp.enabled[dof]) {
                continue;
            }

            Profile& p = profiles[dof];

            // Calculate brakes (if input exceeds or will exceed limits)
            inp_min_velocity[dof] = inp.min_velocity ? inp.min_velocity.value()[dof] : -inp.max_velocity[dof];
            Brake::get_brake_trajectory(inp.current_velocity[dof], inp.current_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp.max_jerk[dof], p.t_brakes, p.j_brakes);
            
            p.t_brake = p.t_brakes[0] + p.t_brakes[1];

            p0s[dof] = inp.current_position[dof];
            v0s[dof] = inp.current_velocity[dof];
            a0s[dof] = inp.current_acceleration[dof];

            if (p.t_brakes[0] > 0.0) {
                p.p_brakes[0] = p0s[dof];
                p.v_brakes[0] = v0s[dof];
                p.a_brakes[0] = a0s[dof];
                std::tie(p0s[dof], v0s[dof], a0s[dof]) = Profile::integrate(p.t_brakes[0], p0s[dof], v0s[dof], a0s[dof], p.j_brakes[0]);

                if (p.t_brakes[1] > 0.0) {
                    p.p_brakes[1] = p0s[dof];
                    p.v_brakes[1] = v0s[dof];
                    p.a_brakes[1] = a0s[dof];
                    std::tie(p0s[dof], v0s[dof], a0s[dof]) = Profile::integrate(p.t_brakes[1], p0s[dof], v0s[dof], a0s[dof], p.j_brakes[1]);
                }
            }

            Step1 step1 {p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp.max_jerk[dof]};
            bool found_profile = step1.get_profile(p);
            if (!found_profile) {
                if constexpr (throw_error) {
                    throw std::runtime_error("[ruckig] error in step 1, dof: " + std::to_string(dof) + " input: " + input.to_string());
                }
                return Result::ErrorExecutionTimeCalculation;
            }

            blocks[dof] = step1.block;
            output.independent_min_durations[dof] = step1.block.t_min;
        }

        size_t limiting_dof; // The DoF that doesn't need step 2
        bool found_synchronization = synchronize(blocks, inp.minimum_duration, tf, limiting_dof, profiles);
        if (!found_synchronization) {
            if constexpr (throw_error) {
                throw std::runtime_error("[ruckig] error in time synchronization: " + std::to_string(tf));
            }
            return Result::ErrorSynchronizationCalculation;
        }
        
        if (tf > 0.0) {
            for (size_t dof = 0; dof < DOFs; ++dof) {
                if (!inp.enabled[dof] || dof == limiting_dof) {
                    continue;
                }

                const double t_profile = tf - profiles[dof].t_brake.value_or(0.0);

                Step2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], inp.target_position[dof], inp.target_velocity[dof], inp.target_acceleration[dof], inp.max_velocity[dof], inp_min_velocity[dof], inp.max_acceleration[dof], inp.max_jerk[dof]};
                bool found_time_synchronization = step2.get_profile(profiles[dof]);
                if (!found_time_synchronization) {
                    if constexpr (throw_error) {
                        throw std::runtime_error("[ruckig] error in step 2 in dof: " + std::to_string(dof) + " for t sync: " + std::to_string(tf) + " input: " + input.to_string());
                    }
                    return Result::ErrorSynchronizationCalculation;
                }
            }
        }

        t = 0.0;
        output.duration = tf;
        output.new_calculation = true;
        return Result::Working;
    }

public:
    // Just a shorter notation
    using Input = InputParameter<DOFs>;
    using Output = OutputParameter<DOFs>;

    //! Time step between updates (cycle time) in [s]
    const double delta_time;

    explicit Ruckig(double delta_time): delta_time(delta_time) { }

    bool validate_input(const InputParameter<DOFs>& input) {
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

            if (input.target_acceleration[dof] > input.max_acceleration[dof]) {
                return false;
            }

            double max_target_acceleration = std::sqrt(2 * input.max_jerk[dof] * (input.max_velocity[dof] - std::abs(input.target_velocity[dof])));
            if (std::abs(input.target_acceleration[dof]) > max_target_acceleration) {
                return false;
            }
        }

        return true;
    }

    Result update(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        auto start = std::chrono::high_resolution_clock::now();

        t += delta_time;
        output.new_calculation = false;

        if (input != current_input && Result::Working != calculate(input, output)) {
            return Result::Error;
        }

        at_time(t, output);

        auto stop = std::chrono::high_resolution_clock::now();
        output.calculation_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000.0;

        if (t + delta_time > tf) {
            return Result::Finished;
        }

        current_input.current_position = output.new_position;
        current_input.current_velocity = output.new_velocity;
        current_input.current_acceleration = output.new_acceleration;
        return Result::Working;
    }

    void at_time(double time, OutputParameter<DOFs>& output) const {
        if (time + delta_time > tf) {
            // Keep constant acceleration
            for (size_t dof = 0; dof < DOFs; ++dof) {
                std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(time - tf, normalized_input.target_position[dof], normalized_input.target_velocity[dof], normalized_input.target_acceleration[dof], 0, 1./scale);
            }
            return;
        }

        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (!normalized_input.enabled[dof]) {
                // Keep constant acceleration
                std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(time, normalized_input.current_position[dof], normalized_input.current_velocity[dof], normalized_input.current_acceleration[dof], 0, 1./scale);
            }

            const auto& p = profiles[dof];

            double t_diff = time;
            if (p.t_brake) {
                if (t_diff < p.t_brake.value()) {
                    const size_t index = (t_diff < p.t_brakes[0]) ? 0 : 1;
                    if (index > 0) {
                        t_diff -= p.t_brakes[index - 1];
                    }

                    std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(t_diff, p.p_brakes[index], p.v_brakes[index], p.a_brakes[index], p.j_brakes[index], 1./scale);
                    continue;
                } else {
                    t_diff -= p.t_brake.value();
                }
            }

            // Non-time synchronization
            if (t_diff >= p.t_sum[6]) {
                // Keep constant acceleration
                std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(t_diff - p.t_sum[6], normalized_input.target_position[dof], normalized_input.target_velocity[dof], normalized_input.target_acceleration[dof], 0, 1./scale);
                continue;
            }

            const auto index_ptr = std::upper_bound(p.t_sum.begin(), p.t_sum.end(), t_diff);
            const size_t index = std::distance(p.t_sum.begin(), index_ptr);

            if (index > 0) {
                t_diff -= p.t_sum[index - 1];
            }

            std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(t_diff, p.p[index], p.v[index], p.a[index], p.j[index], 1./scale);
        }
    }

    double get_time() const {
        return t;
    }
};

} // namespace ruckig
