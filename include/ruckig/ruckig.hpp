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


namespace ruckig {

//! Which times are possible for synchronization?
struct Block {
    struct Interval {
        double left, right; // [s]
    };

    double t_min; // [s]
    Profile p_min; // Save min profile so that it doesn't need to be recalculated in Step2

    std::optional<Interval> a, b; // Max. 2 intervals can be blocked
    std::optional<Profile> p_a, p_b;

    bool is_blocked(double t) const {
        return (t < t_min) || (a && a->left < t && t < a->right) || (b && b->left < t && t < b->right);
    }
};


//! Calculates (pre-) trajectory to get current state below the limits
class Brake {
    static constexpr double eps {3e-15};

    static void acceleration_brake(double v0, double a0, double vMax, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);
    static void velocity_brake(double v0, double a0, double vMax, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);

public:
    static void get_brake_trajectory(double v0, double a0, double vMax, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);
    static void get_brake_trajectory(double v0, double a0, double vMax, double vMin, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);
};


class Step1 {
    using Limits = Profile::Limits;
    using Teeth = Profile::Teeth;

    double p0, v0, a0;
    double pf, vf, af;
    double vMax, vMin, aMax, jMax;

    // Pre-calculated expressions
    double pd;
    double v0_v0, vf_vf;
    double a0_a0, af_af, aMax_aMax;
    double jMax_jMax;

    double a0_p3, a0_p4, a0_p5, a0_p6;
    double af_p3, af_p4, af_p5, af_p6;

    // Max 6 valid profiles
    std::array<Profile, 6> valid_profiles;
    size_t valid_profile_counter;

    void add_profile(Profile profile, Limits limits, double jMax);

    void time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax);
    void time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax);
    void time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax);
    void time_up_vel(Profile& profile, double vMax, double aMax, double jMax);
    void time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax);
    void time_up_acc1(Profile& profile, double vMax, double aMax, double jMax);
    void time_up_acc0(Profile& profile, double vMax, double aMax, double jMax);
    void time_up_none(Profile& profile, double vMax, double aMax, double jMax);

    void time_down_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax);
    void time_down_acc1_vel(Profile& profile, double vMax, double aMax, double jMax);
    void time_down_acc0_vel(Profile& profile, double vMax, double aMax, double jMax);
    void time_down_vel(Profile& profile, double vMax, double aMax, double jMax);
    void time_down_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax);
    void time_down_acc1(Profile& profile, double vMax, double aMax, double jMax);
    void time_down_acc0(Profile& profile, double vMax, double aMax, double jMax);
    void time_down_none(Profile& profile, double vMax, double aMax, double jMax);

    bool calculate_block();

public:
    Block block;

    explicit Step1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double vMin, double aMax, double jMax);

    bool get_profile(const Profile& input);
};


class Step2 {
    using Teeth = Profile::Teeth;

    double tf;
    double p0, v0, a0;
    double pf, vf, af;
    double vMax, vMin, aMax, jMax;

    // Pre-calculated expressions
    double pd;
    double tf_tf, tf_p3, tf_p4;
    double vd, vd_vd, v0_v0, vf_vf;
    double ad, ad_ad, a0_a0, af_af, aMax_aMax;
    double jMax_jMax;

    double a0_p3, a0_p4, a0_p5, a0_p6;
    double af_p3, af_p4, af_p5, af_p6;

    bool time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax);
    bool time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax);
    bool time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax);
    bool time_up_vel(Profile& profile, double vMax, double aMax, double jMax);
    bool time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax);
    bool time_up_acc1(Profile& profile, double vMax, double aMax, double jMax);
    bool time_up_acc0(Profile& profile, double vMax, double aMax, double jMax);
    bool time_up_none(Profile& profile, double vMax, double aMax, double jMax);

    bool time_down_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax);
    bool time_down_acc1_vel(Profile& profile, double vMax, double aMax, double jMax);
    bool time_down_acc0_vel(Profile& profile, double vMax, double aMax, double jMax);
    bool time_down_vel(Profile& profile, double vMax, double aMax, double jMax);
    bool time_down_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax);
    bool time_down_acc1(Profile& profile, double vMax, double aMax, double jMax);
    bool time_down_acc0(Profile& profile, double vMax, double aMax, double jMax);
    bool time_down_none(Profile& profile, double vMax, double aMax, double jMax);

public:
    explicit Step2(double tf, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double vMin, double aMax, double jMax);

    bool get_profile(Profile& profile);
};


template<size_t DOFs, bool THROW_ERROR = false>
class Ruckig {
    InputParameter<DOFs> current_input;

    double t, tf;
    std::array<Profile, DOFs> profiles;

    bool synchronize(const std::array<Block, DOFs>& blocks, std::optional<double> t_min, double& t_sync, int& limiting_dof, std::array<Profile, DOFs>& profiles) {
        if (DOFs == 1 && !t_min) {
            limiting_dof = 0;
            t_sync = blocks[0].t_min;
            profiles[limiting_dof] = blocks[0].p_min;
            return true;
        }

        // Possible t_syncs are the start times of the intervals
        std::array<double, 3*DOFs> possible_t_syncs;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            auto& block = blocks[dof];
            possible_t_syncs[3 * dof] = block.t_min;
            possible_t_syncs[3 * dof + 1] = block.a ? block.a->right : std::numeric_limits<double>::infinity();
            possible_t_syncs[3 * dof + 2] = block.b ? block.b->right : std::numeric_limits<double>::infinity();
        }

        // Test them in sorted order
        std::array<size_t, 3*DOFs> idx;
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&possible_t_syncs](size_t i, size_t j) { return possible_t_syncs[i] < possible_t_syncs[j]; });

        for (size_t i: idx) {
            double possible_t_sync = possible_t_syncs[i];
            if (std::any_of(blocks.begin(), blocks.end(), [possible_t_sync](auto block){ return block.is_blocked(possible_t_sync); }) || possible_t_sync < t_min.value_or(0.0)) {
                continue;
            }

            t_sync = possible_t_sync;
            limiting_dof = std::ceil((i + 1.0) / 3) - 1;
            // std::cout << "sync: " << limiting_dof << " " << i % 3 << " " << t_sync << std::endl;
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
        // auto start = std::chrono::high_resolution_clock::now();
        current_input = input;

        // std::cout << "reference: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000.0 << std::endl;

        if (!validate_input(input)) {
            return Result::ErrorInvalidInput;
        }

        // std::cout << "validate: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000.0 << std::endl;

        std::array<Block, DOFs> blocks;
        std::array<double, DOFs> p0s, v0s, a0s; // Starting point of profiles without brake trajectory
        std::array<double, DOFs> input_min_velocity;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (!input.enabled[dof]) {
                continue;
            }

            // Calculate brakes (if input exceeds or will exceed limits)
            input_min_velocity[dof] = input.min_velocity ? input.min_velocity.value()[dof] : -input.max_velocity[dof];
            Brake::get_brake_trajectory(input.current_velocity[dof], input.current_acceleration[dof], input.max_velocity[dof], input_min_velocity[dof], input.max_acceleration[dof], input.max_jerk[dof], profiles[dof].t_brakes, profiles[dof].j_brakes);
            
            profiles[dof].t_brake = profiles[dof].t_brakes[0] + profiles[dof].t_brakes[1];

            p0s[dof] = input.current_position[dof];
            v0s[dof] = input.current_velocity[dof];
            a0s[dof] = input.current_acceleration[dof];

            if (profiles[dof].t_brakes[0] > 0.0) {
                profiles[dof].p_brakes[0] = p0s[dof];
                profiles[dof].v_brakes[0] = v0s[dof];
                profiles[dof].a_brakes[0] = a0s[dof];
                std::tie(p0s[dof], v0s[dof], a0s[dof]) = Profile::integrate(profiles[dof].t_brakes[0], p0s[dof], v0s[dof], a0s[dof], profiles[dof].j_brakes[0]);

                if (profiles[dof].t_brakes[1] > 0.0) {
                    profiles[dof].p_brakes[1] = p0s[dof];
                    profiles[dof].v_brakes[1] = v0s[dof];
                    profiles[dof].a_brakes[1] = a0s[dof];
                    std::tie(p0s[dof], v0s[dof], a0s[dof]) = Profile::integrate(profiles[dof].t_brakes[1], p0s[dof], v0s[dof], a0s[dof], profiles[dof].j_brakes[1]);
                }
            }

            Step1 step1 {p0s[dof], v0s[dof], a0s[dof], input.target_position[dof], input.target_velocity[dof], input.target_acceleration[dof], input.max_velocity[dof], input_min_velocity[dof], input.max_acceleration[dof], input.max_jerk[dof]};
            bool found_profile = step1.get_profile(profiles[dof]);
            if (!found_profile) {
                if constexpr (THROW_ERROR) {
                    throw std::runtime_error("[ruckig] error in step 1, dof: " + std::to_string(dof) + " input: " + input.to_string());
                }
                return Result::ErrorExecutionTimeCalculation;
            }

            blocks[dof] = step1.block;
            output.independent_min_durations[dof] = step1.block.t_min;
        }

        // std::cout << "step1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000.0 << std::endl;

        int limiting_dof; // The DoF that doesn't need step 2
        bool found_synchronization = synchronize(blocks, input.minimum_duration, tf, limiting_dof, profiles);
        if (!found_synchronization) {
            if constexpr (THROW_ERROR) {
                throw std::runtime_error("[ruckig] error in time synchronization: " + std::to_string(tf));
            }
            return Result::ErrorSynchronizationCalculation;
        }
        
        if (tf > 0.0) {
            for (size_t dof = 0; dof < DOFs; ++dof) {
                if (!input.enabled[dof] || dof == limiting_dof) {
                    continue;
                }

                double t_profile = tf - profiles[dof].t_brake.value_or(0.0);

                Step2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], input.target_position[dof], input.target_velocity[dof], input.target_acceleration[dof], input.max_velocity[dof], input_min_velocity[dof], input.max_acceleration[dof], input.max_jerk[dof]};
                bool found_time_synchronization = step2.get_profile(profiles[dof]);
                if (!found_time_synchronization) {
                    if constexpr (THROW_ERROR) {
                        throw std::runtime_error("[ruckig] error in step 2 in dof: " + std::to_string(dof) + " for t sync: " + std::to_string(tf) + " input: " + input.to_string());
                    }
                    return Result::ErrorSynchronizationCalculation;
                }
            }
        }

        // std::cout << "step2: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000.0 << std::endl;

        t = 0.0;
        output.duration = tf;
        output.new_calculation = true;
        return Result::Working;
    }

public:
    //! Just a shorter notation
    using Input = InputParameter<DOFs>;
    using Output = OutputParameter<DOFs>;

    //! Time step between updates (cycle time) in [s]
    const double delta_time;

    explicit Ruckig(double delta_time): delta_time(delta_time) { }

    bool validate_input(const InputParameter<DOFs>& input) {
        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (input.max_velocity[dof] <= std::numeric_limits<double>::min()) {
                std::cerr << "[ruckig] velocity limit needs to be positive." << std::endl;
                return false;
            }

            if (input.min_velocity && input.min_velocity.value()[dof] >= -std::numeric_limits<double>::min()) {
                std::cerr << "[ruckig] minimum velocity limit needs to be negative." << std::endl;
                return false;
            }

            if (input.max_acceleration[dof] <= std::numeric_limits<double>::min()) {
                std::cerr << "[ruckig] acceleration limit needs to be positive." << std::endl;
                return false;
            }

            if (input.max_jerk[dof] <= std::numeric_limits<double>::min()) {
                std::cerr << "[ruckig] jerk limit needs to be positive." << std::endl;
                return false;
            }

            if (std::isnan(input.target_position[dof])) {
                std::cerr << "[ruckig] target position is not a number." << std::endl;
                return false;
            }

            if (input.min_velocity) {
                if (input.target_velocity[dof] > input.max_velocity[dof] || input.target_velocity[dof] < input.min_velocity.value()[dof]) {
                    std::cerr << "[ruckig] target velocity exceeds velocity limit." << std::endl;
                    return false;
                }
            
            } else {
                if (std::abs(input.target_velocity[dof]) > input.max_velocity[dof]) {
                    std::cerr << "[ruckig] target velocity exceeds velocity limit." << std::endl;
                    return false;
                }
            }

            if (input.target_acceleration[dof] > input.max_acceleration[dof]) {
                std::cerr << "[ruckig] target acceleration exceeds acceleration limit." << std::endl;
                return false;
            }

            double max_target_acceleration = std::sqrt(2 * input.max_jerk[dof] * (input.max_velocity[dof] - std::abs(input.target_velocity[dof])));
            if (std::abs(input.target_acceleration[dof]) > max_target_acceleration) {
                std::cerr << "[ruckig] target acceleration exceeds maximal possible acceleration." << std::endl;
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
                std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(time - tf, current_input.target_position[dof], current_input.target_velocity[dof], current_input.target_acceleration[dof], 0);
            }
            return;
        }

        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (!current_input.enabled[dof]) {
                std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(time, current_input.current_position[dof], current_input.current_velocity[dof], current_input.current_acceleration[dof], 0);
            }

            auto& p = profiles[dof];

            double t_diff = time;
            if (p.t_brake.has_value()) {
                if (t_diff < p.t_brake.value()) {
                    size_t index = (t_diff < p.t_brakes[0]) ? 0 : 1;
                    if (index > 0) {
                        t_diff -= p.t_brakes[index - 1];
                    }

                    std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(t_diff, p.p_brakes[index], p.v_brakes[index], p.a_brakes[index], p.j_brakes[index]);
                    continue;
                } else {
                    t_diff -= p.t_brake.value();
                }
            }

            // Non-time synchronization
            if (t_diff >= p.t_sum[6]) {
                output.new_position[dof] = p.p[7];
                output.new_velocity[dof] = p.v[7];
                output.new_acceleration[dof] = p.a[7];
                continue;
            }

            auto index_ptr = std::upper_bound(p.t_sum.begin(), p.t_sum.end(), t_diff);
            size_t index = std::distance(p.t_sum.begin(), index_ptr);

            if (index > 0) {
                t_diff -= p.t_sum[index - 1];
            }

            std::tie(output.new_position[dof], output.new_velocity[dof], output.new_acceleration[dof]) = Profile::integrate(t_diff, p.p[index], p.v[index], p.a[index], p.j[index]);
        }
    }

    double get_time() const {
        return t;
    }
};

} // namespace ruckig
