#pragma once

#include <array>
#include <chrono>
#include <iostream>
#include <optional>

#include <ruckig/parameter.hpp>


namespace ruckig {

struct Profile {
    //! Profile names indicate which limits are reached.
    enum class Type {
        UP_ACC0_ACC1_VEL, UP_VEL, UP_ACC0, UP_ACC1, UP_ACC0_ACC1, UP_ACC0_VEL, UP_ACC1_VEL, UP_NONE,
        DOWN_ACC0_ACC1_VEL, DOWN_VEL, DOWN_ACC0, DOWN_ACC1, DOWN_ACC0_ACC1, DOWN_ACC0_VEL, DOWN_ACC1_VEL, DOWN_NONE
    };

    Type type;
    std::array<double, 7> t, t_sum, j;
    std::array<double, 8> a, v, p;

    //! Total time of the braking segments
    std::optional<double> t_brake;

    //! Allow up to two segments of braking before the "correct" profile starts
    std::array<double, 2> t_brakes, j_brakes, a_brakes, v_brakes, p_brakes;

    void set(double p0, double v0, double a0, const std::array<double, 7>& j);
    bool check(double pf, double vf, double af, double vMax, double aMax) const;
    bool check(double tf, double pf, double vf, double af, double vMax, double aMax) const;

    //! Integrate with constant jerk for duration t. Returns new position, new velocity, and new acceleration.
    static std::tuple<double, double, double> integrate(double t, double p0, double v0, double a0, double j);
};


struct RuckigStep1 {
    double p0, v0, a0;
    double pf, vf, af;

    explicit RuckigStep1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax);

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

    bool get_profile(Profile& profile, double vMax, double aMax, double jMax);

    static void get_brake_trajectory(double v0, double a0, double vMax, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);
};


struct RuckigStep2 {
    double tf;
    double p0, v0, a0;
    double pf, vf, af;

    explicit RuckigStep2(double tf, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax);

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

    bool get_profile(Profile& profile, double vMax, double aMax, double jMax);
};


template<size_t DOFs>
class Ruckig {
    InputParameter<DOFs> current_input;

    double t, tf;
    std::array<Profile, DOFs> profiles;

    bool calculate(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        current_input = input;

        // Check input
        if ((input.max_velocity.array() <= 0.0).any() || (input.max_acceleration.array() <= 0.0).any() || (input.max_jerk.array() <= 0.0).any()) {
            return false;
        }

        if ((input.target_velocity.array().abs() > input.max_velocity.array()).any()) {
            std::cerr << "Target velocity exceeds maximal velocity." << std::endl;
            return false;
        }

        if ((input.target_acceleration.array() != 0.0).any()) {
            std::cerr << "Ruckig does not support a target acceleration." << std::endl;
            return false;
        }

        if (input.minimum_duration.has_value()) {
            std::cerr << "Ruckig does not support a minimum duration." << std::endl;
            return false;
        }

        // Calculate brakes (if input exceeds or will exceed limits)
        for (size_t dof = 0; dof < DOFs; dof += 1) {
            if (!input.enabled[dof]) {
                continue;
            }

            RuckigStep1::get_brake_trajectory(input.current_velocity[dof], input.current_acceleration[dof], input.max_velocity[dof], input.max_acceleration[dof], input.max_jerk[dof], profiles[dof].t_brakes, profiles[dof].j_brakes);
            profiles[dof].t_brake = profiles[dof].t_brakes[0] + profiles[dof].t_brakes[1];
        }

        std::array<double, DOFs> tfs; // Profile duration
        std::array<double, DOFs> p0s, v0s, a0s; // Starting point of profiles without brake trajectory
        for (size_t dof = 0; dof < DOFs; dof += 1) {
            if (!input.enabled[dof]) {
                tfs[dof] = 0.0;
                continue;
            }

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

            RuckigStep1 step1 {p0s[dof], v0s[dof], a0s[dof], input.target_position[dof], input.target_velocity[dof], input.target_acceleration[dof], input.max_velocity[dof], input.max_acceleration[dof], input.max_jerk[dof]};
            bool found_profile = step1.get_profile(profiles[dof], input.max_velocity[dof], input.max_acceleration[dof], input.max_jerk[dof]);
            if (!found_profile) {
                throw std::runtime_error("[ruckig] error in step 1: " + input.to_string(dof) + " all: " + input.to_string());
            }
            tfs[dof] = profiles[dof].t_sum[6] + profiles[dof].t_brake.value_or(0.0);
        }

        auto tf_max_pointer = std::max_element(tfs.begin(), tfs.end());
        size_t limiting_dof = std::distance(tfs.begin(), tf_max_pointer);
        tf = *tf_max_pointer;

        if (tf > 0.0) {
            for (size_t dof = 0; dof < DOFs; dof += 1) {
                if (!input.enabled[dof] || dof == limiting_dof) {
                    continue;
                }

                double t_profile = tf - profiles[dof].t_brake.value_or(0.0);

                RuckigStep2 step2 {t_profile, p0s[dof], v0s[dof], a0s[dof], input.target_position[dof], input.target_velocity[dof], input.target_acceleration[dof], input.max_velocity[dof], input.max_acceleration[dof], input.max_jerk[dof]};
                bool found_time_synchronization = step2.get_profile(profiles[dof], input.max_velocity[dof], input.max_acceleration[dof], input.max_jerk[dof]);
                if (!found_time_synchronization) {
                    throw std::runtime_error("[ruckig] error in step 2: " + input.to_string(dof) + " all: " + input.to_string());
                }
            }
        }

        t = 0.0;
        output.duration = tf;
        output.new_calculation = true;
        return true;
    }

public:
    //! Time step between updates (cycle time) in [s]
    const double delta_time;

    explicit Ruckig(double delta_time): delta_time(delta_time) { }

    Result update(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        t += delta_time;
        output.new_calculation = false;

        auto start = std::chrono::high_resolution_clock::now();

        if (input != current_input && !calculate(input, output)) {
            return Result::Error;
        }

        atTime(t, output);

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

    void atTime(double time, OutputParameter<DOFs>& output) {
        if (time + delta_time > tf) {
            // Keep velocity
            output.new_position = current_input.target_position + current_input.target_velocity * (time - tf);
            output.new_velocity = current_input.target_velocity;
            output.new_acceleration = current_input.target_acceleration;
            return;
        }

        for (size_t dof = 0; dof < DOFs; dof += 1) {
            if (!current_input.enabled[dof]) {
                output.new_acceleration[dof] = current_input.current_acceleration[dof];
                output.new_velocity[dof] = current_input.current_velocity[dof];
                output.new_position[dof] = current_input.current_position[dof];
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
};

} // namespace ruckig
