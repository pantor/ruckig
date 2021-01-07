# pragma once

#include <chrono>

#include <ruckig/parameter.hpp>


namespace ruckig {

/**
 * Adapted from: Wisama Khalil and Etienne Dombre. 2002. Modeling, Identification and Control of Robots (Kogan Page Science Paper edition).
 */
template<size_t DOFs>
class Smoothie {
    using Vector = std::array<double, DOFs>;

    static constexpr double q_delta_motion_finished {1e-6};

    InputParameter<DOFs> current_input;
    double time;

    Vector q_initial, q_delta;
    Vector dq_max_sync_, q_1_;
    Vector t_1_sync, t_2_sync, t_f_sync;
    Vector dq_max_, ddq_max_initial, ddq_max_target;

    void calculateSynchronizedValues(const InputParameter<DOFs>& input) {
        Vector dq_max_reach;
        Vector t_f {};
        Vector delta_t_2 {};
        Vector t_1 {};
        Vector delta_t_2_sync {};
        Vector sign_delta_q;
        
        time = 0.0;

        for (size_t dof = 0; dof < DOFs; dof++) {
            dq_max_reach[dof] = dq_max_[dof];
            sign_delta_q[dof] = std::signbit(q_delta[dof]) ? -1 : 1;

            if (std::abs(q_delta[dof]) > q_delta_motion_finished) {
                if (std::abs(q_delta[dof]) < (3.0 / 4 * (std::pow(dq_max_[dof], 2) / ddq_max_initial[dof]) + 3.0 / 4 * (std::pow(dq_max_[dof], 2) / ddq_max_target[dof]))) {
                    dq_max_reach[dof] = std::sqrt(4.0 / 3 * q_delta[dof] * sign_delta_q[dof] * (ddq_max_initial[dof] * ddq_max_target[dof]) / (ddq_max_initial[dof] + ddq_max_target[dof]));
                }
                t_1[dof] = 1.5 * dq_max_reach[dof] / ddq_max_initial[dof];
                delta_t_2[dof] = 1.5 * dq_max_reach[dof] / ddq_max_target[dof];
                t_f[dof] = t_1[dof] / 2 + delta_t_2[dof] / 2 + std::abs(q_delta[dof]) / dq_max_reach[dof];
            }
        }

        double max_t_f = *std::max_element(t_f.cbegin(), t_f.cend()); 
        if (input.minimum_duration.has_value()) {
            max_t_f = std::max<double>({max_t_f, input.minimum_duration.value()});
        }

        for (size_t dof = 0; dof < DOFs; dof++) {
            if (std::abs(q_delta[dof]) > q_delta_motion_finished) {
                double a = 1.5 / 2 * (ddq_max_target[dof] + ddq_max_initial[dof]);
                double b = -1.0 * max_t_f * ddq_max_target[dof] * ddq_max_initial[dof];
                double c = std::abs(q_delta[dof]) * ddq_max_target[dof] * ddq_max_initial[dof];
                double delta = b * b - 4 * a * c;
                if (delta < 0.0) {
                    delta = 0.0;
                }
                dq_max_sync_[dof] = (-1.0 * b - std::sqrt(delta)) / (2 * a);
                t_1_sync[dof] = 1.5 * dq_max_sync_[dof] / ddq_max_initial[dof];
                delta_t_2_sync[dof] = 1.5 * dq_max_sync_[dof] / ddq_max_target[dof];
                t_f_sync[dof] = (t_1_sync)[dof] / 2 + delta_t_2_sync[dof] / 2 + std::abs(q_delta[dof] / dq_max_sync_[dof]);
                t_2_sync[dof] = (t_f_sync)[dof] - delta_t_2_sync[dof];
                q_1_[dof] = (dq_max_sync_)[dof] * sign_delta_q[dof] * (0.5 * (t_1_sync)[dof]);
            }
        }
    }

    bool calculateDesiredValues(double t, Vector& q_delta_d) const {
        Vector sign_delta_q;
        Vector t_d;
        Vector delta_t_2_sync;
        std::array<bool, DOFs> joint_motion_finished {};

        for (size_t dof = 0; dof < DOFs; dof++) {
            t_d[dof] = t_2_sync[dof] - t_1_sync[dof];
            delta_t_2_sync[dof] = t_f_sync[dof] - t_2_sync[dof];
            sign_delta_q[dof] = std::signbit(q_delta[dof]) ? -1 : 1;

            if (std::abs(q_delta[dof]) < q_delta_motion_finished) {
                q_delta_d[dof] = 0;
                joint_motion_finished[dof] = true;
            } else {
                if (t < t_1_sync[dof]) {
                    q_delta_d[dof] = -1.0 / std::pow(t_1_sync[dof], 3) * dq_max_sync_[dof] * sign_delta_q[dof] * (0.5 * t - t_1_sync[dof]) * std::pow(t, 3);
                } else if (t >= t_1_sync[dof] && t < t_2_sync[dof]) {
                    q_delta_d[dof] = q_1_[dof] + (t - t_1_sync[dof]) * dq_max_sync_[dof] * sign_delta_q[dof];
                } else if (t >= t_2_sync[dof] && t < t_f_sync[dof]) {
                    q_delta_d[dof] = q_delta[dof] + 0.5 * (1.0 / std::pow(delta_t_2_sync[dof], 3) * (t - t_1_sync[dof] - 2.0 * delta_t_2_sync[dof] - t_d[dof]) * std::pow((t - t_1_sync[dof] - t_d[dof]), 3) + (2 * t - 2 * t_1_sync[dof] - delta_t_2_sync[dof] - 2 * t_d[dof])) * dq_max_sync_[dof] * sign_delta_q[dof];
                } else {
                    q_delta_d[dof] = q_delta[dof];
                    joint_motion_finished[dof] = true;
                }
            }
        }
        return std::all_of(joint_motion_finished.cbegin(), joint_motion_finished.cend(), [](bool x) { return x; });
    }

public:
    double delta_time;

    explicit Smoothie(double delta_time): delta_time(delta_time) { }

    bool validate_input(const InputParameter<DOFs>& input) {
        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (input.max_velocity[dof] <= 0.0) {
                std::cerr << "Velocity limit needs to be positive." << std::endl;
                return false;
            }

            if (input.max_acceleration[dof] <= 0.0) {
                std::cerr << "Acceleration limit needs to be positive." << std::endl;
                return false;
            }

            if (input.target_velocity[dof] != 0.0) {
                std::cerr << "Target velocity exceeds velocity limit." << std::endl;
                return false;
            }

            if (input.target_acceleration[dof] != 0.0) {
                std::cerr << "Target acceleration exceeds acceleration limit." << std::endl;
                return false;
            }
        }

        return true;
    }

    Result update(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        time += delta_time;
        output.new_calculation = false;

        if (input != current_input) {
            auto start = std::chrono::high_resolution_clock::now();

            current_input = input;

            if (!validate_input(input)) {
                return Result::Error;
            }

            std::copy_n(input.max_velocity.data(), DOFs, dq_max_.begin());
            std::copy_n(input.max_acceleration.data(), DOFs, ddq_max_initial.begin());
            std::copy_n(input.max_acceleration.data(), DOFs, ddq_max_target.begin());
            std::copy_n(input.current_position.data(), DOFs, q_initial.begin());

            for (size_t dof = 0; dof < DOFs; ++dof) {
                q_delta[dof] = input.target_position[dof] - q_initial[dof];
            }
            calculateSynchronizedValues(input);

            output.duration = *std::max_element(t_f_sync.cbegin(), t_f_sync.cend());
            output.new_calculation = true;

            auto stop = std::chrono::high_resolution_clock::now();
            output.calculation_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000.0;
        }

        Vector q_delta_d;
        bool motion_finished = calculateDesiredValues(time, q_delta_d);

        for (size_t dof = 0; dof < DOFs; ++dof) {
            output.new_position[dof] = q_initial[dof] + q_delta_d[dof];
            output.new_velocity[dof] = 0.0;
            output.new_acceleration[dof] = 0.0;
        }

        current_input.current_position = output.new_position;
        current_input.current_velocity = output.new_velocity;
        current_input.current_acceleration = output.new_acceleration;

        if (motion_finished) {
            // Keep position
            for (size_t dof = 0; dof < DOFs; ++dof) {
                output.new_position[dof] = input.target_position[dof];
                output.new_velocity[dof] = 0.0;
                output.new_acceleration[dof] = 0.0;
            }

            return Result::Finished;
        }
        return Result::Working;
    }
};

} // namespace ruckig
