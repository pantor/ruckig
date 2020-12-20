# pragma once

#include <Eigen/Core>

#include <ruckig/parameter.hpp>


namespace ruckig {

/**
 * Adapted from: Wisama Khalil and Etienne Dombre. 2002. Modeling, Identification and Control of Robots (Kogan Page Science Paper edition).
 */
template<size_t DOFs>
class Smoothie {
    using Vector = Eigen::Matrix<double, DOFs, 1, Eigen::ColMajor>;

    static constexpr double q_delta_motion_finished {1e-6};

    InputParameter<DOFs> current_input;
    double time;

    Vector q_initial, q_delta;
    Vector dq_max_sync_, q_1_;
    Vector t_1_sync, t_2_sync, t_f_sync;
    Vector dq_max_, ddq_max_initial, ddq_max_target;

    void calculateSynchronizedValues() {
        Vector dq_max_reach(dq_max_);
        Vector t_f = Vector::Zero();
        Vector delta_t_2 = Vector::Zero();
        Vector t_1 = Vector::Zero();
        Vector delta_t_2_sync = Vector::Zero();
        Vector sign_delta_q = q_delta.cwiseSign();

        time = 0.0;

        for (size_t i = 0; i < DOFs; i++) {
            if (std::abs(q_delta[i]) > q_delta_motion_finished) {
                if (std::abs(q_delta[i]) < (3.0 / 4.0 * (std::pow(dq_max_[i], 2.0) / ddq_max_initial[i]) + 3.0 / 4.0 * (std::pow(dq_max_[i], 2.0) / ddq_max_target[i]))) {
                    dq_max_reach[i] = std::sqrt(4.0 / 3.0 * q_delta[i] * sign_delta_q[i] * (ddq_max_initial[i] * ddq_max_target[i]) / (ddq_max_initial[i] + ddq_max_target[i]));
                }
                t_1[i] = 1.5 * dq_max_reach[i] / ddq_max_initial[i];
                delta_t_2[i] = 1.5 * dq_max_reach[i] / ddq_max_target[i];
                t_f[i] = t_1[i] / 2.0 + delta_t_2[i] / 2.0 + std::abs(q_delta[i]) / dq_max_reach[i];
            }
        }

        double max_t_f = t_f.maxCoeff();
        for (size_t i = 0; i < DOFs; i++) {
            if (std::abs(q_delta[i]) > q_delta_motion_finished) {
                double a = 1.5 / 2.0 * (ddq_max_target[i] + ddq_max_initial[i]);
                double b = -1.0 * max_t_f * ddq_max_target[i] * ddq_max_initial[i];
                double c = std::abs(q_delta[i]) * ddq_max_target[i] * ddq_max_initial[i];
                double delta = b * b - 4.0 * a * c;
                if (delta < 0.0) {
                    delta = 0.0;
                }
                dq_max_sync_[i] = (-1.0 * b - std::sqrt(delta)) / (2.0 * a);
                t_1_sync[i] = 1.5 * dq_max_sync_[i] / ddq_max_initial[i];
                delta_t_2_sync[i] = 1.5 * dq_max_sync_[i] / ddq_max_target[i];
                t_f_sync[i] = (t_1_sync)[i] / 2.0 + delta_t_2_sync[i] / 2.0 + std::abs(q_delta[i] / dq_max_sync_[i]);
                t_2_sync[i] = (t_f_sync)[i] - delta_t_2_sync[i];
                q_1_[i] = (dq_max_sync_)[i] * sign_delta_q[i] * (0.5 * (t_1_sync)[i]);
            }
        }
    }

    bool calculateDesiredValues(double t, Vector& q_delta_d) const {
        Vector sign_delta_q = q_delta.cwiseSign();
        Vector t_d = t_2_sync - t_1_sync;
        Vector delta_t_2_sync = t_f_sync - t_2_sync;
        std::array<bool, DOFs> joint_motion_finished {};

        for (size_t i = 0; i < DOFs; i++) {
            if (std::abs(q_delta[i]) < q_delta_motion_finished) {
                q_delta_d[i] = 0;
                joint_motion_finished[i] = true;
            } else {
                if (t < t_1_sync[i]) {
                    q_delta_d[i] = -1.0 / std::pow(t_1_sync[i], 3.0) * dq_max_sync_[i] * sign_delta_q[i] * (0.5 * t - t_1_sync[i]) * std::pow(t, 3.0);
                } else if (t >= t_1_sync[i] && t < t_2_sync[i]) {
                    q_delta_d[i] = q_1_[i] + (t - t_1_sync[i]) * dq_max_sync_[i] * sign_delta_q[i];
                } else if (t >= t_2_sync[i] && t < t_f_sync[i]) {
                    q_delta_d[i] = q_delta[i] + 0.5 * (1.0 / std::pow(delta_t_2_sync[i], 3.0) * (t - t_1_sync[i] - 2.0 * delta_t_2_sync[i] - t_d[i]) * std::pow((t - t_1_sync[i] - t_d[i]), 3.0) + (2.0 * t - 2.0 * t_1_sync[i] - delta_t_2_sync[i] - 2.0 * t_d[i])) * dq_max_sync_[i] * sign_delta_q[i];
                } else {
                    q_delta_d[i] = q_delta[i];
                    joint_motion_finished[i] = true;
                }
            }
        }
        return std::all_of(joint_motion_finished.cbegin(), joint_motion_finished.cend(), [](bool x) { return x; });
    }

public:
    double delta_time;

    explicit Smoothie(double delta_time): delta_time(delta_time) { }

    Result update(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        time += delta_time;
        output.new_calculation = false;

        if (input != current_input) {
            current_input = input;

            if ((input.max_velocity.array() <= 0.0).any() || (input.max_acceleration.array() <= 0.0).any()) {
                return Result::Error;
            }
            if ((input.target_velocity.array() != 0.0).any() || (input.target_acceleration.array() != 0.0).any()) {
                return Result::Error;
            }

            dq_max_ = input.max_velocity;
            ddq_max_initial = input.max_acceleration;
            ddq_max_target = input.max_acceleration;

            q_initial = input.current_position;
            q_delta = input.target_position - q_initial;
            calculateSynchronizedValues();

            output.duration = t_f_sync.maxCoeff();
            output.new_calculation = true;
        }

        Vector q_delta_d;
        bool motion_finished = calculateDesiredValues(time, q_delta_d);

        output.new_position = q_initial + q_delta_d;
        output.new_velocity = Vector::Zero();
        output.new_acceleration = Vector::Zero();

        current_input.current_position = output.new_position;
        current_input.current_velocity = output.new_velocity;
        current_input.current_acceleration = output.new_acceleration;

        if (motion_finished) {
            output.new_position = input.target_position;
            output.new_velocity = input.target_velocity;
            output.new_acceleration = input.target_acceleration;

            return Result::Finished;
        }
        return Result::Working;
    }
};

} // namespace ruckig
