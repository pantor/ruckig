#pragma once

#include <ruckig/input_parameter.hpp>
#include <ruckig/output_parameter.hpp>


namespace ruckig {

//! Alternative OTG algorithm (for comparison)
template<size_t DOFs>
class Quintic {
    std::array<double, DOFs> a, b, c, d, e, f;
    double t, tf;
    InputParameter<DOFs> current_input;

    bool validate_input(const InputParameter<DOFs>& input) {
        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (input.max_velocity[dof] <= 0.0) {
                std::cerr << "[quintic] velocity limit needs to be positive." << std::endl;
                return false;
            }

            if (input.max_acceleration[dof] <= 0.0) {
                std::cerr << "[quintic] acceleration limit needs to be positive." << std::endl;
                return false;
            }

            if (input.max_jerk[dof] <= 0.0) {
                std::cerr << "[quintic] jerk limit needs to be positive." << std::endl;
                return false;
            }

            if (std::isnan(input.target_position[dof])) {
                std::cerr << "[quintic] target position is not a number." << std::endl;
                return false;
            }

            if (input.target_velocity[dof] > input.max_velocity[dof]) {
                std::cerr << "[quintic] target velocity exceeds velocity limit." << std::endl;
                return false;
            }

            if (input.target_acceleration[dof] > input.max_acceleration[dof]) {
                std::cerr << "[quintic] target acceleration exceeds acceleration limit." << std::endl;
                return false;
            }
        }

        return true;
    }

    bool calculate(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        current_input = input;

        if (!validate_input(input)) {
            return false;
        }

        const auto& x0 = input.current_position;
        const auto& v0 = input.current_velocity;
        const auto& a0 = input.current_acceleration;
        const auto& xf = input.target_position;
        const auto& vf = input.target_velocity;
        const auto& af = input.target_acceleration;
        const auto& v_max = input.max_velocity;
        const auto& a_max = input.max_acceleration;
        const auto& j_max = input.max_jerk;

        // Approximations for v0 == 0, vf == 0, a0 == 0, af == 0
        std::array<double, DOFs> v_max_tfs, a_max_tfs, j_max_tfs;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            v_max_tfs[dof] = (15 * std::abs(x0[dof] - xf[dof])) / (8 * v_max[dof]);
            a_max_tfs[dof] = (std::sqrt(10) * std::pow(std::pow(x0[dof], 2) - 2 * x0[dof] * xf[dof] + std::pow(xf[dof], 2), 1./4)) / (std::pow(3, 1./4) * std::sqrt(a_max[dof]));
            j_max_tfs[dof] = std::cbrt((60 * std::abs(x0[dof] - xf[dof])) / j_max[dof]); // Also solvable for v0 != 0
        }

        double v_max_tfs_max = *std::max_element(v_max_tfs.cbegin(), v_max_tfs.cend());
        double a_max_tfs_max = *std::max_element(a_max_tfs.cbegin(), a_max_tfs.cend());
        double j_max_tfs_max = *std::max_element(j_max_tfs.cbegin(), j_max_tfs.cend());

        tf = std::max<double>({v_max_tfs_max, a_max_tfs_max, j_max_tfs_max});
        if (input.minimum_duration) {
            tf = std::max<double>({tf, input.minimum_duration.value()});
        }

        for (size_t dof = 0; dof < DOFs; ++dof) {
            a[dof] = -((a0[dof] - af[dof]) * std::pow(tf, 2) + 6 * tf * (v0[dof] + vf[dof]) + 12 * (x0[dof] - xf[dof])) / (2 * std::pow(tf, 5));
            b[dof] = -((2 * af[dof] - 3 * a0[dof]) * std::pow(tf, 2) - 16 * tf * v0[dof] - 14 * tf * vf[dof] - 30 * (x0[dof] - xf[dof])) / (2 * std::pow(tf, 4));
            c[dof] = -((3 * a0[dof] - af[dof]) * std::pow(tf, 2) + 12 * tf * v0[dof] + 8 * tf * vf[dof] + 20 * (x0[dof] - xf[dof])) / (2 * std::pow(tf, 3));
            d[dof] = a0[dof] / 2;
            e[dof] = v0[dof];
            f[dof] = x0[dof];
        }

        t = 0.0;
        output.trajectory.duration = tf;
        output.new_calculation = true;
        return true;
    }

public:
    double delta_time;
    static constexpr size_t degrees_of_freedom {DOFs};

    explicit Quintic(double delta_time): delta_time(delta_time) { }

    Result update(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        t += delta_time;
        output.new_calculation = false;

        if (input != current_input && !calculate(input, output)) {
            return Result::Error;
        }

        if (t >= tf) {
            // Keep velocity
            for (size_t dof = 0; dof < DOFs; ++dof) {
                output.new_position[dof] = input.target_position[dof] + (t - tf) * input.target_velocity[dof];
                output.new_velocity[dof] = input.target_velocity[dof];
                output.new_acceleration[dof] = 0.0;
            }

            return Result::Finished;
        }

        for (size_t dof = 0; dof < DOFs; ++dof) {
            output.new_position[dof] = f[dof] + t * (e[dof] + t * (d[dof] + t * (c[dof] + t * (b[dof] + a[dof] * t))));
            output.new_velocity[dof] = e[dof] + t * (2 * d[dof] + t * (3 * c[dof] + t * (4 * b[dof] + 5 * a[dof] * t)));
            output.new_acceleration[dof] = 2 * d[dof] + t * (6 * c[dof] + t * (12 * b[dof] + t * (20 * a[dof])));
        }

        current_input.current_position = output.new_position;
        current_input.current_velocity = output.new_velocity;
        current_input.current_acceleration = output.new_acceleration;
        return Result::Working;
    }
};

} // namespace ruckig
