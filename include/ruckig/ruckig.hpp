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

#include <ruckig/input_parameter.hpp>
#include <ruckig/output_parameter.hpp>
#include <ruckig/trajectory.hpp>


namespace ruckig {

//! Main class for the Ruckig algorithm.
template<size_t DOFs, bool throw_error = false, bool return_error_at_maximal_duration = true>
class Ruckig {
    //! Current input, only for comparison for recalculation
    InputParameter<DOFs> current_input;

    Result calculate(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        current_input = input;

        if (!validate_input(input)) {
            return Result::ErrorInvalidInput;
        }

        Result result = output.trajectory.template calculate<throw_error, return_error_at_maximal_duration>(input, delta_time);
        if (result != Result::Working) {
            return result;
        }

        output.time = 0.0;
        output.new_calculation = true;
        return Result::Working;
    }

public:
    // Just a shorter notation
    using Input = InputParameter<DOFs>;
    using Output = OutputParameter<DOFs>;
    static constexpr size_t degrees_of_freedom {DOFs};

    //! Time step between updates (cycle time) in [s]
    double delta_time;

    explicit Ruckig() { }
    explicit Ruckig(double delta_time): delta_time(delta_time) { }

    bool validate_input(const InputParameter<DOFs>& input) const {
        for (size_t dof = 0; dof < DOFs; ++dof) {
            if (input.interface == Interface::Position && input.max_velocity[dof] <= std::numeric_limits<double>::min()) {
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

            if (input.interface == Interface::Position && std::isnan(input.target_position[dof])) {
                return false;
            }

            if (input.interface == Interface::Position) {
                if (input.min_velocity) {
                    if (input.target_velocity[dof] > input.max_velocity[dof] || input.target_velocity[dof] < input.min_velocity.value()[dof]) {
                        return false;
                    }

                } else {
                    if (std::abs(input.target_velocity[dof]) > input.max_velocity[dof]) {
                        return false;
                    }
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

            if (input.interface == Interface::Position) {
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

        current_input.current_position = output.new_position;
        current_input.current_velocity = output.new_velocity;
        current_input.current_acceleration = output.new_acceleration;

        if (output.time > output.trajectory.duration) {
            return Result::Finished;
        }

        return Result::Working;
    }
};

} // namespace ruckig
