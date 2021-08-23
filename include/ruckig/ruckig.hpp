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

constexpr static size_t DynamicDOFs {0};

//! Main class for the Ruckig algorithm.
template<size_t DOFs = 0, bool throw_error = false, bool return_error_at_maximal_duration = true>
class Ruckig {
    //! Current input, only for comparison for recalculation
    InputParameter<DOFs> current_input;

public:
    size_t degrees_of_freedom;

    //! Time step between updates (cycle time) in [s]
    const double delta_time;

    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    explicit Ruckig(): degrees_of_freedom(DOFs), delta_time(-1.0) {
    }

    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    explicit Ruckig(double delta_time): degrees_of_freedom(DOFs), delta_time(delta_time) {
    }


    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    explicit Ruckig(size_t dofs): degrees_of_freedom(dofs), delta_time(-1.0), current_input(InputParameter<0>(dofs)) {
    }

    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    explicit Ruckig(size_t dofs, double delta_time): degrees_of_freedom(dofs), delta_time(delta_time), current_input(InputParameter<0>(dofs)) {
    }


    //! Validate the input for the trajectory calculation
    bool validate_input(const InputParameter<DOFs>& input) const {
        for (size_t dof = 0; dof < degrees_of_freedom; ++dof) {
            if (input.control_interface == ControlInterface::Position && std::isnan(input.current_position[dof])) {
                return false;
            }

            if (input.control_interface == ControlInterface::Position && std::isnan(input.max_velocity[dof])) {
                return false;
            }

            if (input.control_interface == ControlInterface::Position && input.max_velocity[dof] <= std::numeric_limits<double>::min()) {
                return false;
            }

            if (input.min_velocity && input.min_velocity.value()[dof] >= -std::numeric_limits<double>::min()) {
                return false;
            }

            if (std::isnan(input.max_acceleration[dof])) {
                return false;
            }

            if (input.max_acceleration[dof] <= std::numeric_limits<double>::min()) {
                return false;
            }

            if (input.min_acceleration && input.min_acceleration.value()[dof] >= -std::numeric_limits<double>::min()) {
                return false;
            }

            if (std::isnan(input.max_jerk[dof])) {
                return false;
            }

            if (input.max_jerk[dof] <= std::numeric_limits<double>::min()) {
                return false;
            }

            if (input.control_interface == ControlInterface::Position && std::isnan(input.target_position[dof])) {
                return false;
            }

            if (input.control_interface == ControlInterface::Position) {
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

            // Target acceleration needs to be accessible from "above" and "below"
            if (input.control_interface == ControlInterface::Position) {
                const double min_velocity = input.min_velocity ? input.min_velocity.value()[dof] : -input.max_velocity[dof];
                const double v_diff = std::min(std::abs(input.max_velocity[dof] - input.target_velocity[dof]), std::abs(min_velocity - input.target_velocity[dof]));
                const double max_target_acceleration = std::sqrt(2 * input.max_jerk[dof] * v_diff);
                if (std::abs(input.target_acceleration[dof]) > max_target_acceleration) {
                    return false;
                }
            }
        }

        // Check for intermediate waypoints here
        if (!input.intermediate_positions.empty()) {
            return false;
        }

        return true;
    }

    //! Calculate a new trajectory for the given input
    Result calculate(const InputParameter<DOFs>& input, Trajectory<DOFs>& trajectory) {
        if (!validate_input(input)) {
            return Result::ErrorInvalidInput;
        }

        return trajectory.template calculate<throw_error, return_error_at_maximal_duration>(input, delta_time);
    }

    //! Get the next output state (with step delta_time) along the calculated trajectory for the given input
    Result update(const InputParameter<DOFs>& input, OutputParameter<DOFs>& output) {
        const auto start = std::chrono::high_resolution_clock::now();

        if constexpr (DOFs == 0 && throw_error) {
            if (degrees_of_freedom != input.degrees_of_freedom || degrees_of_freedom != output.degrees_of_freedom) {
                throw std::runtime_error("[ruckig] mismatch in degrees of freedom (vector size).");
            }
        }

        output.new_calculation = false;

        if (input != current_input) {
            Result result = calculate(input, output.trajectory);
            if (result != Result::Working) {
                return result;
            }

            current_input = input;
            output.time = 0.0;
            output.new_calculation = true;
        }

        output.time += delta_time;
        output.trajectory.at_time(output.time, output.new_position, output.new_velocity, output.new_acceleration);

        const auto stop = std::chrono::high_resolution_clock::now();
        output.calculation_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000.0;

        current_input.current_position = output.new_position;
        current_input.current_velocity = output.new_velocity;
        current_input.current_acceleration = output.new_acceleration;

        if (output.time > output.trajectory.get_duration()) {
            return Result::Finished;
        }

        return Result::Working;
    }
};

} // namespace ruckig
