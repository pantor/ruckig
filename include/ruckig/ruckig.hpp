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

#include <ruckig/calculator.hpp>
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

    //! Calculator for new trajectories
    Calculator<DOFs> calculator;

    //! Max number of intermediate waypoints
    const size_t max_number_of_waypoints;

    inline static double v_at_a_zero(double v0, double a0, double j) {
        return v0 + (a0 * a0)/(2 * j);
    }

public:
    size_t degrees_of_freedom;

    //! Time step between updates (cycle time) in [s]
    const double delta_time;

    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    explicit Ruckig(): degrees_of_freedom(DOFs), delta_time(-1.0), max_number_of_waypoints(0) {
    }

    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    explicit Ruckig(double delta_time): degrees_of_freedom(DOFs), delta_time(delta_time), max_number_of_waypoints(0) {
    }

#if defined WITH_ONLINE_CLIENT
    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    explicit Ruckig(double delta_time, size_t max_number_of_waypoints): degrees_of_freedom(DOFs), delta_time(delta_time), calculator(Calculator<DOFs>(max_number_of_waypoints)), max_number_of_waypoints(max_number_of_waypoints), current_input(InputParameter<DOFs>(max_number_of_waypoints)) {
    }
#endif

    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    explicit Ruckig(size_t dofs): degrees_of_freedom(dofs), delta_time(-1.0), calculator(Calculator<0>(dofs)), max_number_of_waypoints(0), current_input(InputParameter<0>(dofs)) {
    }

    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    explicit Ruckig(size_t dofs, double delta_time): degrees_of_freedom(dofs), delta_time(delta_time), calculator(Calculator<0>(dofs)), max_number_of_waypoints(0), current_input(InputParameter<0>(dofs)) {
    }

#if defined WITH_ONLINE_CLIENT
    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    explicit Ruckig(size_t dofs, double delta_time, size_t max_number_of_waypoints): degrees_of_freedom(dofs), delta_time(delta_time), calculator(Calculator<0>(dofs, max_number_of_waypoints)), max_number_of_waypoints(max_number_of_waypoints), current_input(InputParameter<0>(dofs, max_number_of_waypoints)) {
    }
#endif

    //! Validate the input for trajectory calculation and kinematic limits
    bool validate_input(const InputParameter<DOFs>& input, bool check_current_state_within_limits=false, bool check_target_state_within_limits=true) const {
        for (size_t dof = 0; dof < degrees_of_freedom; ++dof) {
            const double jMax = input.max_jerk[dof];
            if (std::isnan(jMax) || jMax <= std::numeric_limits<double>::min()) {
                return false;
            }

            const double aMax = input.max_acceleration[dof];
            if (std::isnan(aMax) || aMax <= std::numeric_limits<double>::min()) {
                return false;
            }

            const double aMin = input.min_acceleration ? input.min_acceleration.value()[dof] : -input.max_acceleration[dof];
            if (std::isnan(aMin) || aMin >= -std::numeric_limits<double>::min()) {
                return false;
            }

            const double a0 = input.current_acceleration[dof];
            const double af = input.target_acceleration[dof];
            if (std::isnan(a0) || std::isnan(af)) {
                return false;
            }

            if (check_current_state_within_limits && (a0 > aMax || a0 < aMin)) {
                return false;
            }
            if (check_target_state_within_limits && (af > aMax || af < aMin)) {
                return false;  
            }

            const double v0 = input.current_velocity[dof];
            const double vf = input.target_velocity[dof];
            if (std::isnan(v0) || std::isnan(vf)) {
                return false;
            }

            auto control_interface = input.per_dof_control_interface ? input.per_dof_control_interface.value()[dof] : input.control_interface;
            if (control_interface == ControlInterface::Position) {
                if (std::isnan(input.current_position[dof]) || std::isnan(input.target_position[dof])) {
                    return false;
                }

                const double vMax = input.max_velocity[dof];
                if (std::isnan(vMax) || vMax <= std::numeric_limits<double>::min()) {
                    return false;
                }

                const double vMin = input.min_velocity ? input.min_velocity.value()[dof] : -input.max_velocity[dof];
                if (std::isnan(vMin) || vMin >= -std::numeric_limits<double>::min()) {
                    return false;
                }

                if (check_current_state_within_limits && (v0 > vMax || v0 < vMin)) {
                    return false;
                }
                if (check_target_state_within_limits && (vf > vMax || vf < vMin)) {
                    return false;
                }

                if (check_current_state_within_limits && ((a0 > 0 && v_at_a_zero(v0, a0, jMax) > vMax) || (a0 < 0 && v_at_a_zero(v0, a0, -jMax) < vMin))) {
                    return false;
                }
                if (check_target_state_within_limits && ((af < 0 && v_at_a_zero(vf, af, jMax) > vMax) || (af > 0 && v_at_a_zero(vf, af, -jMax) < vMin))) {
                    return false;
                }
            }
        }

        if (!input.intermediate_positions.empty() && input.control_interface == ControlInterface::Position) {
            if (input.intermediate_positions.size() > max_number_of_waypoints) {
                return false;
            }
            
            if (input.minimum_duration || input.duration_discretization != DurationDiscretization::Continuous) {
                return false;
            }

            if (input.per_dof_control_interface || input.per_dof_synchronization) {
                return false;
            }
        }

        return true;
    }

    //! Calculate a new trajectory for the given input
    Result calculate(const InputParameter<DOFs>& input, Trajectory<DOFs>& trajectory) {
        bool was_interrupted {false};
        return calculate(input, trajectory, was_interrupted);
    }

    //! Calculate a new trajectory for the given input and check for interruption
    Result calculate(const InputParameter<DOFs>& input, Trajectory<DOFs>& trajectory, bool& was_interrupted) {
        if (!validate_input(input, false, true)) {
            return Result::ErrorInvalidInput;
        }

        return calculator.template calculate<throw_error, return_error_at_maximal_duration>(input, trajectory, delta_time, was_interrupted);
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
            Result result = calculate(input, output.trajectory, output.was_calculation_interrupted);
            if (result != Result::Working) {
                return result;
            }

            current_input = input;
            output.time = 0.0;
            output.new_calculation = true;
        }

        const size_t old_section = output.new_section;
        output.time += delta_time;
        output.trajectory.at_time(output.time, output.new_position, output.new_velocity, output.new_acceleration, output.new_section);
        output.did_section_change = (output.new_section != old_section);

        const auto stop = std::chrono::high_resolution_clock::now();
        output.calculation_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / 1000.0;

        output.pass_to_input(current_input);

        if (output.time > output.trajectory.get_duration()) {
            return Result::Finished;
        }

        return Result::Working;
    }
};

} // namespace ruckig
