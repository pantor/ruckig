#pragma once

#include <ruckig/calculator_target.hpp>
#ifdef WITH_ONLINE_CLIENT
#include <ruckig/calculator_online.hpp>
#endif
#include <ruckig/input_parameter.hpp>
#include <ruckig/trajectory.hpp>


namespace ruckig {

//! Calculation interface.
template<size_t DOFs, template<class, size_t> class CustomVector = StandardVector>
class Calculator {
    using InputParameter = InputParameter<DOFs, CustomVector>;
    using Trajectory = Trajectory<DOFs, CustomVector>;
    using TargetCalculator = TargetCalculator<DOFs, CustomVector>;
    using WaypointsCalculator = WaypointsCalculator<DOFs, CustomVector>;

    TargetCalculator target_calculator;
#if defined WITH_ONLINE_CLIENT
    WaypointsCalculator waypoints_calculator;
#endif

    inline bool use_waypoints_trajectory(const InputParameter& input) {
        return !input.intermediate_positions.empty() && input.control_interface == ControlInterface::Position;
    }

public:
    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    explicit Calculator() { }

#if defined WITH_ONLINE_CLIENT
    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    explicit Calculator(size_t max_waypoints): waypoints_calculator(WaypointsCalculator(max_waypoints)) { }

    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    explicit Calculator(size_t dofs): target_calculator(TargetCalculator(dofs)), waypoints_calculator(WaypointsCalculator(dofs)) { }

    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    explicit Calculator(size_t dofs, size_t max_waypoints): target_calculator(TargetCalculator(dofs)), waypoints_calculator(WaypointsCalculator(dofs, max_waypoints)) { }
#else
    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    explicit Calculator(size_t dofs): target_calculator(TargetCalculator(dofs)) { }
#endif

    //! Calculate the time-optimal waypoint-based trajectory
    template<bool throw_error>
    Result calculate(const InputParameter& input, Trajectory& trajectory, double delta_time, bool& was_interrupted) {
        Result result;
#if defined WITH_ONLINE_CLIENT
        if (use_waypoints_trajectory(input)) {
            result = waypoints_calculator.template calculate<throw_error>(input, trajectory, delta_time, was_interrupted);
        } else {
            result = target_calculator.template calculate<throw_error>(input, trajectory, delta_time, was_interrupted);
        }
#else
        result = target_calculator.template calculate<throw_error>(input, trajectory, delta_time, was_interrupted);
#endif

        return result;
    }

    //! Continue the trajectory calculation
    template<bool throw_error>
    Result continue_calculation(const InputParameter& input, Trajectory& trajectory, double delta_time, bool& was_interrupted) {
        Result result;
#if defined WITH_ONLINE_CLIENT
        if (use_waypoints_trajectory(input)) {
            result = waypoints_calculator.template continue_calculation<throw_error>(input, trajectory, delta_time, was_interrupted);
        } else {
            result = target_calculator.template continue_calculation<throw_error>(input, trajectory, delta_time, was_interrupted);
        }
#else
        result = target_calculator.template continue_calculation<throw_error>(input, trajectory, delta_time, was_interrupted);
#endif

        return result;
    }
};

} // namespace ruckig
