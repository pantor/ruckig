#include <vector>

#include <emscripten/bind.h>

#include <ruckig/ruckig.hpp>


namespace em = emscripten;
using namespace ruckig;

struct TrajectoryState {
    size_t section;
    std::vector<double> position;
    std::vector<double> velocity;
    std::vector<double> acceleration;
    std::vector<double> jerk;
};


EMSCRIPTEN_BINDINGS(ruckig) {
    em::register_vector<double>("Vector");

    em::enum_<ControlInterface>("ControlInterface")
        .value("Position", ControlInterface::Position)
        .value("Velocity", ControlInterface::Velocity);

    em::enum_<Synchronization>("Synchronization")
        .value("Phase", Synchronization::Phase)
        .value("Time", Synchronization::Time)
        .value("TimeIfNecessary", Synchronization::TimeIfNecessary)
        .value("No", Synchronization::None);

    em::enum_<DurationDiscretization>("DurationDiscretization")
        .value("Continuous", DurationDiscretization::Continuous)
        .value("Discrete", DurationDiscretization::Discrete);

    em::enum_<Result>("Result")
        .value("Working", Result::Working)
        .value("Finished", Result::Finished)
        .value("Error", Result::Error)
        .value("ErrorInvalidInput", Result::ErrorInvalidInput)
        .value("ErrorPositionalLimits", Result::ErrorPositionalLimits)
        .value("ErrorExecutionTimeCalculation", Result::ErrorExecutionTimeCalculation)
        .value("ErrorSynchronizationCalculation", Result::ErrorSynchronizationCalculation);

    em::class_<TrajectoryState>("TrajectoryState")
        .property("section", &TrajectoryState::section)
        .property("position", &TrajectoryState::position)
        .property("velocity", &TrajectoryState::velocity)
        .property("acceleration", &TrajectoryState::acceleration)
        .property("jerk", &TrajectoryState::jerk);

    em::class_<Trajectory<DynamicDOFs>>("Trajectory")
        .constructor<size_t>()
        .property("degrees_of_freedom", &Trajectory<DynamicDOFs>::degrees_of_freedom)
        .function("get_duration", &Trajectory<DynamicDOFs>::get_duration)
        .function("get_intermediate_durations", &Trajectory<DynamicDOFs>::get_intermediate_durations)
        .function("get_independent_min_durations", &Trajectory<DynamicDOFs>::get_independent_min_durations)
        .function("get_position_extrema", &Trajectory<DynamicDOFs>::get_position_extrema)
        .function("at_time", em::select_overload<TrajectoryState(const Trajectory<DynamicDOFs>&, double)>([](const Trajectory<DynamicDOFs>& traj, double time) {
            TrajectoryState result;
            result.position.resize(traj.degrees_of_freedom);
            result.velocity.resize(traj.degrees_of_freedom);
            result.acceleration.resize(traj.degrees_of_freedom);
            result.jerk.resize(traj.degrees_of_freedom);
            traj.at_time(time, result.position, result.velocity, result.acceleration, result.jerk, result.section);
            return result;
        }));

    em::class_<InputParameter<DynamicDOFs>>("InputParameter")
        .constructor<size_t>()
        .property("degrees_of_freedom", &InputParameter<DynamicDOFs>::degrees_of_freedom)
        .property("current_position", &InputParameter<DynamicDOFs>::current_position)
        .property("current_velocity", &InputParameter<DynamicDOFs>::current_velocity)
        .property("current_acceleration", &InputParameter<DynamicDOFs>::current_acceleration)
        .property("target_position", &InputParameter<DynamicDOFs>::target_position)
        .property("target_velocity", &InputParameter<DynamicDOFs>::target_velocity)
        .property("target_acceleration", &InputParameter<DynamicDOFs>::target_acceleration)
        .property("max_velocity", &InputParameter<DynamicDOFs>::max_velocity)
        .property("max_acceleration", &InputParameter<DynamicDOFs>::max_acceleration)
        .property("max_jerk", &InputParameter<DynamicDOFs>::max_jerk)
        .property("min_velocity", &InputParameter<DynamicDOFs>::min_velocity)
        .property("min_acceleration", &InputParameter<DynamicDOFs>::min_acceleration)
        .property("intermediate_positions", &InputParameter<DynamicDOFs>::intermediate_positions)
        .property("per_section_max_velocity", &InputParameter<DynamicDOFs>::per_section_max_velocity)
        .property("per_section_max_acceleration", &InputParameter<DynamicDOFs>::per_section_max_acceleration)
        .property("per_section_max_jerk", &InputParameter<DynamicDOFs>::per_section_max_jerk)
        .property("per_section_min_velocity", &InputParameter<DynamicDOFs>::per_section_min_velocity)
        .property("per_section_min_acceleration", &InputParameter<DynamicDOFs>::per_section_min_acceleration)
        .property("max_position", &InputParameter<DynamicDOFs>::max_position)
        .property("min_position", &InputParameter<DynamicDOFs>::min_position)
        .property("enabled", &InputParameter<DynamicDOFs>::enabled)
        .property("control_interface", &InputParameter<DynamicDOFs>::control_interface)
        .property("synchronization", &InputParameter<DynamicDOFs>::synchronization)
        .property("duration_discretization", &InputParameter<DynamicDOFs>::duration_discretization)
        .property("per_dof_control_interface", &InputParameter<DynamicDOFs>::per_dof_control_interface)
        .property("per_dof_synchronization", &InputParameter<DynamicDOFs>::per_dof_synchronization)
        .property("minimum_duration", &InputParameter<DynamicDOFs>::minimum_duration)
        .property("per_section_minimum_duration", &InputParameter<DynamicDOFs>::per_section_minimum_duration)
        .property("interrupt_calculation_duration", &InputParameter<DynamicDOFs>::interrupt_calculation_duration)
        .function("validate", &InputParameter<DynamicDOFs>::validate<true>);

    em::class_<Ruckig<DynamicDOFs>>("Ruckig")
        .constructor<size_t>()
        .constructor<size_t, double>()
        .property("max_number_of_waypoints", &Ruckig<DynamicDOFs>::max_number_of_waypoints)
        .property("degrees_of_freedom", &Ruckig<DynamicDOFs>::degrees_of_freedom)
        .property("delta_time", &Ruckig<DynamicDOFs>::delta_time)
        .function("reset", &Ruckig<DynamicDOFs>::reset)
        .function("validate_input", &Ruckig<DynamicDOFs>::validate_input<true>)
        .function("calculate", static_cast<Result (Ruckig<DynamicDOFs>::*)(const InputParameter<DynamicDOFs>&, Trajectory<DynamicDOFs>&)>(&Ruckig<DynamicDOFs>::calculate))
        .function("update", static_cast<Result (Ruckig<DynamicDOFs>::*)(const InputParameter<DynamicDOFs>&, OutputParameter<DynamicDOFs>&)>(&Ruckig<DynamicDOFs>::update));
}
