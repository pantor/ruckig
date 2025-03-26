#include <array>
#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/operators.h>

#include <ruckig/ruckig.hpp>


namespace nb = nanobind;
using namespace nb::literals; // to bring in the `_a` literal
using namespace ruckig;


NB_MODULE(ruckig, m) {
    m.doc() = "Instantaneous Motion Generation for Robots and Machines. Real-time and time-optimal trajectory calculation \
given a target waypoint with position, velocity, and acceleration, starting from any initial state \
limited by velocity, acceleration, and jerk constraints.";
    m.attr("__version__")  = "0.14.0";

    nb::enum_<ControlInterface>(m, "ControlInterface")
        .value("Position", ControlInterface::Position)
        .value("Velocity", ControlInterface::Velocity)
        .export_values();

    nb::enum_<Synchronization>(m, "Synchronization")
        .value("Phase", Synchronization::Phase)
        .value("Time", Synchronization::Time)
        .value("TimeIfNecessary", Synchronization::TimeIfNecessary)
        .value("No", Synchronization::None)
        .export_values();

    nb::enum_<DurationDiscretization>(m, "DurationDiscretization")
        .value("Continuous", DurationDiscretization::Continuous)
        .value("Discrete", DurationDiscretization::Discrete)
        .export_values();

    nb::enum_<Result>(m, "Result", nb::is_arithmetic())
        .value("Working", Result::Working)
        .value("Finished", Result::Finished)
        .value("Error", Result::Error)
        .value("ErrorInvalidInput", Result::ErrorInvalidInput)
        .value("ErrorPositionalLimits", Result::ErrorPositionalLimits)
        .value("ErrorExecutionTimeCalculation", Result::ErrorExecutionTimeCalculation)
        .value("ErrorSynchronizationCalculation", Result::ErrorSynchronizationCalculation)
        .export_values();

    nb::exception<RuckigError>(m, "RuckigError");

    nb::class_<Bound>(m, "Bound")
        .def_ro("min", &Bound::min)
        .def_ro("max", &Bound::max)
        .def_ro("t_min", &Bound::t_min)
        .def_ro("t_max", &Bound::t_max)
        .def("__repr__", [](const Bound& ext) {
            return "[" + std::to_string(ext.min) + ", " + std::to_string(ext.max) + "]";
        });

    nb::class_<Trajectory<DynamicDOFs>>(m, "Trajectory")
        .def(nb::init<size_t>(), "dofs"_a)
#if defined WITH_CLOUD_CLIENT
        .def(nb::init<size_t, size_t>(), "dofs"_a, "max_number_of_waypoints"_a)
#endif
        .def_ro("degrees_of_freedom", &Trajectory<DynamicDOFs>::degrees_of_freedom)
        .def_prop_ro("profiles", &Trajectory<DynamicDOFs>::get_profiles)
        .def_prop_ro("duration", &Trajectory<DynamicDOFs>::get_duration)
        .def_prop_ro("intermediate_durations", &Trajectory<DynamicDOFs>::get_intermediate_durations)
        .def_prop_ro("independent_min_durations", &Trajectory<DynamicDOFs>::get_independent_min_durations)
        .def_prop_ro("position_extrema", &Trajectory<DynamicDOFs>::get_position_extrema)
        .def("at_time", [](const Trajectory<DynamicDOFs>& traj, double time, bool return_section=false) {
            std::vector<double> new_position(traj.degrees_of_freedom), new_velocity(traj.degrees_of_freedom), new_acceleration(traj.degrees_of_freedom), new_jerk(traj.degrees_of_freedom);
            size_t new_section;
            traj.at_time(time, new_position, new_velocity, new_acceleration, new_jerk, new_section);
            if (return_section) {
                return nb::make_tuple(new_position, new_velocity, new_acceleration, new_section);
            }
            return nb::make_tuple(new_position, new_velocity, new_acceleration);
        }, "time"_a, "return_section"_a=false)
        .def("get_first_time_at_position", &Trajectory<DynamicDOFs>::get_first_time_at_position, "dof"_a, "position"_a, "time_after"_a=0.0);

    nb::class_<InputParameter<DynamicDOFs>>(m, "InputParameter")
        .def(nb::init<size_t>(), "dofs"_a)
#if defined WITH_CLOUD_CLIENT
        .def(nb::init<size_t, size_t>(), "dofs"_a, "max_number_of_waypoints"_a)
#endif
        .def_ro("degrees_of_freedom", &InputParameter<DynamicDOFs>::degrees_of_freedom)
        .def_rw("current_position", &InputParameter<DynamicDOFs>::current_position)
        .def_rw("current_velocity", &InputParameter<DynamicDOFs>::current_velocity)
        .def_rw("current_acceleration", &InputParameter<DynamicDOFs>::current_acceleration)
        .def_rw("target_position", &InputParameter<DynamicDOFs>::target_position)
        .def_rw("target_velocity", &InputParameter<DynamicDOFs>::target_velocity)
        .def_rw("target_acceleration", &InputParameter<DynamicDOFs>::target_acceleration)
        .def_rw("max_velocity", &InputParameter<DynamicDOFs>::max_velocity)
        .def_rw("max_acceleration", &InputParameter<DynamicDOFs>::max_acceleration)
        .def_rw("max_jerk", &InputParameter<DynamicDOFs>::max_jerk)
        .def_rw("min_velocity", &InputParameter<DynamicDOFs>::min_velocity, nb::arg().none())
        .def_rw("min_acceleration", &InputParameter<DynamicDOFs>::min_acceleration, nb::arg().none())
        .def_rw("intermediate_positions", &InputParameter<DynamicDOFs>::intermediate_positions)
        .def_rw("per_section_max_velocity", &InputParameter<DynamicDOFs>::per_section_max_velocity, nb::arg().none())
        .def_rw("per_section_max_acceleration", &InputParameter<DynamicDOFs>::per_section_max_acceleration, nb::arg().none())
        .def_rw("per_section_max_jerk", &InputParameter<DynamicDOFs>::per_section_max_jerk, nb::arg().none())
        .def_rw("per_section_min_velocity", &InputParameter<DynamicDOFs>::per_section_min_velocity, nb::arg().none())
        .def_rw("per_section_min_acceleration", &InputParameter<DynamicDOFs>::per_section_min_acceleration, nb::arg().none())
        .def_rw("per_section_max_position", &InputParameter<DynamicDOFs>::per_section_max_position, nb::arg().none())
        .def_rw("per_section_min_position", &InputParameter<DynamicDOFs>::per_section_min_position, nb::arg().none())
        .def_rw("max_position", &InputParameter<DynamicDOFs>::max_position, nb::arg().none())
        .def_rw("min_position", &InputParameter<DynamicDOFs>::min_position, nb::arg().none())
        .def_rw("enabled", &InputParameter<DynamicDOFs>::enabled)
        .def_rw("control_interface", &InputParameter<DynamicDOFs>::control_interface)
        .def_rw("synchronization", &InputParameter<DynamicDOFs>::synchronization)
        .def_rw("duration_discretization", &InputParameter<DynamicDOFs>::duration_discretization)
        .def_rw("per_dof_control_interface", &InputParameter<DynamicDOFs>::per_dof_control_interface, nb::arg().none())
        .def_rw("per_dof_synchronization", &InputParameter<DynamicDOFs>::per_dof_synchronization, nb::arg().none())
        .def_rw("minimum_duration", &InputParameter<DynamicDOFs>::minimum_duration, nb::arg().none())
        .def_rw("per_section_minimum_duration", &InputParameter<DynamicDOFs>::per_section_minimum_duration, nb::arg().none())
        .def_rw("interrupt_calculation_duration", &InputParameter<DynamicDOFs>::interrupt_calculation_duration, nb::arg().none())
        .def("validate", &InputParameter<DynamicDOFs>::validate<true>, "check_current_state_within_limits"_a=false, "check_target_state_within_limits"_a=true)
        .def(nb::self != nb::self)
        .def("__repr__", &InputParameter<DynamicDOFs>::to_string);

    nb::class_<OutputParameter<DynamicDOFs>>(m, "OutputParameter")
        .def(nb::init<size_t>(), "dofs"_a)
#if defined WITH_CLOUD_CLIENT
        .def(nb::init<size_t, size_t>(), "dofs"_a, "max_number_of_waypoints"_a)
#endif
        .def_ro("degrees_of_freedom", &OutputParameter<DynamicDOFs>::degrees_of_freedom)
        .def_ro("new_position", &OutputParameter<DynamicDOFs>::new_position)
        .def_ro("new_velocity", &OutputParameter<DynamicDOFs>::new_velocity)
        .def_ro("new_acceleration", &OutputParameter<DynamicDOFs>::new_acceleration)
        .def_ro("new_section", &OutputParameter<DynamicDOFs>::new_section)
        .def_ro("did_section_change", &OutputParameter<DynamicDOFs>::did_section_change)
        .def_ro("trajectory", &OutputParameter<DynamicDOFs>::trajectory)
        .def_rw("time", &OutputParameter<DynamicDOFs>::time)
        .def_ro("new_calculation", &OutputParameter<DynamicDOFs>::new_calculation)
        .def_ro("was_calculation_interrupted", &OutputParameter<DynamicDOFs>::was_calculation_interrupted)
        .def_ro("calculation_duration", &OutputParameter<DynamicDOFs>::calculation_duration)
        .def("pass_to_input", &OutputParameter<DynamicDOFs>::pass_to_input, "input"_a)
        .def("__repr__", &OutputParameter<DynamicDOFs>::to_string)
        .def("__copy__",  [](const OutputParameter<DynamicDOFs> &self) {
            return OutputParameter<DynamicDOFs>(self);
        });

    nb::class_<RuckigThrow<DynamicDOFs>>(m, "Ruckig")
        .def(nb::init<size_t>(), "dofs"_a)
        .def(nb::init<size_t, double>(), "dofs"_a, "delta_time"_a)
#if defined WITH_CLOUD_CLIENT
        .def(nb::init<size_t, double, size_t>(), "dofs"_a, "delta_time"_a, "max_number_of_waypoints"_a=0)
        .def("filter_intermediate_positions", &RuckigThrow<DynamicDOFs>::filter_intermediate_positions, "input"_a, "threshold_distance"_a)
#endif
        .def_ro("max_number_of_waypoints", &RuckigThrow<DynamicDOFs>::max_number_of_waypoints)
        .def_ro("degrees_of_freedom", &RuckigThrow<DynamicDOFs>::degrees_of_freedom)
        .def_rw("delta_time", &RuckigThrow<DynamicDOFs>::delta_time)
        .def("reset", &RuckigThrow<DynamicDOFs>::reset)
        .def("validate_input", &RuckigThrow<DynamicDOFs>::validate_input<true>, "input"_a, "check_current_state_within_limits"_a=false, "check_target_state_within_limits"_a=true)
        .def("calculate", static_cast<Result (RuckigThrow<DynamicDOFs>::*)(const InputParameter<DynamicDOFs>&, Trajectory<DynamicDOFs>&)>(&RuckigThrow<DynamicDOFs>::calculate), "input"_a, "trajectory"_a)
        .def("calculate", static_cast<Result (RuckigThrow<DynamicDOFs>::*)(const InputParameter<DynamicDOFs>&, Trajectory<DynamicDOFs>&, bool&)>(&RuckigThrow<DynamicDOFs>::calculate), "input"_a, "trajectory"_a, "was_interrupted"_a)
        .def("update", static_cast<Result (RuckigThrow<DynamicDOFs>::*)(const InputParameter<DynamicDOFs>&, OutputParameter<DynamicDOFs>&)>(&RuckigThrow<DynamicDOFs>::update), "input"_a, "output"_a);

    nb::class_<BrakeProfile>(m, "BrakeProfile")
        .def_ro("duration", &BrakeProfile::duration)
        .def_ro("t", &BrakeProfile::t)
        .def_ro("j", &BrakeProfile::j)
        .def_ro("a", &BrakeProfile::a)
        .def_ro("v", &BrakeProfile::v)
        .def_ro("p", &BrakeProfile::p);

    nb::class_<Profile>(m, "Profile")
        .def_ro("t", &Profile::t)
        .def_ro("j", &Profile::j)
        .def_ro("a", &Profile::a)
        .def_ro("v", &Profile::v)
        .def_ro("p", &Profile::p)
        .def_ro("brake", &Profile::brake)
        .def_ro("accel", &Profile::accel)
        .def_ro("pf", &Profile::pf)
        .def_ro("vf", &Profile::vf)
        .def_ro("af", &Profile::af)
        .def_ro("limits", &Profile::limits)
        .def_ro("direction", &Profile::direction)
        .def_ro("control_signs", &Profile::control_signs)
        .def("__repr__", &Profile::to_string);
}
