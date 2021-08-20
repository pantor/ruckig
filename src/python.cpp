#include <array>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <ruckig/ruckig.hpp>

#ifdef WITH_REFLEXXES
    #include <ruckig/reflexxes_comparison.hpp>
#endif


namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal
using namespace ruckig;


PYBIND11_MODULE(ruckig, m) {
    m.doc() = "Online Trajectory Generation. Real-time and time-optimal trajectory calculation \
given a target waypoint with position, velocity, and acceleration, starting from any initial state \
limited by velocity, acceleration, and jerk constraints.";

    py::enum_<ControlInterface>(m, "ControlInterface")
        .value("Position", ControlInterface::Position)
        .value("Velocity", ControlInterface::Velocity)
        .export_values();

    py::enum_<Synchronization>(m, "Synchronization")
        .value("Phase", Synchronization::Phase)
        .value("Time", Synchronization::Time)
        .value("TimeIfNecessary", Synchronization::TimeIfNecessary)
        .value("No", Synchronization::None)
        .export_values();

    py::enum_<DurationDiscretization>(m, "DurationDiscretization")
        .value("Continuous", DurationDiscretization::Continuous)
        .value("Discrete", DurationDiscretization::Discrete)
        .export_values();

    py::enum_<Result>(m, "Result", py::arithmetic())
        .value("Working", Result::Working)
        .value("Finished", Result::Finished)
        .value("Error", Result::Error)
        .value("ErrorInvalidInput", Result::ErrorInvalidInput)
        .value("ErrorExecutionTimeCalculation", Result::ErrorExecutionTimeCalculation)
        .value("ErrorSynchronizationCalculation", Result::ErrorSynchronizationCalculation)
        .export_values();

    py::class_<PositionExtrema>(m, "PositionExtrema")
        .def_readonly("min", &PositionExtrema::min)
        .def_readonly("max", &PositionExtrema::max)
        .def_readonly("t_min", &PositionExtrema::t_min)
        .def_readonly("t_max", &PositionExtrema::t_max)
        .def("__repr__", [](const PositionExtrema& ext) {
            return "[" + std::to_string(ext.min) + ", " + std::to_string(ext.max) + "]";
        });

    py::class_<Trajectory<DynamicDOFs>>(m, "Trajectory")
        .def(py::init<size_t>(), "dofs"_a)
        .def_readonly("degrees_of_freedom", &Trajectory<DynamicDOFs>::degrees_of_freedom)
        .def_property_readonly("duration", &Trajectory<DynamicDOFs>::get_duration)
        .def_property_readonly("intermediate_durations", &Trajectory<DynamicDOFs>::get_intermediate_durations)
        .def_property_readonly("independent_min_durations", &Trajectory<DynamicDOFs>::get_independent_min_durations)
        .def_property_readonly("position_extrema", &Trajectory<DynamicDOFs>::get_position_extrema)
        .def("at_time", [](const Trajectory<DynamicDOFs>& traj, double time) {
            std::vector<double> new_position(traj.degrees_of_freedom), new_velocity(traj.degrees_of_freedom), new_acceleration(traj.degrees_of_freedom);
            traj.at_time(time, new_position, new_velocity, new_acceleration);
            return py::make_tuple(new_position, new_velocity, new_acceleration);
        })
        .def("get_first_time_at_position", [](const Trajectory<DynamicDOFs>& traj, size_t dof, double position) -> py::object {
            double time;
            if (traj.get_first_time_at_position(dof, position, time)) {
                return py::cast(time);
            }
            return py::none();
        });

    py::class_<InputParameter<DynamicDOFs>>(m, "InputParameter")
        .def(py::init<size_t>(), "dofs"_a)
        .def_readonly("degrees_of_freedom", &InputParameter<DynamicDOFs>::degrees_of_freedom)
        .def_readwrite("current_position", &InputParameter<DynamicDOFs>::current_position)
        .def_readwrite("current_velocity", &InputParameter<DynamicDOFs>::current_velocity)
        .def_readwrite("current_acceleration", &InputParameter<DynamicDOFs>::current_acceleration)
        .def_readwrite("target_position", &InputParameter<DynamicDOFs>::target_position)
        .def_readwrite("target_velocity", &InputParameter<DynamicDOFs>::target_velocity)
        .def_readwrite("target_acceleration", &InputParameter<DynamicDOFs>::target_acceleration)
        .def_readwrite("max_velocity", &InputParameter<DynamicDOFs>::max_velocity)
        .def_readwrite("max_acceleration", &InputParameter<DynamicDOFs>::max_acceleration)
        .def_readwrite("max_jerk", &InputParameter<DynamicDOFs>::max_jerk)
        .def_readwrite("min_velocity", &InputParameter<DynamicDOFs>::min_velocity)
        .def_readwrite("min_acceleration", &InputParameter<DynamicDOFs>::min_acceleration)
        .def_readwrite("intermediate_positions", &InputParameter<DynamicDOFs>::intermediate_positions)
        .def_readwrite("enabled", &InputParameter<DynamicDOFs>::enabled)
        .def_readwrite("control_interface", &InputParameter<DynamicDOFs>::control_interface)
        .def_readwrite("synchronization", &InputParameter<DynamicDOFs>::synchronization)
        .def_readwrite("duration_discretization", &InputParameter<DynamicDOFs>::duration_discretization)
        .def_readwrite("minimum_duration", &InputParameter<DynamicDOFs>::minimum_duration)
        .def_readwrite("interrupt_calculation_duration", &InputParameter<DynamicDOFs>::interrupt_calculation_duration)
        .def(py::self != py::self)
        .def("__repr__", static_cast<std::string (InputParameter<DynamicDOFs>::*)() const>(&InputParameter<DynamicDOFs>::to_string));

    py::class_<OutputParameter<DynamicDOFs>>(m, "OutputParameter")
        .def(py::init<size_t>(), "dofs"_a)
        .def_readonly("degrees_of_freedom", &OutputParameter<DynamicDOFs>::degrees_of_freedom)
        .def_readonly("new_position", &OutputParameter<DynamicDOFs>::new_position)
        .def_readonly("new_velocity", &OutputParameter<DynamicDOFs>::new_velocity)
        .def_readonly("new_acceleration", &OutputParameter<DynamicDOFs>::new_acceleration)
        .def_readonly("trajectory", &OutputParameter<DynamicDOFs>::trajectory)
        .def_readonly("time", &OutputParameter<DynamicDOFs>::time)
        .def_readonly("new_calculation", &OutputParameter<DynamicDOFs>::new_calculation)
        .def_readonly("was_calculation_interrupted", &OutputParameter<DynamicDOFs>::was_calculation_interrupted)
        .def_readonly("calculation_duration", &OutputParameter<DynamicDOFs>::calculation_duration)
        .def("__copy__",  [](const OutputParameter<DynamicDOFs> &self) {
            return OutputParameter<DynamicDOFs>(self);
        });

    py::class_<Ruckig<0, true>>(m, "Ruckig")
        .def(py::init<size_t>(), "dofs"_a)
        .def(py::init<size_t, double>(), "dofs"_a, "delta_time"_a)
        .def_readonly("delta_time", &Ruckig<0, true>::delta_time)
        .def_readonly("degrees_of_freedom", &Ruckig<0, true>::degrees_of_freedom)
        .def("validate_input", &Ruckig<0, true>::validate_input)
        .def("calculate", &Ruckig<0, true>::calculate)
        .def("update", &Ruckig<0, true>::update);

#ifdef WITH_REFLEXXES
    py::class_<Reflexxes<DynamicDOFs>>(m, "Reflexxes")
        .def(py::init<size_t, double>(), "dofs"_a, "delta_time"_a)
        .def_readonly("degrees_of_freedom", &Reflexxes<DynamicDOFs>::degrees_of_freedom)
        .def_readonly("delta_time", &Reflexxes<DynamicDOFs>::delta_time)
        .def("update", &Reflexxes<DynamicDOFs>::update);
#endif
}
