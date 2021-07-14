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

    py::enum_<Interface>(m, "Interface")
        .value("Position", Interface::Position)
        .value("Velocity", Interface::Velocity)
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

    py::class_<Trajectory<0>>(m, "Trajectory")
        .def(py::init<size_t>(), "dofs"_a)
        .def_property_readonly("duration", &Trajectory<0>::get_duration)
        .def_property_readonly("independent_min_durations", &Trajectory<0>::get_independent_min_durations)
        .def_property_readonly("position_extrema", &Trajectory<0>::get_position_extrema)
        .def("at_time", [](const Trajectory<0>& traj, double time) {
            std::vector<double> new_position(traj.degrees_of_freedom), new_velocity(traj.degrees_of_freedom), new_acceleration(traj.degrees_of_freedom);
            traj.at_time(time, new_position, new_velocity, new_acceleration);
            return py::make_tuple(new_position, new_velocity, new_acceleration);
        })
        .def("get_first_time_at_position", [](const Trajectory<0>& traj, size_t dof, double position) -> py::object {
            double time;
            if (traj.get_first_time_at_position(dof, position, time)) {
                return py::cast(time);
            }
            return py::none();
        });

    py::class_<InputParameter<0>>(m, "InputParameter")
        .def(py::init<size_t>(), "dofs"_a)
        .def_readonly("degrees_of_freedom", &InputParameter<0>::degrees_of_freedom)
        .def_readwrite("current_position", &InputParameter<0>::current_position)
        .def_readwrite("current_velocity", &InputParameter<0>::current_velocity)
        .def_readwrite("current_acceleration", &InputParameter<0>::current_acceleration)
        .def_readwrite("target_position", &InputParameter<0>::target_position)
        .def_readwrite("target_velocity", &InputParameter<0>::target_velocity)
        .def_readwrite("target_acceleration", &InputParameter<0>::target_acceleration)
        .def_readwrite("max_velocity", &InputParameter<0>::max_velocity)
        .def_readwrite("max_acceleration", &InputParameter<0>::max_acceleration)
        .def_readwrite("max_jerk", &InputParameter<0>::max_jerk)
        .def_readwrite("min_velocity", &InputParameter<0>::min_velocity)
        .def_readwrite("min_acceleration", &InputParameter<0>::min_acceleration)
        .def_readwrite("enabled", &InputParameter<0>::enabled)
        .def_readwrite("interface", &InputParameter<0>::interface)
        .def_readwrite("synchronization", &InputParameter<0>::synchronization)
        .def_readwrite("duration_discretization", &InputParameter<0>::duration_discretization)
        .def_readwrite("minimum_duration", &InputParameter<0>::minimum_duration)
        .def(py::self != py::self)
        .def("__repr__", static_cast<std::string (InputParameter<0>::*)() const>(&InputParameter<0>::to_string));

    py::class_<OutputParameter<0>>(m, "OutputParameter")
        .def(py::init<size_t>(), "dofs"_a)
        .def_readonly("degrees_of_freedom", &OutputParameter<0>::degrees_of_freedom)
        .def_readonly("new_position", &OutputParameter<0>::new_position)
        .def_readonly("new_velocity", &OutputParameter<0>::new_velocity)
        .def_readonly("new_acceleration", &OutputParameter<0>::new_acceleration)
        .def_readonly("trajectory", &OutputParameter<0>::trajectory)
        .def_readonly("time", &OutputParameter<0>::time)
        .def_readonly("new_calculation", &OutputParameter<0>::new_calculation)
        .def_readonly("calculation_duration", &OutputParameter<0>::calculation_duration)
        .def("__copy__",  [](const OutputParameter<0> &self) {
            return OutputParameter<0>(self);
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
    py::class_<Reflexxes<0>>(m, "Reflexxes")
        .def(py::init<size_t, double>(), "dofs"_a, "delta_time"_a)
        .def_readonly("degrees_of_freedom", &Reflexxes<0>::degrees_of_freedom)
        .def_readonly("delta_time", &Reflexxes<0>::delta_time)
        .def("update", &Reflexxes<0>::update);
#endif
}
