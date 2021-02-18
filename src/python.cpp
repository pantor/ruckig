#include <array>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <ruckig/ruckig.hpp>

#include <ruckig/alternative/quintic.hpp>
#include <ruckig/alternative/smoothie.hpp>

#ifdef WITH_REFLEXXES
    #include <ruckig/alternative/reflexxes.hpp>
#endif


namespace py = pybind11;
using namespace pybind11::literals; // to bring in the `_a` literal
using namespace ruckig;


PYBIND11_MODULE(_ruckig, m) {
    m.doc() = "Online Trajectory Generation. Real-time and time-optimal trajectory calculation \
given a target waypoint with position, velocity, and acceleration, starting from any initial state \
limited by velocity, acceleration, and jerk constraints.";

    constexpr size_t DOFs {3};
    using IP = InputParameter<DOFs>;
    using OP = OutputParameter<DOFs>;

    py::class_<IP> input_parameter(m, "InputParameter");
    input_parameter.def(py::init<>())
        .def_readonly_static("degrees_of_freedom", &IP::degrees_of_freedom)
        .def_readwrite("current_position", &IP::current_position)
        .def_readwrite("current_velocity", &IP::current_velocity)
        .def_readwrite("current_acceleration", &IP::current_acceleration)
        .def_readwrite("target_position", &IP::target_position)
        .def_readwrite("target_velocity", &IP::target_velocity)
        .def_readwrite("target_acceleration", &IP::target_acceleration)
        .def_readwrite("max_velocity", &IP::max_velocity)
        .def_readwrite("max_acceleration", &IP::max_acceleration)
        .def_readwrite("max_jerk", &IP::max_jerk)
        .def_readwrite("enabled", &IP::enabled)
        .def_readwrite("minimum_duration", &IP::minimum_duration)
        .def_readwrite("min_velocity", &IP::min_velocity)
        .def_readwrite("min_acceleration", &IP::min_acceleration)
        .def_readwrite("synchronization", &IP::synchronization)
        .def(py::self != py::self)
        .def("__repr__", static_cast<std::string (IP::*)() const>(&IP::to_string));

    py::enum_<IP::Type>(input_parameter, "ParameterType")
        .value("Position", IP::Type::Position)
        .value("Velocity", IP::Type::Velocity)
        .export_values();

    py::enum_<IP::Synchronization>(input_parameter, "Synchronization")
        .value("Time", IP::Synchronization::Time)
        .value("TimeIfNecessary", IP::Synchronization::TimeIfNecessary)
        .value("No", IP::Synchronization::None)
        .export_values();

    py::class_<Trajectory<DOFs>>(m, "Trajectory")
        .def_readonly("duration", &Trajectory<DOFs>::duration)
        .def_readonly("independent_min_durations", &Trajectory<DOFs>::independent_min_durations)
        .def("at_time", &Trajectory<DOFs>::at_time)
        .def("get_position_extrema", &Trajectory<DOFs>::get_position_extrema);

    py::class_<OP>(m, "OutputParameter")
        .def(py::init<>())
        .def_readonly_static("degrees_of_freedom", &IP::degrees_of_freedom)
        .def_readonly("new_position", &OP::new_position)
        .def_readonly("new_velocity", &OP::new_velocity)
        .def_readonly("new_acceleration", &OP::new_acceleration)
        .def_readonly("trajectory", &OP::trajectory)
        .def_readonly("time", &OP::time)
        .def_readonly("new_calculation", &OP::new_calculation)
        .def_readonly("calculation_duration", &OP::calculation_duration)
        .def("__copy__",  [](const OP &self) {
            return OP(self);
        });

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
        .def_readonly("t_max", &PositionExtrema::t_max);

    py::class_<Quintic<DOFs>>(m, "Quintic")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Quintic<DOFs>::delta_time)
        .def_readonly_static("degrees_of_freedom", &Quintic<DOFs>::degrees_of_freedom)
        .def("update", &Quintic<DOFs>::update);

    py::class_<Smoothie<DOFs>>(m, "Smoothie")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Smoothie<DOFs>::delta_time)
        .def_readonly_static("degrees_of_freedom", &Smoothie<DOFs>::degrees_of_freedom)
        .def("update", &Smoothie<DOFs>::update);

    py::class_<Ruckig<DOFs, true>>(m, "Ruckig")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Ruckig<DOFs, true>::delta_time)
        .def_readonly_static("degrees_of_freedom", &Ruckig<DOFs, true>::degrees_of_freedom)
        .def("update", &Ruckig<DOFs, true>::update);

#ifdef WITH_REFLEXXES
    py::class_<Reflexxes<DOFs>>(m, "Reflexxes")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Reflexxes<DOFs>::delta_time)
        .def_readonly_static("degrees_of_freedom", &Reflexxes<DOFs>::degrees_of_freedom)
        .def("update", &Reflexxes<DOFs>::update);
#endif
}
