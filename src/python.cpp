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

    py::enum_<InputParameter<DOFs>::Type>(m, "InputParameterType", py::arithmetic())
        .value("Position", InputParameter<DOFs>::Type::Position)
        .value("Velocity", InputParameter<DOFs>::Type::Velocity)
        .export_values();

    py::class_<InputParameter<DOFs>>(m, "InputParameter")
        .def(py::init<>())
        .def_readonly_static("degrees_of_freedom", &InputParameter<DOFs>::degrees_of_freedom)
        .def_readwrite("current_position", &InputParameter<DOFs>::current_position)
        .def_readwrite("current_velocity", &InputParameter<DOFs>::current_velocity)
        .def_readwrite("current_acceleration", &InputParameter<DOFs>::current_acceleration)
        .def_readwrite("target_position", &InputParameter<DOFs>::target_position)
        .def_readwrite("target_velocity", &InputParameter<DOFs>::target_velocity)
        .def_readwrite("target_acceleration", &InputParameter<DOFs>::target_acceleration)
        .def_readwrite("max_velocity", &InputParameter<DOFs>::max_velocity)
        .def_readwrite("max_acceleration", &InputParameter<DOFs>::max_acceleration)
        .def_readwrite("max_jerk", &InputParameter<DOFs>::max_jerk)
        .def_readwrite("enabled", &InputParameter<DOFs>::enabled)
        .def_readwrite("minimum_duration", &InputParameter<DOFs>::minimum_duration)
        .def_readwrite("min_velocity", &InputParameter<DOFs>::min_velocity)
        .def_readwrite("min_acceleration", &InputParameter<DOFs>::min_acceleration)
        .def(py::self != py::self)
        .def("__repr__", static_cast<std::string (InputParameter<DOFs>::*)() const>(&InputParameter<DOFs>::to_string));

    py::class_<Trajectory<DOFs>>(m, "Trajectory")
        .def_readonly("duration", &Trajectory<DOFs>::duration)
        .def_readonly("independent_min_durations", &Trajectory<DOFs>::independent_min_durations)
        .def("at_time", &Trajectory<DOFs>::at_time)
        .def("get_position_extrema", &Trajectory<DOFs>::get_position_extrema);

    py::class_<OutputParameter<DOFs>>(m, "OutputParameter")
        .def(py::init<>())
        .def_readonly_static("degrees_of_freedom", &InputParameter<DOFs>::degrees_of_freedom)
        .def_readonly("new_position", &OutputParameter<DOFs>::new_position)
        .def_readonly("new_velocity", &OutputParameter<DOFs>::new_velocity)
        .def_readonly("new_acceleration", &OutputParameter<DOFs>::new_acceleration)
        .def_readonly("trajectory", &OutputParameter<DOFs>::trajectory)
        .def_readonly("new_calculation", &OutputParameter<DOFs>::new_calculation)
        .def_readonly("calculation_duration", &OutputParameter<DOFs>::calculation_duration)
        .def("__copy__",  [](const OutputParameter<DOFs> &self) {
            return OutputParameter<DOFs>(self);
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
        .def("update", &Quintic<DOFs>::update);

    py::class_<Smoothie<DOFs>>(m, "Smoothie")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Smoothie<DOFs>::delta_time)
        .def("update", &Smoothie<DOFs>::update);

    py::class_<Ruckig<DOFs, true>>(m, "Ruckig")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Ruckig<DOFs, true>::delta_time)
        .def("update", &Ruckig<DOFs, true>::update);

#ifdef WITH_REFLEXXES
    py::class_<Reflexxes<DOFs>>(m, "Reflexxes")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Reflexxes<DOFs>::delta_time)
        .def("update", &Reflexxes<DOFs>::update);
#endif
}
