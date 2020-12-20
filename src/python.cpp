#include <array>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
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
    m.doc() = "Online Trajectory Generation";

    constexpr size_t DOFs {1};

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
        .def_readwrite("minimum_duration", &InputParameter<DOFs>::minimum_duration);

    py::class_<OutputParameter<DOFs>>(m, "OutputParameter")
        .def(py::init<>())
        .def_readonly_static("degrees_of_freedom", &InputParameter<DOFs>::degrees_of_freedom)
        .def_readwrite("new_position", &OutputParameter<DOFs>::new_position)
        .def_readwrite("new_velocity", &OutputParameter<DOFs>::new_velocity)
        .def_readwrite("new_acceleration", &OutputParameter<DOFs>::new_acceleration)
        .def_readwrite("duration", &OutputParameter<DOFs>::duration)
        .def_readwrite("new_calculation", &OutputParameter<DOFs>::new_calculation)
        .def_readwrite("calculation_duration", &OutputParameter<DOFs>::calculation_duration)
        .def("__copy__",  [](const OutputParameter<DOFs> &self) {
            return OutputParameter<DOFs>(self);
        });

    py::enum_<Result>(m, "Result")
        .value("Working", Result::Working)
        .value("Finished", Result::Finished)
        .value("Error", Result::Error)
        .export_values();

    py::class_<Quintic<DOFs>>(m, "Quintic")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Quintic<DOFs>::delta_time)
        .def("update", &Quintic<DOFs>::update);

    py::class_<Smoothie<DOFs>>(m, "Smoothie")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Smoothie<DOFs>::delta_time)
        .def("update", &Smoothie<DOFs>::update);

    py::class_<Ruckig<DOFs>>(m, "Ruckig")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Ruckig<DOFs>::delta_time)
        .def("update", &Ruckig<DOFs>::update)
        .def("at_time", &Ruckig<DOFs>::atTime);

#ifdef WITH_REFLEXXES
    py::class_<Reflexxes<DOFs>>(m, "Reflexxes")
        .def(py::init<double>(), "delta_time"_a)
        .def_readonly("delta_time", &Reflexxes<DOFs>::delta_time)
        .def("update", &Reflexxes<DOFs>::update)
        .def("at_time", &Reflexxes<DOFs>::atTime);
#endif
}
