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


template<size_t MAX>
struct PerDOF {
    static const char* append(std::string name, size_t DOFs) {
        return (name.append(std::to_string(DOFs))).c_str();
    }

    template<size_t DOFs>
    static void define(py::module& m) {
        using IP = InputParameter<DOFs>;
        using OP = OutputParameter<DOFs>;

        py::class_<Trajectory<DOFs>>(m, append("Trajectory", DOFs))
            .def_property_readonly("duration", &Trajectory<DOFs>::get_duration)
            .def_property_readonly("independent_min_durations", &Trajectory<DOFs>::get_independent_min_durations)
            .def_property_readonly("position_extrema", &Trajectory<DOFs>::get_position_extrema)
            .def("at_time", [](const Trajectory<DOFs>& traj, double time) {
                std::array<double, DOFs> new_position, new_velocity, new_acceleration;
                traj.at_time(time, new_position, new_velocity, new_acceleration);
                return py::make_tuple(new_position, new_velocity, new_acceleration);
            });

        py::class_<IP>(m, append("InputParameter", DOFs))
            .def(py::init<>())
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
            .def_readwrite("min_velocity", &IP::min_velocity)
            .def_readwrite("min_acceleration", &IP::min_acceleration)
            .def_readwrite("enabled", &IP::enabled)
            .def_readwrite("interface", &IP::interface)
            .def_readwrite("synchronization", &IP::synchronization)
            .def_readwrite("duration_discretization", &IP::duration_discretization)
            .def_readwrite("minimum_duration", &IP::minimum_duration)
            .def(py::self != py::self)
            .def("__repr__", static_cast<std::string (IP::*)() const>(&IP::to_string));

        py::class_<OP>(m, append("OutputParameter", DOFs))
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

        py::class_<Ruckig<DOFs, true>>(m, append("Ruckig", DOFs))
            .def(py::init<double>(), "delta_time"_a)
            .def_readonly("delta_time", &Ruckig<DOFs, true>::delta_time)
            .def_readonly_static("degrees_of_freedom", &Ruckig<DOFs, true>::degrees_of_freedom)
            .def("validate_input", &Ruckig<DOFs, true>::validate_input)
            .def("update", &Ruckig<DOFs, true>::update);

    #ifdef WITH_REFLEXXES
        py::class_<Reflexxes<DOFs>>(m, append("Reflexxes", DOFs))
            .def(py::init<double>(), "delta_time"_a)
            .def_readonly("delta_time", &Reflexxes<DOFs>::delta_time)
            .def_readonly_static("degrees_of_freedom", &Reflexxes<DOFs>::degrees_of_freedom)
            .def("update", &Reflexxes<DOFs>::update);
    #endif

        // Recurse upwards
        if constexpr (DOFs < MAX - 1) PerDOF<MAX>::define<DOFs+1>(m);
    }
};


template<class T>
py::object cast_unique() {
    return py::cast(std::make_unique<T>());
}

template<class T>
py::object cast_unique(double delta_time) {
    return py::cast(std::make_unique<T>(delta_time));
}

py::object handle_dof_error(size_t dofs) {
    throw std::runtime_error("For Python, the number of DOFs needs to be between 1 and 10, but is " + std::to_string(dofs) +  ".");
    return py::none();
}


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
        .def_readonly("t_max", &PositionExtrema::t_max);

    PerDOF<11>::define<1>(m);

    m.def("InputParameter", [](size_t dofs) {
        switch (dofs) {
            case 1: return cast_unique<InputParameter<1>>();
            case 2: return cast_unique<InputParameter<2>>();
            case 3: return cast_unique<InputParameter<3>>();
            case 4: return cast_unique<InputParameter<4>>();
            case 5: return cast_unique<InputParameter<5>>();
            case 6: return cast_unique<InputParameter<6>>();
            case 7: return cast_unique<InputParameter<7>>();
            case 8: return cast_unique<InputParameter<8>>();
            case 9: return cast_unique<InputParameter<9>>();
            case 10: return cast_unique<InputParameter<10>>();
            default: return handle_dof_error(dofs);
        }
    }, "dofs"_a);

    m.def("OutputParameter", [](size_t dofs) {
        switch (dofs) {
            case 1: return cast_unique<OutputParameter<1>>();
            case 2: return cast_unique<OutputParameter<2>>();
            case 3: return cast_unique<OutputParameter<3>>();
            case 4: return cast_unique<OutputParameter<4>>();
            case 5: return cast_unique<OutputParameter<5>>();
            case 6: return cast_unique<OutputParameter<6>>();
            case 7: return cast_unique<OutputParameter<7>>();
            case 8: return cast_unique<OutputParameter<8>>();
            case 9: return cast_unique<OutputParameter<9>>();
            case 10: return cast_unique<OutputParameter<10>>();
            default: return handle_dof_error(dofs);
        }
    }, "dofs"_a);

    m.def("Ruckig", [](size_t dofs, double delta_time) {
        switch (dofs) {
            case 1: return cast_unique<Ruckig<1, true>>(delta_time);
            case 2: return cast_unique<Ruckig<2, true>>(delta_time);
            case 3: return cast_unique<Ruckig<3, true>>(delta_time);
            case 4: return cast_unique<Ruckig<4, true>>(delta_time);
            case 5: return cast_unique<Ruckig<5, true>>(delta_time);
            case 6: return cast_unique<Ruckig<6, true>>(delta_time);
            case 7: return cast_unique<Ruckig<7, true>>(delta_time);
            case 8: return cast_unique<Ruckig<8, true>>(delta_time);
            case 9: return cast_unique<Ruckig<9, true>>(delta_time);
            case 10: return cast_unique<Ruckig<10, true>>(delta_time);
            default: return handle_dof_error(dofs);
        }
    }, "dofs"_a, "delta_time"_a);

#ifdef WITH_REFLEXXES
    m.def("Reflexxes", [](size_t dofs, double delta_time) {
        switch (dofs) {
            case 1: return cast_unique<Reflexxes<1>>(delta_time);
            case 2: return cast_unique<Reflexxes<2>>(delta_time);
            case 3: return cast_unique<Reflexxes<3>>(delta_time);
            case 4: return cast_unique<Reflexxes<4>>(delta_time);
            case 5: return cast_unique<Reflexxes<5>>(delta_time);
            case 6: return cast_unique<Reflexxes<6>>(delta_time);
            case 7: return cast_unique<Reflexxes<7>>(delta_time);
            case 8: return cast_unique<Reflexxes<8>>(delta_time);
            case 9: return cast_unique<Reflexxes<9>>(delta_time);
            case 10: return cast_unique<Reflexxes<10>>(delta_time);
            default: return handle_dof_error(dofs);
        }
    }, "dofs"_a, "delta_time"_a);
#endif
}
