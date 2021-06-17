#pragma once

#include <array>

#include <ruckig/trajectory.hpp>


namespace ruckig {

//! Output type of the OTG
template<size_t DOFs>
struct OutputParameter {
    template<class T> using Vector = std::array<T, DOFs>;
    static constexpr size_t degrees_of_freedom {DOFs};

    Vector<double> new_position, new_velocity, new_acceleration;

    //! Was a new trajectory calculation performed in the last cycle?
    bool new_calculation {false};

    //! Computational duration of the last update call
    double calculation_duration; // [µs]

    //! Current trajectory
    Trajectory<DOFs> trajectory;

    //! Current time on trajectory
    double time;
};

} // namespace ruckig
