#pragma once

#include <array>

#include <ruckig/trajectory.hpp>


namespace ruckig {

//! Output type of the OTG
template<size_t DOFs>
struct OutputParameter {
    template<class T> using Vector = std::array<T, DOFs>;
    size_t degrees_of_freedom;

    Vector<double> new_position, new_velocity, new_acceleration;

    //! Was a new trajectory calculation performed in the last cycle?
    bool new_calculation {false};

    //! Computational duration of the last update call
    double calculation_duration; // [µs]

    //! Current trajectory
    Trajectory<DOFs> trajectory;

    //! Current time on trajectory
    double time;

    template <class = typename std::enable_if<DOFs >= 1>::type>
    OutputParameter(): degrees_of_freedom(DOFs) { }
};

} // namespace ruckig
