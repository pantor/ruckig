#pragma once

#include <array>

#include <ruckig/trajectory.hpp>


namespace ruckig {

//! Output type of the OTG
template<size_t DOFs>
struct OutputParameter {
    static constexpr size_t degrees_of_freedom {DOFs};

    std::array<double, DOFs> new_position, new_velocity, new_acceleration;

    bool new_calculation {false};
    double calculation_duration; // [µs]

    // Current trajectory
    Trajectory<DOFs> trajectory;

    // Current time on trajectory
    double time;
};

} // namespace ruckig
