#pragma once

#include <array>
#include <cmath>
#include <iostream>


namespace ruckig {

//! Calculates (pre-) trajectory to get current state below the limits
class Brake {
    static constexpr double eps {2.2e-14};

    static void acceleration_brake(double v0, double a0, double vMax, double vMin, double aMax, double aMin, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);
    static void velocity_brake(double v0, double a0, double vMax, double vMin, double aMax, double aMin, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);

public:
    static void get_position_brake_trajectory(double v0, double a0, double vMax, double vMin, double aMax, double aMin, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);
    static void get_velocity_brake_trajectory(double a0, double aMax, double aMin, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake);
};

} // namespace ruckig
