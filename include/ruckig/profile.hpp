#pragma once

#include <array>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <optional>
#include <tuple>


namespace ruckig {

//! The state profile for position, velocity, acceleration and jerk for a single DoF
struct Profile {
    enum class Limits { ACC0_ACC1_VEL, VEL, ACC0, ACC1, ACC0_ACC1, ACC0_VEL, ACC1_VEL, NONE } limits;
    enum class Direction { UP, DOWN } direction;
    enum class Teeth { UDDU, UDUD } teeth;

    std::array<double, 7> t, t_sum, j;
    std::array<double, 8> a, v, p;

    //! Total time of the braking segments
    std::optional<double> t_brake;

    //! Allow up to two segments of braking before the "correct" profile starts
    std::array<double, 2> t_brakes, j_brakes, a_brakes, v_brakes, p_brakes;

    template<Teeth teeth, Limits limits>
    bool check(double pf, double vf, double af, double jf, double vMax, double aMax) {
        if constexpr (teeth == Teeth::UDDU) {
            j = {jf, 0, -jf, 0, -jf, 0, jf};
        } else {
            j = {jf, 0, -jf, 0, jf, 0, -jf};
        }

        t_sum[0] = t[0];
        if (t[0] < 0) {
            return false;
        }

        for (size_t i = 0; i < 6; i += 1) {
            if (t[i+1] < 0) {
                return false;
            }
            t_sum[i+1] = t_sum[i] + t[i+1];
        }
        for (size_t i = 0; i < 7; i += 1) {
            a[i+1] = a[i] + t[i] * j[i];
            v[i+1] = v[i] + t[i] * (a[i] + t[i] * j[i] / 2);
            p[i+1] = p[i] + t[i] * (v[i] + t[i] * (a[i] / 2 + t[i] * j[i] / 6));

            if constexpr (limits == Limits::ACC0_ACC1_VEL || limits == Limits::ACC0_VEL || limits == Limits::ACC1_VEL || limits == Limits::VEL) {
                if (i == 2) {
                    a[i+1] = 0.0;
                }
            }
        }

        if (t_sum[6] > 1e12) { // For numerical reasons, equal to around 32000 years in SI units...
            return false;
        }

        this->teeth = teeth;
        this->limits = limits;

        // Velocity limit can be broken in the beginning if both initial velocity and acceleration are too high
        // std::cout << std::setprecision(15) << "target: " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << " T: " << t_sum[6] << " " << to_string() << std::endl;
        const double vMaxAbs = std::abs(vMax) + 1e-9;
        const double aMaxAbs = std::abs(aMax) + 1e-9;
        return std::abs(p[7] - pf) < 1e-8
            && std::abs(v[7] - vf) < 1e-8
            && std::abs(a[7] - af) < 1e-12 // This is not really needed, but we want to double check
            && std::abs(v[3]) < vMaxAbs
            && std::abs(v[4]) < vMaxAbs
            && std::abs(v[5]) < vMaxAbs
            && std::abs(v[6]) < vMaxAbs
            && std::abs(a[1]) < aMaxAbs
            && std::abs(a[3]) < aMaxAbs
            && std::abs(a[5]) < aMaxAbs;
    }
    
    template<Teeth teeth, Limits limits>
    inline bool check(double tf, double pf, double vf, double af, double jf, double vMax, double aMax) {
        // std::cout << std::setprecision(15) << "target: " << std::abs(t_sum[6]-tf) << " " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << std::endl;
        // Time doesn't need to be checked as every profile has one tf - ...
        return check<teeth, limits>(pf, vf, af, jf, vMax, aMax); // && (std::abs(t_sum[6] - tf) < 1e-8);
    }
    
    template<Teeth teeth, Limits limits>
    inline bool check(double tf, double pf, double vf, double af, double jf, double vMax, double aMax, double jMax) {
        return (std::abs(jf) < std::abs(jMax) + 1e-12) && check<teeth, limits>(tf, pf, vf, af, jf, vMax, aMax);
    }

    //! Integrate with constant jerk for duration t. Returns new position, new velocity, and new acceleration.
    static std::tuple<double, double, double> integrate(double t, double p0, double v0, double a0, double j, double scale = 1.0)  {
        return {
            scale * (p0 + t * (v0 + t * (a0 / 2 + t * j / 6))),
            scale * (v0 + t * (a0 + t * j / 2)),
            scale * (a0 + t * j),
        };
    }

    std::string to_string() const {
        std::string result;
        switch (direction) {
            case Direction::UP: result += "UP_"; break;
            case Direction::DOWN: result += "DOWN_"; break;
        }
        switch (limits) {
            case Limits::ACC0_ACC1_VEL: result += "ACC0_ACC1_VEL"; break;
            case Limits::VEL: result += "VEL"; break;
            case Limits::ACC0: result += "ACC0"; break;
            case Limits::ACC1: result += "ACC1"; break;
            case Limits::ACC0_ACC1: result += "ACC0_ACC1"; break;
            case Limits::ACC0_VEL: result += "ACC0_VEL"; break;
            case Limits::ACC1_VEL: result += "ACC1_VEL"; break;
            case Limits::NONE: result += "NONE"; break;
        }
        switch (teeth) {
            case Teeth::UDDU: result += "_UDDU"; break;
            case Teeth::UDUD: result += "_UDUD"; break;
        }
        return result; 
    }
};

} // namespace ruckig
