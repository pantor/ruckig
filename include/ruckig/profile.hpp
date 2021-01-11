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

    template<Teeth teeth>
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
        }

        this->teeth = teeth;

        // Velocity and acceleration limits can be broken in the beginning if the initial velocity and acceleration are too high
        // std::cout << std::setprecision(15) << "target: " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << std::endl;
        return std::all_of(v.begin() + 3, v.end(), [vMax](double vm){ return std::abs(vm) < std::abs(vMax) + 1e-9; })
            && std::all_of(a.begin() + 2, a.end(), [aMax](double am){ return std::abs(am) < std::abs(aMax) + 1e-9; })
            && std::abs(p[7] - pf) < 1e-8 && std::abs(v[7] - vf) < 1e-8 && std::abs(a[7] - af) < 1e-8;
    }
    
    template<Teeth teeth>
    inline bool check(double tf, double pf, double vf, double af, double jf, double vMax, double aMax) {
        // std::cout << std::setprecision(15) << "target: " << std::abs(t_sum[6]-tf) << " " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << std::endl;
        return check<teeth>(pf, vf, af, jf, vMax, aMax) && (std::abs(t_sum[6] - tf) < 1e-8);
    }
    
    template<Teeth teeth>
    inline bool check(double tf, double pf, double vf, double af, double jf, double vMax, double aMax, double jMax) {
        return (std::abs(jf) < std::abs(jMax) + 1e-12) && check<teeth>(tf, pf, vf, af, jf, vMax, aMax);
    }

    //! Integrate with constant jerk for duration t. Returns new position, new velocity, and new acceleration.
    static std::tuple<double, double, double> integrate(double t, double p0, double v0, double a0, double j)  {
        const double p_new = p0 + t * (v0 + t * (a0 / 2 + t * j / 6));
        const double v_new = v0 + t * (a0 + t * j / 2);
        const double a_new = a0 + t * j;
        return {p_new, v_new, a_new};
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
