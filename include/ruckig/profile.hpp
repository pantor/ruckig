#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <tuple>

#include <ruckig/roots.hpp>


namespace ruckig {

//! Information about the position extrema
struct PositionExtrema {
    double min, max;
    double t_min, t_max;
};


//! The state profile for position, velocity, acceleration and jerk for a single DoF
struct Profile {
    enum class Limits { ACC0_ACC1_VEL, VEL, ACC0, ACC1, ACC0_ACC1, ACC0_VEL, ACC1_VEL, NONE } limits;
    enum class Direction { UP, DOWN } direction;
    enum class JerkSigns { UDDU, UDUD } jerk_signs;

    std::array<double, 7> t, t_sum, j;
    std::array<double, 8> a, v, p;

    //! Target (final) kinematic state
    double pf, vf, af;

    //! Total time of the braking segments
    std::optional<double> t_brake;

    //! Allow up to two segments of braking before the "correct" profile starts
    std::array<double, 2> t_brakes, j_brakes, a_brakes, v_brakes, p_brakes;

    // For velocity interface
    template<JerkSigns jerk_signs, Limits limits>
    bool check(double jf, double aMax, double aMin) {
        if (t[0] < 0) {
            return false;
        }

        t_sum[0] = t[0];
        for (size_t i = 0; i < 6; ++i) {
            if (t[i+1] < 0) {
                return false;
            }

            t_sum[i+1] = t_sum[i] + t[i+1];
        }

        if constexpr (limits == Limits::ACC0) {
            if (t[1] < std::numeric_limits<double>::epsilon()) {
                return false;
            }
        }

        if (t_sum[6] > 1e12) { // For numerical reasons, is that needed?
            return false;
        }

        if constexpr (jerk_signs == JerkSigns::UDDU) {
            j = {jf, 0, -jf, 0, -jf, 0, jf};
        } else {
            j = {jf, 0, -jf, 0, jf, 0, -jf};
        }

        for (size_t i = 0; i < 7; ++i) {
            a[i+1] = a[i] + t[i] * j[i];
            v[i+1] = v[i] + t[i] * (a[i] + t[i] * j[i] / 2);
            p[i+1] = p[i] + t[i] * (v[i] + t[i] * (a[i] / 2 + t[i] * j[i] / 6));
        }

        this->jerk_signs = jerk_signs;
        this->limits = limits;

        const double aUppLim = ((aMax > 0) ? aMax : aMin) + 1e-12;
        const double aLowLim = ((aMax > 0) ? aMin : aMax) - 1e-12;

        // Velocity limit can be broken in the beginning if both initial velocity and acceleration are too high
        // std::cout << std::setprecision(15) << "target: " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << " T: " << t_sum[6] << " " << to_string() << std::endl;
        return std::abs(v[7] - vf) < 1e-8
            && std::abs(a[7] - af) < 1e-12 // This is not really needed, but we want to double check
            && a[1] >= aLowLim && a[3] >= aLowLim && a[5] >= aLowLim
            && a[1] <= aUppLim && a[3] <= aUppLim && a[5] <= aUppLim;
    }

    // For position interface
    template<JerkSigns jerk_signs, Limits limits>
    bool check(double jf, double vMax, double vMin, double aMax, double aMin) {
        if (t[0] < 0) {
            return false;
        }

        t_sum[0] = t[0];
        for (size_t i = 0; i < 6; ++i) {
            if (t[i+1] < 0) {
                return false;
            }

            t_sum[i+1] = t_sum[i] + t[i+1];
        }

        if constexpr (limits == Limits::ACC0_ACC1_VEL || limits == Limits::ACC0_VEL || limits == Limits::ACC1_VEL || limits == Limits::VEL) {
            if (t[3] < std::numeric_limits<double>::epsilon()) {
                return false;
            }
        }

        if constexpr (limits == Limits::ACC0 || limits == Limits::ACC0_ACC1) {
            if (t[1] < std::numeric_limits<double>::epsilon()) {
                return false;
            }
        }

        if constexpr (limits == Limits::ACC1 || limits == Limits::ACC0_ACC1) {
            if (t[5] < std::numeric_limits<double>::epsilon()) {
                return false;
            }
        }

        if (t_sum[6] > 1e12) { // For numerical reasons, is that needed?
            return false;
        }

        if constexpr (jerk_signs == JerkSigns::UDDU) {
            j = {jf, 0, -jf, 0, -jf, 0, jf};
        } else {
            j = {jf, 0, -jf, 0, jf, 0, -jf};
        }

        const double vUppLim = ((vMax > 0) ? vMax : vMin) + 1e-12;
        const double vLowLim = ((vMax > 0) ? vMin : vMax) - 1e-12;

        for (size_t i = 0; i < 7; ++i) {
            a[i+1] = a[i] + t[i] * j[i];
            v[i+1] = v[i] + t[i] * (a[i] + t[i] * j[i] / 2);
            p[i+1] = p[i] + t[i] * (v[i] + t[i] * (a[i] / 2 + t[i] * j[i] / 6));

            if constexpr (limits == Limits::ACC0_ACC1_VEL || limits == Limits::ACC0_VEL || limits == Limits::ACC1_VEL || limits == Limits::VEL) {
                if (i == 2) {
                    a[3] = 0.0;
                }
            }

            if (i > 1 && a[i+1] * a[i] < -std::numeric_limits<double>::epsilon()) {
                const double v_a_zero = v[i] - (a[i] * a[i]) / (2 * j[i]);
                if (v_a_zero > vUppLim || v_a_zero < vLowLim) {
                    return false;
                }
            }
        }

        this->jerk_signs = jerk_signs;
        this->limits = limits;

        const double aUppLim = ((aMax > 0) ? aMax : aMin) + 1e-12;
        const double aLowLim = ((aMax > 0) ? aMin : aMax) - 1e-12;

        // Velocity limit can be broken in the beginning if both initial velocity and acceleration are too high
        // std::cout << std::setprecision(16) << "target: " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << " T: " << t_sum[6] << " " << to_string() << std::endl;
        return std::abs(p[7] - pf) < 1e-8
            && std::abs(v[7] - vf) < 1e-8
            && std::abs(a[7] - af) < 1e-12 // This is not really needed, but we want to double check
            && a[1] >= aLowLim && a[3] >= aLowLim && a[5] >= aLowLim
            && a[1] <= aUppLim && a[3] <= aUppLim && a[5] <= aUppLim
            && v[3] <= vUppLim && v[4] <= vUppLim && v[5] <= vUppLim && v[6] <= vUppLim
            && v[3] >= vLowLim && v[4] >= vLowLim && v[5] >= vLowLim && v[6] >= vLowLim; // This is not really needed, but we want to double check
    }

    template<JerkSigns jerk_signs, Limits limits>
    inline bool check(double, double jf, double vMax, double vMin, double aMax, double aMin) {
        // Time doesn't need to be checked as every profile has a: tf - ... equation
        return check<jerk_signs, limits>(jf, vMax, vMin, aMax, aMin); // && (std::abs(t_sum[6] - tf) < 1e-8);
    }

    template<JerkSigns jerk_signs, Limits limits>
    inline bool check(double tf, double jf, double vMax, double vMin, double aMax, double aMin, double jMax) {
        return (std::abs(jf) < std::abs(jMax) + 1e-12) && check<jerk_signs, limits>(tf, jf, vMax, vMin, aMax, aMin);
    }

    //! Integrate with constant jerk for duration t. Returns new position, new velocity, and new acceleration.
    inline static std::tuple<double, double, double> integrate(double t, double p0, double v0, double a0, double j) {
        return std::make_tuple(
            p0 + t * (v0 + t * (a0 / 2 + t * j / 6)),
            v0 + t * (a0 + t * j / 2),
            a0 + t * j
        );
    }

    //! Set boundary values for the position interface
    inline void set_boundary(double p0_new, double v0_new, double a0_new, double pf_new, double vf_new, double af_new) {
        a[0] = a0_new;
        v[0] = v0_new;
        p[0] = p0_new;
        af = af_new;
        vf = vf_new;
        pf = pf_new;
    }

    //! Set boundary values for the velocity interface
    inline void set_boundary(double p0_new, double v0_new, double a0_new, double vf_new, double af_new) {
        a[0] = a0_new;
        v[0] = v0_new;
        p[0] = p0_new;
        af = af_new;
        vf = vf_new;
    }

    inline static void check_position_extremum(double t_ext, double t_sum, double t, double p, double v, double a, double j, PositionExtrema& ext) {
        if (0 < t_ext && t_ext < t) {
            double p_ext, a_ext;
            std::tie(p_ext, std::ignore, a_ext) = integrate(t_ext, p, v, a, j);
            if (a_ext > 0 && p_ext < ext.min) {
                ext.min = p_ext;
                ext.t_min = t_sum + t_ext;
            } else if (a_ext < 0 && p_ext > ext.max) {
                ext.max = p_ext;
                ext.t_max = t_sum + t_ext;
            }
        }
    }

    static void check_step_for_position_extremum(double t_sum, double t, double p, double v, double a, double j, PositionExtrema& ext) {
        if (p < ext.min) {
            ext.min = p;
            ext.t_min = t_sum;
        }
        if (p > ext.max) {
            ext.max = p;
            ext.t_max = t_sum;
        }

        if (j != 0) {
            const double D = a * a - 2 * j * v;
            if (std::abs(D) < std::numeric_limits<double>::epsilon()) {
                check_position_extremum(-a / j, t_sum, t, p, v, a, j, ext);

            } else if (D > 0) {
                const double D_sqrt = std::sqrt(D);
                check_position_extremum((-a - D_sqrt) / j, t_sum, t, p, v, a, j, ext);
                check_position_extremum((-a + D_sqrt) / j, t_sum, t, p, v, a, j, ext);
            }
        }
    }

    PositionExtrema get_position_extrema() const {
        PositionExtrema extrema;
        extrema.min = std::numeric_limits<double>::infinity();
        extrema.max = -std::numeric_limits<double>::infinity();

        if (t_brake) {
            if (t_brakes[0] > 0.0) {
                check_step_for_position_extremum(0.0, t_brakes[0], p_brakes[0], v_brakes[0], a_brakes[0], j_brakes[0], extrema);

                if (t_brakes[1] > 0.0) {
                    check_step_for_position_extremum(t_brakes[0], t_brakes[1], p_brakes[1], v_brakes[1], a_brakes[1], j_brakes[1], extrema);
                }
            }
        }

        double t_current_sum {0.0};
        for (size_t i = 0; i < 7; ++i) {
            if (i > 0) {
                t_current_sum = t_sum[i - 1];
            }
            check_step_for_position_extremum(t_current_sum + t_brake.value_or(0.0), t[i], p[i], v[i], a[i], j[i], extrema);
        }

        if (pf < extrema.min) {
            extrema.min = pf;
            extrema.t_min = t_sum[6] + t_brake.value_or(0.0);
        }
        if (pf > extrema.max) {
            extrema.max = pf;
            extrema.t_max = t_sum[6] + t_brake.value_or(0.0);
        }

        return extrema;
    }

    bool get_first_state_at_position(double pt, double& time, double& vt, double& at, double offset = 0.0) const {
        for (size_t i = 0; i < 7; ++i) {
            if (std::abs(p[i] - pt) < std::numeric_limits<double>::epsilon()) {
                time = offset + ((i > 0) ? t_sum[i-1] : 0.0);
                vt = v[i];
                at = a[i];
                return true;
            }

            if (t[i] == 0.0) {
                continue;
            }

            for (const double _t: Roots::solveCub(j[i]/6, a[i]/2, v[i], p[i]-pt)) {
                if (0 < _t && _t <= t[i]) {
                    time = offset + _t + ((i > 0) ? t_sum[i-1] : 0.0);
                    std::tie(std::ignore, vt, at) = integrate(_t, p[i], v[i], a[i], j[i]);
                    return true;
                }
            }
        }

        if (std::abs(pf - pt) < 1e-9) {
            time = offset + t_sum[6];
            vt = vf;
            at = af;
            return true;
        }

        return false;
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
        switch (jerk_signs) {
            case JerkSigns::UDDU: result += "_UDDU"; break;
            case JerkSigns::UDUD: result += "_UDUD"; break;
        }
        return result;
    }
};

} // namespace ruckig
