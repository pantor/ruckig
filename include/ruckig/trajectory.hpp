#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <tuple>


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

    template<JerkSigns jerk_signs, Limits limits>
    bool check(double pf, double vf, double af, double jf, double vMax, double aMax) {
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

        const double vMaxAbs = std::abs(vMax) + 1e-12;
        const double aMaxAbs = std::abs(aMax) + 1e-12;

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
                if (std::abs(v_a_zero) > vMaxAbs) {
                    return false;
                }
            }
        }

        this->jerk_signs = jerk_signs;
        this->limits = limits;

        // Velocity limit can be broken in the beginning if both initial velocity and acceleration are too high
        // std::cout << std::setprecision(15) << "target: " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << " T: " << t_sum[6] << " " << to_string() << std::endl;
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
    
    template<JerkSigns jerk_signs, Limits limits>
    inline bool check([[maybe_unused]] double tf, double pf, double vf, double af, double jf, double vMax, double aMax) {
        // Time doesn't need to be checked as every profile has a: tf - ... equation
        return check<jerk_signs, limits>(pf, vf, af, jf, vMax, aMax); // && (std::abs(t_sum[6] - tf) < 1e-8);
    }
    
    template<JerkSigns jerk_signs, Limits limits>
    inline bool check(double tf, double pf, double vf, double af, double jf, double vMax, double aMax, double jMax) {
        return (std::abs(jf) < std::abs(jMax) + 1e-12) && check<jerk_signs, limits>(tf, pf, vf, af, jf, vMax, aMax);
    }

    //! Integrate with constant jerk for duration t. Returns new position, new velocity, and new acceleration.
    inline static std::tuple<double, double, double> integrate(double t, double p0, double v0, double a0, double j) {
        return {
            p0 + t * (v0 + t * (a0 / 2 + t * j / 6)),
            v0 + t * (a0 + t * j / 2),
            a0 + t * j,
        };
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
                ext.max = t_sum + t_ext;
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
        for (size_t i = 0; i < 7; ++i) {
            check_step_for_position_extremum(t_sum[i] + t_brake.value_or(0.0), t[i], p[i], v[i], a[i], j[i], extrema);
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


template<size_t DOFs>
struct Trajectory {
    //! Duration of the synchronized trajectory
    double duration;

    //! Set of current profiles for each DoF
    std::array<Profile, DOFs> profiles;

    //! Minimum duration of each independent DoF
    std::array<double, DOFs> independent_min_durations;

    //! Get the output parameter for the given time
    void at_time(double time, std::array<double, DOFs>& new_position, std::array<double, DOFs>& new_velocity, std::array<double, DOFs>& new_acceleration) const {
        if (time > duration) {
            // Keep constant acceleration
            for (size_t dof = 0; dof < DOFs; ++dof) {
                std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(time - duration, profiles[dof].pf, profiles[dof].vf, profiles[dof].af, 0);
            }
            return;
        }

        for (size_t dof = 0; dof < DOFs; ++dof) {
            const Profile& p = profiles[dof];

            double t_diff = time;
            if (p.t_brake) {
                if (t_diff < p.t_brake.value()) {
                    const size_t index = (t_diff < p.t_brakes[0]) ? 0 : 1;
                    if (index > 0) {
                        t_diff -= p.t_brakes[index - 1];
                    }

                    std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff, p.p_brakes[index], p.v_brakes[index], p.a_brakes[index], p.j_brakes[index]);
                    continue;
                } else {
                    t_diff -= p.t_brake.value();
                }
            }

            // Non-time synchronization
            if (t_diff >= p.t_sum[6]) {
                // Keep constant acceleration
                std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff - p.t_sum[6], p.pf, p.vf, p.af, 0);
                continue;
            }

            const auto index_ptr = std::upper_bound(p.t_sum.begin(), p.t_sum.end(), t_diff);
            const size_t index = std::distance(p.t_sum.begin(), index_ptr);

            if (index > 0) {
                t_diff -= p.t_sum[index - 1];
            }

            std::tie(new_position[dof], new_velocity[dof], new_acceleration[dof]) = Profile::integrate(t_diff, p.p[index], p.v[index], p.a[index], p.j[index]);
        }
    }

    //! Get the min/max values of the position for each DoF and the current trajectory
    std::array<PositionExtrema, DOFs> get_position_extrema() {
        std::array<PositionExtrema, DOFs> result;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            result[dof] = profiles[dof].get_position_extrema();
        }
        return result;
    }
};

} // namespace ruckig
