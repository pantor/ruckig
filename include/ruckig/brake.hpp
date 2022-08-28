#pragma once

#include <array>
#include <cmath>
#include <iostream>

#ifdef WITH_SERIALIZATION
#include <json/json.hpp>
#endif

#include <ruckig/utils.hpp>


namespace ruckig {

//! Calculates (pre- or post-) profile to get current or final state below the limits
class BrakeProfile {
    static constexpr double eps {2.2e-14};

    void acceleration_brake(double v0, double a0, double vMax, double vMin, double aMax, double aMin, double jMax);
    void velocity_brake(double v0, double a0, double vMax, double vMin, double aMax, double aMin, double jMax);

public:
    //! Overall duration
    double duration {0.0};

    //! Profile information for a two-step profile
    std::array<double, 2> t, j, a, v, p;

    //! Calculate brake trajectory for position interface
    void get_position_brake_trajectory(double v0, double a0, double vMax, double vMin, double aMax, double aMin, double jMax);

    //! Calculate brake trajectory for velocity interface
    void get_velocity_brake_trajectory(double a0, double aMax, double aMin, double jMax);

    //! Finalize by integrating along kinematic state
    void finalize(double& ps, double& vs, double& as) {
        if (t[0] <= 0.0 && t[1] <= 0.0) {
            duration = 0.0;
            return;
        }

        duration = t[0];
        p[0] = ps;
        v[0] = vs;
        a[0] = as;
        std::tie(ps, vs, as) = integrate(t[0], ps, vs, as, j[0]);

        if (t[1] > 0.0) {
            duration += t[1];
            p[1] = ps;
            v[1] = vs;
            a[1] = as;
            std::tie(ps, vs, as) = integrate(t[1], ps, vs, as, j[1]);
        }
    }

#ifdef WITH_SERIALIZATION
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(BrakeProfile, duration, t, j, a, v, p)
#endif
};

} // namespace ruckig
