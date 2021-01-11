#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>


namespace ruckig {

Step1::Step1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax): p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af) {
    // max values needs to be invariant to plus minus sign change
    pd = pf - p0;

    v0_v0 = v0 * v0;
    vf_vf = vf * vf;
    vMax_vMax = vMax * vMax;

    a0_a0 = a0 * a0;
    af_af = af * af;
    aMax_aMax = aMax * aMax;

    a0_p3 = a0 * a0 * a0;
    a0_p4 = a0_a0 * a0_a0;
    af_p3 = af * af * af;
    af_p4 = af_af * af_af;

    jMax_jMax = jMax * jMax;
}

void Step1::add_profile(Profile profile, Profile::Limits limits, double jMax) {
    profile.limits = limits;
    profile.direction = (jMax > 0) ? Profile::Direction::UP : Profile::Direction::DOWN;
    valid_profiles[valid_profile_counter] = profile;
    ++valid_profile_counter;
}

void Step1::time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax_aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*(a0_p4 + af_p4) + 8*aMax*(af_p3 - a0_p3) + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 6*af_af*(aMax_aMax - 2*jMax*vf) - 12*jMax*(aMax_aMax*(v0 + vf + 2*vMax) - jMax*(2*aMax*pd + v0_v0 + vf_vf - 2*vMax_vMax)))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = profile.t[2];
    profile.t[5] = (af_af/2 - aMax_aMax - jMax*(vf - vMax))/(aMax*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, Limits::ACC0_ACC1_VEL, jMax);
    }
}

void Step1::time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(a0_a0/2 + jMax*(vMax - v0))/Abs(jMax);

    profile.t[0] = (-a0 + h1*jMax)/jMax;
    profile.t[1] = 0;
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*af_p4 + 8*aMax*(af_p3 - a0_p3) + 24*aMax*jMax*(a0*v0 - af*vf) + 6*af_af*(aMax_aMax - 2*jMax*vf) - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(vf + vMax) + jMax*(-vf_vf + vMax_vMax) - aMax*h1*(a0_a0 - 2*jMax*(v0 + vMax))))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = aMax/jMax;
    profile.t[5] = (af_af/2 - aMax_aMax + jMax*(vMax - vf))/(aMax*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, Limits::ACC1_VEL, jMax);
    }
}

void Step1::time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(af_af/2 + jMax*(vMax - vf))/Abs(jMax);

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax_aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*a0_p4 + 8*(af_p3 - a0_p3)*aMax + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 12*af_af*aMax*h1*jMax - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vMax) + jMax*(-v0_v0 + vMax_vMax) + 2*aMax*(vf + vMax)*h1*jMax))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, Limits::ACC0_VEL, jMax);
    }
}

void Step1::time_up_vel(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(af_af/2 + jMax*(vMax - vf))/Abs(jMax);
    const double h2 = Sqrt(a0_a0/2 + jMax*(vMax - v0))/Abs(jMax);

    // Solution 3/4
    profile.t[0] = (-a0 + h2*jMax)/jMax;
    profile.t[1] = 0;
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (2*(af_p3 - a0_p3) + 6*jMax*(a0*v0 - af*vf) + 3*af_af*h1*jMax + 3*jMax_jMax*(h2*(a0_a0/jMax - 2*(v0 + vMax)) + 2*pd - 2*(vf + vMax)*h1))/(6*jMax_jMax*vMax);
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, Limits::VEL, jMax);
    }
}

void Step1::time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt((a0_p4 + af_p4)/2 + 4./3*(af_p3 - a0_p3)*aMax + 4*aMax*jMax*(a0*v0 - af*vf) + a0_a0*(aMax_aMax - 2*jMax*v0) + af_af*(aMax_aMax - 2*jMax*vf) + (aMax_aMax*(aMax_aMax - 2*jMax*(v0 + vf)) + 2*jMax_jMax*(v0_v0 + vf_vf + 2*aMax*pd)));
    const double h2 = a0_a0 - 3*aMax_aMax - 2*jMax*v0;

    // UDDU: Solution 2
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (h2 - h1)/(2*aMax*jMax);
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = profile.t[2];
        profile.t[5] = profile.t[1] + ((af_af - a0_a0)/2 - jMax*(vf - v0))/(aMax*jMax);
        profile.t[6] = profile.t[4] + af/jMax;

        if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, Limits::ACC0_ACC1, jMax);
        }
    }

    // UDDU: Solution 1
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (h2 + h1)/(2*aMax*jMax);
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = profile.t[2];
        profile.t[5] = profile.t[1] + ((af_af - a0_a0)/2 - jMax*(vf - v0))/(aMax*jMax);
        profile.t[6] = profile.t[4] + af/jMax;

        if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, Limits::ACC0_ACC1, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        const double h3 = (3*(a0_p4 - af_p4) + 8*aMax*(af_p3 - a0_p3) + 24*aMax*jMax*(a0*v0 + af*vf) + 12*jMax_jMax*(2*aMax*pd - (vf_vf - v0_v0)) - 12*jMax*(a0_a0*v0 + af_af*vf))/(12*aMax_aMax) + (a0_a0 - af_af)/2  - jMax*(v0 + vf);

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (h3 + a0_a0 - 2*aMax_aMax - 2*jMax*v0)/(2*aMax*jMax);
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = profile.t[2];
        profile.t[5] = -(h3 - af_af + 2*aMax_aMax - 2*jMax*vf)/(2*aMax*jMax);
        profile.t[6] = profile.t[4] - af/jMax;

        if (profile.check<Teeth::UDUD>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, Limits::ACC0_ACC1, jMax);
        }
    }
}

void Step1::time_up_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = 2*(2*a0 + aMax)/jMax;
    polynom[2] = (5*a0_a0 + 6*a0*aMax + aMax_aMax + 2*jMax*v0)/jMax_jMax;
    polynom[3] = (2*(a0 + aMax)*(a0_a0 + a0*aMax + 2*jMax*v0))/Power(jMax,3);
    polynom[4] = (3*(a0_p4 - af_p4) + 8*(a0_p3 - af_p3)*aMax + 24*aMax*jMax*(a0*v0 + af*vf) + 6*a0_a0*(aMax_aMax + 2*jMax*v0) - 6*af_af*(aMax_aMax - 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf) + jMax*(v0_v0 - vf_vf)))/(12*jMax_jMax*jMax_jMax);
    
    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        profile.t[0] = t;
        profile.t[1] = 0;
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = aMax/jMax;
        profile.t[5] = ((a0_a0 + af_af)/2 - aMax_aMax + 2*a0*jMax*t + jMax_jMax*t*t - jMax*(vf - v0))/(aMax*jMax);
        profile.t[6] = profile.t[4] + af/jMax;

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, Limits::ACC1, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = 2*(2*a0 - aMax)/jMax;
        polynom[2] = (5*a0_a0 - 6*a0*aMax - aMax_aMax + 2*jMax*v0)/jMax_jMax;
        polynom[3] = (2*(a0_p3 - a0*aMax*(2*a0 + aMax) + 2*(a0 - aMax)*jMax*v0))/Power(jMax,3);
        polynom[4] = -(3*(af_p4 - a0_p4) + 8*(a0_p3 - af_p3)*aMax + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 6*af_af*(aMax_aMax + 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf) + jMax*(vf_vf - v0_v0)))/(12*jMax_jMax*jMax_jMax);

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0) {
                continue;
            }

            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = 0;
            profile.t[4] = aMax/jMax;
            profile.t[5] = -((a0_a0 - af_af)/2 + aMax_aMax + 2*a0*jMax*t+ jMax*(jMax*t*t - (vf - v0)))/(aMax*jMax);
            profile.t[6] = profile.t[4] - af/jMax;

            if (profile.check<Teeth::UDUD>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, Limits::ACC1, jMax);
            }
        }
    }
}

void Step1::time_up_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (2*aMax)/jMax;
    polynom[2] = (-af_af + aMax_aMax + 2*jMax*vf)/jMax_jMax;
    polynom[3] = (-2*aMax*(af_af - 2*jMax*vf))/Power(jMax,3);
    polynom[4] = (3*(af_p4 - a0_p4) + 8*(a0_p3 - af_p3)*aMax - 24*aMax*jMax*(a0*v0 - af*vf) - 6*a0_a0*(aMax_aMax - 2*jMax*v0) - 6*af_af*(aMax_aMax + 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf) + jMax*(vf_vf - v0_v0)))/(12*jMax_jMax*jMax_jMax);

    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (a0_a0/2 - af_af/2 - aMax_aMax + jMax_jMax*t*t + jMax*(vf - v0))/(aMax*jMax);
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = t;
        profile.t[5] = 0;
        profile.t[6] = profile.t[4] + af/jMax;

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, Limits::ACC0, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (-2*aMax)/jMax;
        polynom[2] = -((af_af + aMax_aMax + 2*jMax*vf)/jMax_jMax);
        polynom[3] = (2*aMax*(af_af + 2*jMax*vf))/Power(jMax,3);
        polynom[4] = (3*(af_p4 - a0_p4) + 8*(a0_p3 - af_p3)*aMax - 24*aMax*jMax*(a0*v0 + af*vf) - 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 6*af_af*(aMax_aMax + 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf) + jMax*(vf_vf - v0_v0)))/(12*jMax_jMax*jMax_jMax);

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0) {
                continue;
            }

            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (a0_a0 + af_af - 2*(aMax_aMax + jMax*(jMax*t*t - (vf - v0))))/(2*aMax*jMax);
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = 0;
            profile.t[4] = t;
            profile.t[5] = 0;
            profile.t[6] = profile.t[4] - af/jMax;
            
            if (profile.check<Teeth::UDUD>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, Limits::ACC0, jMax);
            }
        }
    }
}

void Step1::time_up_none(Profile& profile, double vMax, double aMax, double jMax) {
    if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        profile.t[0] = std::cbrt(pd/(2*jMax));
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, Limits::NONE, jMax);
        }
        return;
    }

    const double p1 = af_af - a0_a0 - 2*jMax*(vf - v0);
    const double p2 = af_p4 - 4*af_af*jMax*(vf - v0) - 4*jMax_jMax*(v0_v0 + 2*v0*vf - vf_vf);
    const double p3 = af_p3 + 3*jMax_jMax*pd - 3*af*jMax*vf;

    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (-4*(2*a0_p3 + af_p3 + 3*jMax_jMax*pd - 3*a0*(af_af + jMax*(v0 - 2*vf)) - 3*af*jMax*vf))/(3*jMax*p1);
    polynom[2] = (-3*a0_p4 + Power(af_af - 2*jMax*(vf - v0),2) - 8*a0*p3 + 2*a0_a0*(5*af_af - 2*jMax*(v0 + 5*vf)))/(2*jMax_jMax*p1);
    polynom[3] = -(Power(a0,5) + 8*(a0_a0 + jMax*v0)*p3 + a0_p3*(-6*af_af + 4*jMax*(v0 + 3*vf)) - 3*a0*p2)/(3*Power(jMax,3)*p1);
    polynom[4] = -(Power(a0,6) + Power(af,6) - 48*af_p3*jMax_jMax*pd + 144*af*Power(jMax,3)*pd*vf - 6*af_p4*jMax*(3*v0 + vf) + 16*a0_p3*p3 + 48*a0*jMax*v0*p3 - 36*af_af*jMax_jMax*(v0_v0 - 2*v0*vf - vf_vf) - 72*Power(jMax,3)*(jMax*Power(pd,2) + (vf - v0)*Power(v0 + vf,2)) + a0_p4*(-9*af_af + 6*jMax*(v0 + 3*vf)) - 9*a0_a0*p2)/(36*jMax_jMax*jMax_jMax*p1);

    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        // Max two Newton steps
        if (std::abs(Roots::polyEval(polynom, t)) > 1e-14) {
            t -= Roots::polyEval(polynom, t) / Roots::polyEval(Roots::polyDeri(polynom), t);

            if (std::abs(Roots::polyEval(polynom, t)) > 1e-14) {
                t -= Roots::polyEval(polynom, t) / Roots::polyEval(Roots::polyDeri(polynom), t);
            }
        }

        profile.t[0] = t;
        profile.t[1] = 0;
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = Sqrt((a0_a0 + af_af)/2 + 2*a0*jMax*t + jMax*(jMax*t*t - (vf - v0)))/Abs(jMax);
        profile.t[5] = 0;
        profile.t[6] = profile.t[4] + af/jMax;
        
        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        if (profile.check<Teeth::UDDU>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, Limits::NONE, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        const double p1 = af_p3 + 3*jMax_jMax*pd + 3*af*jMax*vf;
        const double p2 = af_p4 + 4*af_af*jMax*(vf - v0) + 4*jMax_jMax*(-5*v0_v0 - 2*v0*vf + vf_vf);
        const double p3 = 3*af_af + 34*jMax*v0 + 6*jMax*vf;

        std::array<double, 7> polynom;
        polynom[0] = 1.0;
        polynom[1] = (6*a0)/jMax;
        polynom[2] = (53*a0_a0 + af_af + 2*jMax*(7*v0 + vf))/(4*jMax_jMax);
        polynom[3] = (40*a0_p3 - af_p3 - 3*jMax_jMax*pd - 3*af*jMax*vf + 3*a0*(af_af + 13*jMax*v0 + 2*jMax*vf))/(3*Power(jMax,3));
        polynom[4] = (51*a0_p4 - af_p4 - 4*af_af*jMax*(vf - v0) + 2*a0_a0*(5*af_af + 58*jMax*v0 + 10*jMax*vf) - 8*a0*p1 + 4*jMax_jMax*(7*v0_v0 + 2*v0*vf - vf_vf))/(8*jMax_jMax*jMax_jMax);
        polynom[5] = (17*Power(a0,5) + 2*a0_p3*p3 - 8*(a0_a0 + jMax*v0)*p1 - 3*a0*p2)/(12*Power(jMax,5));
        polynom[6] = -(-17*Power(a0,6) + Power(af,6) - 48*af_p3*jMax_jMax*pd - 144*af*Power(jMax,3)*pd*vf + 6*af_p4*jMax*(3*v0 + vf) - 3*a0_p4*p3 + 16*a0_p3*p1 + 48*a0*jMax*v0*p1 - 36*af_af*jMax_jMax*(v0_v0 - 2*v0*vf - vf_vf) - 72*Power(jMax,3)*(jMax*Power(pd,2) - (vf - v0)*Power(v0 + vf,2)) + 9*a0_a0*p2)/(144*Power(jMax,6));

        std::array<double, 6> deriv;
        deriv[0] = 1.0;
        deriv[1] = 5./6 * polynom[1];
        deriv[2] = 4./6 * polynom[2];
        deriv[3] = 3./6 * polynom[3];
        deriv[4] = 2./6 * polynom[4];
        deriv[5] = 1./6 * polynom[5];

        auto dd_extremas = Roots::solveQuartMonic(4./5 * deriv[1], 3./5 * deriv[2], 2./5 * deriv[3], 1./5 * deriv[4]);
        std::set<std::tuple<double, double>> dd_tz_intervals;

        double tz_min {0.0};
        double tz_max = {1e5};
        double dd_tz_current {tz_min};

        for (double tz: dd_extremas) {
            if (tz <= 0.0 || tz >= tz_max) {
                continue;
            }

            // Check that polynom(lower) and polynom(upper) have different signs (should only happen at first and last boundary)
            if (Roots::polyEval(deriv, dd_tz_current) * Roots::polyEval(deriv, tz) < 0) {
                dd_tz_intervals.insert({dd_tz_current, tz});
            }
            dd_tz_current = tz;
        }
        if (Roots::polyEval(deriv, dd_tz_current) * Roots::polyEval(deriv, tz_max) < 0) {
            dd_tz_intervals.insert({dd_tz_current, tz_max});
        }

        std::set<std::tuple<double, double>> tz_intervals;
        double tz_current {tz_min};

        for (auto interval: dd_tz_intervals) {
            double lower = std::get<0>(interval);
            double upper = std::get<1>(interval);
            double tz = Roots::shrinkInterval(deriv, lower, upper, 1e-14);

            if (tz <= 0.0) {
                continue;
            }
            // Check that polynom(lower) and polynom(upper) have different signs (should only happen at first and last boundary)
            if (Roots::polyEval(polynom, tz_current) * Roots::polyEval(polynom, tz) < 0) {
                tz_intervals.insert({tz_current, tz});
            }
            tz_current = tz;
        }
        if (Roots::polyEval(polynom, tz_current) * Roots::polyEval(polynom, tz_max) < 0) {
            tz_intervals.insert({tz_current, tz_max});
        }

        for (auto interval: tz_intervals) {
            // Use safe Newton method
            double lower = std::get<0>(interval);
            double upper = std::get<1>(interval);
            double t = Roots::shrinkInterval(polynom, lower, upper, 1e-14);

            if (t < 0.0) {
                continue;
            }

            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = 0;
            profile.t[4] = Sqrt((af_af - a0_a0)/2 - 2*a0*jMax*t - jMax*(jMax*t*t - (vf - v0)))/Abs(jMax);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4] - af/jMax;

            if (profile.check<Teeth::UDUD>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, Limits::NONE, jMax);
            }
        }
    }
}

void Step1::time_down_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, -vMax, -aMax, -jMax);
}

void Step1::time_down_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc1_vel(profile, -vMax, -aMax, -jMax);
}

void Step1::time_down_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_vel(profile, -vMax, -aMax, -jMax);
}

void Step1::time_down_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_vel(profile, -vMax, -aMax, -jMax);
}

void Step1::time_down_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, -vMax, -aMax, -jMax);
}

void Step1::time_down_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc1(profile, -vMax, -aMax, -jMax);
}

void Step1::time_down_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0(profile, -vMax, -aMax, -jMax);
}

void Step1::time_down_none(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_none(profile, -vMax, -aMax, -jMax);
}

bool Step1::get_profile(const Profile& input, double vMax, double aMax, double jMax) {
    Profile profile = input;
    profile.a[0] = a0;
    profile.v[0] = v0;
    profile.p[0] = p0;
    valid_profile_counter = 0;

    if (pf > p0) {
        time_up_acc0_acc1_vel(profile, vMax, aMax, jMax);
        time_down_acc0_acc1_vel(profile, vMax, aMax, jMax);
        time_up_acc1_vel(profile, vMax, aMax, jMax);
        time_down_acc1_vel(profile, vMax, aMax, jMax);
        time_up_acc0_vel(profile, vMax, aMax, jMax);
        time_down_acc0_vel(profile, vMax, aMax, jMax);
        time_up_vel(profile, vMax, aMax, jMax);
        time_down_vel(profile, vMax, aMax, jMax);
        time_up_none(profile, vMax, aMax, jMax);
        time_up_acc0(profile, vMax, aMax, jMax);
        time_up_acc1(profile, vMax, aMax, jMax);
        time_up_acc0_acc1(profile, vMax, aMax, jMax);
        time_down_none(profile, vMax, aMax, jMax);
        time_down_acc0(profile, vMax, aMax, jMax);
        time_down_acc1(profile, vMax, aMax, jMax);
        time_down_acc0_acc1(profile, vMax, aMax, jMax);

    } else {
        time_down_acc0_acc1_vel(profile, vMax, aMax, jMax);
        time_up_acc0_acc1_vel(profile, vMax, aMax, jMax);
        time_down_acc1_vel(profile, vMax, aMax, jMax);
        time_up_acc1_vel(profile, vMax, aMax, jMax);
        time_down_acc0_vel(profile, vMax, aMax, jMax);
        time_up_acc0_vel(profile, vMax, aMax, jMax);
        time_down_vel(profile, vMax, aMax, jMax);
        time_up_vel(profile, vMax, aMax, jMax);
        time_down_none(profile, vMax, aMax, jMax);
        time_down_acc0(profile, vMax, aMax, jMax);
        time_down_acc1(profile, vMax, aMax, jMax);
        time_down_acc0_acc1(profile, vMax, aMax, jMax);
        time_up_none(profile, vMax, aMax, jMax);
        time_up_acc0(profile, vMax, aMax, jMax);
        time_up_acc1(profile, vMax, aMax, jMax);
        time_up_acc0_acc1(profile, vMax, aMax, jMax);
    }

    if (!valid_profile_counter) {
        // std::cout << "HERE" << std::endl;
        return false;
    }

    // auto min_profile = std::min_element(valid_profiles.begin(), valid_profiles.end(), [](const Profile& a, const Profile& b){ return a.t_sum[6] < b.t_sum[6]; }); 

    double t_brake = valid_profiles[0].t_brake.value_or(0.0);
    block = Block {valid_profiles[0].t_sum[6] + t_brake, valid_profiles[0]};

    if (valid_profile_counter == 2) {
        {
            auto& p_left = valid_profiles[0];
            auto& p_right = valid_profiles[1];

            if (p_left.t_sum[6] < p_right.t_sum[6]) {
                block.a = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_a = p_right;
            }
        }
        
    } else if (valid_profile_counter == 3) {
        if (valid_profiles[2].t_sum[6] < valid_profiles[0].t_sum[6]) {
            block = Block {valid_profiles[2].t_sum[6] + t_brake, valid_profiles[2]};
            
            auto& p_left = valid_profiles[0];
            auto& p_right = valid_profiles[1];
            if (p_left.direction == p_right.direction && p_left.t_sum[6] < p_right.t_sum[6]) {
                block.a = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_a = p_right;
            }

        } else {
            auto& p_left = valid_profiles[1];
            auto& p_right = valid_profiles[2];
            if (p_left.direction == p_right.direction && p_left.t_sum[6] < p_right.t_sum[6]) {
                block.a = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_a = p_right;

            } else if (p_left.direction == p_right.direction && p_right.t_sum[6] < p_left.t_sum[6]) {
                block.a = Block::Interval {p_right.t_sum[6] + t_brake, p_left.t_sum[6] + t_brake};
                block.p_a = p_left;
            }
        }

    } else if (valid_profile_counter == 4) {
        if (valid_profiles[0].direction == valid_profiles[1].direction) {
            auto& p_left = valid_profiles[0];
            auto& p_right = valid_profiles[1];
            if (p_left.direction == p_right.direction && p_left.t_sum[6] < p_right.t_sum[6]) {
                block.a = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_a = p_right;

            } else if (p_left.direction == p_right.direction && p_right.t_sum[6] < p_left.t_sum[6]) {
                block.a = Block::Interval {p_right.t_sum[6] + t_brake, p_left.t_sum[6] + t_brake};
                block.p_a = p_left;
            }

        } else if (valid_profiles[1].direction == valid_profiles[2].direction) {
            auto& p_left = valid_profiles[1];
            auto& p_right = valid_profiles[2];
            if (p_left.direction == p_right.direction && p_left.t_sum[6] < p_right.t_sum[6]) {
                block.a = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_a = p_right;

            } else if (p_left.direction == p_right.direction && p_right.t_sum[6] < p_left.t_sum[6]) {
                block.a = Block::Interval {p_right.t_sum[6] + t_brake, p_left.t_sum[6] + t_brake};
                block.p_a = p_left;
            }
        }

        {
            auto& p_left = valid_profiles[2];
            auto& p_right = valid_profiles[3];
            if (p_left.direction == p_right.direction && p_left.t_sum[6] < p_right.t_sum[6]) {
                block.b = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_b = p_right;
            } else if (p_left.direction == p_right.direction && p_left.t_sum[6] > p_right.t_sum[6]) {
                block.b = Block::Interval {p_right.t_sum[6] + t_brake, p_left.t_sum[6] + t_brake};
                block.p_b = p_left;
            }
        }
        
    } else if (valid_profile_counter == 5) {
        {
            auto& p_left = valid_profiles[1];
            auto& p_right = valid_profiles[2];
            if (p_left.direction == p_right.direction && p_left.t_sum[6] < p_right.t_sum[6]) {
                block.a = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_a = p_right;
            } else if (p_left.direction == p_right.direction && p_left.t_sum[6] > p_right.t_sum[6]) {
                block.a = Block::Interval {p_right.t_sum[6] + t_brake, p_left.t_sum[6] + t_brake};
                block.p_a = p_left;
            }
        }

        {
            auto& p_left = valid_profiles[3];
            auto& p_right = valid_profiles[4];
            if (p_left.direction == p_right.direction && p_left.t_sum[6] < p_right.t_sum[6]) {
                block.b = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_b = p_right;
            } else if (p_left.direction == p_right.direction && p_left.t_sum[6] > p_right.t_sum[6]) {
                block.b = Block::Interval {p_right.t_sum[6] + t_brake, p_left.t_sum[6] + t_brake};
                block.p_b = p_left;
            }
        }
        
    } else if (valid_profile_counter == 5) {
        
    } else if (valid_profile_counter == 6) {
        {
            auto& p_left = valid_profiles[4];
            auto& p_right = valid_profiles[5];
            if (p_left.direction == p_right.direction && p_left.t_sum[6] < p_right.t_sum[6]) {
                block.b = Block::Interval {p_left.t_sum[6] + t_brake, p_right.t_sum[6] + t_brake};
                block.p_b = p_right;
            } else if (p_left.direction == p_right.direction && p_left.t_sum[6] > p_right.t_sum[6]) {
                block.b = Block::Interval {p_right.t_sum[6] + t_brake, p_left.t_sum[6] + t_brake};
                block.p_b = p_left;
            }
        }
        
    } else {
        // std::cerr << "MORE THAN 6 PROFILES " << valid_profile_counter << std::endl;
        // return false;
    }

    // std::cout << "---\n " << valid_profile_counter << std::endl;

    if (valid_profile_counter > 1) {
        std::cout << "---\n " << valid_profile_counter << std::endl;
        for (size_t i = 0; i < valid_profile_counter; ++i) {
            auto& p = valid_profiles[i];
            std::cout << p.t_sum[6] + t_brake << " " << p.to_string() << std::endl;
        }
    }

    return true;
}

} // namespace ruckig
