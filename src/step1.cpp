#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>


namespace ruckig {

Step1::Step1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax): p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af) {
    valid_profiles.reserve(3);

    pd = pf - p0;

    v0_v0 = v0 * v0;
    vf_vf = vf * vf;
    vMax_vMax = vMax * vMax;

    a0_a0 = a0 * a0;
    af_af = af * af;
    aMax_aMax = aMax * aMax;

    jMax_jMax = jMax * jMax;
}

void Step1::add_profile(Profile profile, Profile::Limits limits, double jMax) {
    profile.limits = limits;
    profile.direction = (jMax > 0) ? Profile::Direction::UP : Profile::Direction::DOWN;
    valid_profiles.emplace_back(profile);
}

void Step1::time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax_aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*(a0_a0*a0_a0 + af_af*af_af) + 8*aMax*(af_af*af - a0_a0*a0) + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 6*af_af*(aMax_aMax - 2*jMax*vf) - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf + 2*vMax) - jMax*(v0_v0 + vf_vf - 2*vMax_vMax)))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = profile.t[2];
    profile.t[5] = (af_af/2 - aMax_aMax - jMax*(vf - vMax))/(aMax*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
        add_profile(profile, Limits::ACC0_ACC1_VEL, jMax);
    }
}

void Step1::time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    profile.t[0] = (-a0*jMax + Sqrt(a0_a0/2 + jMax*(vMax - v0))*Abs(jMax))/(jMax_jMax);
    profile.t[1] = 0;
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (jMax*(3*Power(af,4) + 8*aMax*(Power(af,3) - Power(a0,3)) + 24*aMax*jMax*(a0*v0 - af*vf) + 6*af_af*(aMax_aMax - 2*jMax*vf) - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(vf + vMax) + jMax*(-vf_vf + Power(vMax,2)))) + 6*Sqrt(2)*aMax*Sqrt(a0_a0 + 2*jMax*(-v0 + vMax))*(a0_a0 - 2*jMax*(v0 + vMax))*Abs(jMax))/(24.*aMax*Power(jMax,3)*vMax);
    profile.t[4] = aMax/jMax;
    profile.t[5] = (af_af/2 - aMax_aMax + jMax*(vMax - vf))/(aMax*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
        add_profile(profile, Limits::ACC1_VEL, jMax);
    }
}

void Step1::time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax_aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*Power(a0,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*a0_a0*(aMax_aMax - 2*jMax*v0) - 24*af*aMax*jMax*vf + 6*Sqrt(2)*af_af*aMax*Sqrt(af_af + 2*jMax*(-vf + vMax))*jMax/Abs(jMax) - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vMax) + jMax*(-v0_v0 + Power(vMax,2)) + Sqrt(2)*aMax*(vf + vMax)*Sqrt(af_af + 2*jMax*(-vf + vMax))*jMax/Abs(jMax)))/(24.*aMax*jMax_jMax*vMax);
    profile.t[4] = Sqrt(af_af/2 + jMax*(-vf + vMax))/Abs(jMax);
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
        add_profile(profile, Limits::ACC0_VEL, jMax);
    }
}

void Step1::time_up_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Solution 3/4
    profile.t[0] = (-a0*jMax + Sqrt(a0_a0/2 + jMax*(-v0 + vMax))*Abs(jMax))/jMax_jMax;
    profile.t[1] = 0;
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (jMax*(-4*Power(a0,3) + 4*Power(af,3) + 12*jMax*(a0*v0 - af*vf) + 3*Sqrt(2)*af_af*Sqrt(af_af + 2*jMax*(-vf + vMax))*jMax/Abs(jMax) - 6*Abs(jMax)*(-2*Abs(jMax)*pd + Sqrt(2)*(vf + vMax)*Sqrt(af_af - 2*jMax*vf + 2*jMax*vMax))) + 3*Sqrt(2)*Sqrt(a0_a0 + 2*jMax*(-v0 + vMax))*(a0_a0 - 2*jMax*(v0 + vMax))*Abs(jMax))/(12*Power(jMax,3)*vMax);
    profile.t[4] = Sqrt(af_af/2 + jMax*(-vf + vMax))/Abs(jMax);
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
        add_profile(profile, Limits::VEL, jMax);
    }
}

void Step1::time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(6)*Sqrt(3*Power(a0,4) + 3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 6*af_af*(aMax_aMax - 2*jMax*vf) + 6*(Power(aMax,4) + 4*aMax*jMax_jMax*pd - 2*aMax_aMax*jMax*(v0 + vf) + 2*jMax_jMax*(v0_v0 + vf_vf)));

    // UDDU: Solution 2
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (6*a0_a0*aMax*jMax - 18*Power(aMax,3)*jMax - 12*aMax*jMax_jMax*v0 - h1*Abs(aMax)*Abs(jMax))/(12*aMax_aMax*jMax_jMax);
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = profile.t[2];
        profile.t[5] = profile.t[1] + (af_af/2 - a0_a0/2 + jMax*(v0 - vf))/(aMax*jMax);
        profile.t[6] = profile.t[4] + af/jMax;

        if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
            add_profile(profile, Limits::ACC0_ACC1, jMax);
        }
    }

    // UDDU: Solution 1
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (6*a0_a0*aMax*jMax - 18*Power(aMax,3)*jMax - 12*aMax*jMax_jMax*v0 + h1*Abs(aMax)*Abs(jMax))/(12*aMax_aMax*jMax_jMax);
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = profile.t[2];
        profile.t[5] = profile.t[1] + (af_af/2 - a0_a0/2 + jMax*(v0 - vf))/(aMax*jMax);
        profile.t[6] = profile.t[4] + af/jMax;

        if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
            add_profile(profile, Limits::ACC0_ACC1, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        const double h3 = 3*Power(a0,4) - 3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*aMax*jMax*(a0*v0 + af*vf);

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (h3 + 6*a0_a0*(3*aMax_aMax - 2*jMax*v0) - 6*af_af*(aMax_aMax + 2*jMax*vf) - 12*(2*Power(aMax,4) - 2*aMax*jMax_jMax*pd + aMax_aMax*jMax*(3*v0 + vf) + jMax_jMax*(-v0_v0 + vf_vf)))/(24*Power(aMax,3)*jMax);
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = profile.t[2];
        profile.t[5] = -(h3 + 6*a0_a0*(aMax_aMax - 2*jMax*v0) - 6*af_af*(3*aMax_aMax + 2*jMax*vf) + 12*(2*Power(aMax,4) + 2*aMax*jMax_jMax*pd - aMax_aMax*jMax*(v0 + 3*vf) + jMax_jMax*(v0_v0 - vf_vf)))/(24*Power(aMax,3)*jMax);
        profile.t[6] = profile.t[4] - af/jMax;

        if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, jMax, 0, -jMax})) {
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
    polynom[4] = (3*Power(a0,4) - 3*Power(af,4) + 8*Power(a0,3)*aMax - 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*a0_a0*(aMax_aMax + 2*jMax*v0) + 24*af*aMax*jMax*vf - 6*af_af*(aMax_aMax - 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf) + jMax*(v0_v0 - vf_vf)))/(12*Power(jMax,4));
    
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
        profile.t[5] = (a0_a0/2 + af_af/2 - aMax_aMax + 2*a0*jMax*t + jMax_jMax*t*t + jMax*(v0 - vf))/(aMax*jMax);
        profile.t[6] = profile.t[4] + af/jMax;

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
            add_profile(profile, Limits::ACC1, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (4*a0 - 2*aMax)/jMax;
        polynom[2] = (5*a0_a0 - 6*a0*aMax - aMax_aMax + 2*jMax*v0)/jMax_jMax;
        polynom[3] = (2*(Power(a0,3) - 2*a0_a0*aMax - a0*aMax_aMax + 2*a0*jMax*v0 - 2*aMax*jMax*v0))/Power(jMax,3);
        polynom[4] = -(-3*Power(a0,4) + 3*Power(af,4) + 8*Power(a0,3)*aMax - 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*a0_a0*(aMax_aMax - 2*jMax*v0) - 24*af*aMax*jMax*vf + 6*af_af*(aMax_aMax + 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf) + jMax*(-v0_v0 + vf_vf)))/(12.*Power(jMax,4));

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
            profile.t[5] = -(a0_a0 - af_af + 4*a0*jMax*t + 2*(aMax_aMax + jMax*(jMax*t*t + v0 - vf)))/(2*aMax*jMax);
            profile.t[6] = profile.t[4] - af/jMax;

            if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, jMax, 0, -jMax})) {
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
    polynom[4] = (-3*Power(a0,4) + 3*Power(af,4) + 8*Power(a0,3)*aMax - 8*Power(af,3)*aMax - 24*a0*aMax*jMax*v0 - 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 24*af*aMax*jMax*vf - 6*af_af*(aMax_aMax + 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf) + jMax*(-v0_v0 + vf_vf)))/(12*Power(jMax,4));

    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (a0_a0/2 - af_af/2 - aMax_aMax + jMax_jMax*t*t - jMax*(v0 - vf))/(aMax*jMax);
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = t;
        profile.t[5] = 0;
        profile.t[6] = profile.t[4] + af/jMax;

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
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
        polynom[4] = (-3*Power(a0,4) + 3*Power(af,4) + 8*Power(a0,3)*aMax - 8*Power(af,3)*aMax - 24*a0*aMax*jMax*v0 - 6*a0_a0*(aMax_aMax - 2*jMax*v0) - 24*af*aMax*jMax*vf + 6*af_af*(aMax_aMax + 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vf) + jMax*(-v0_v0 + vf_vf)))/(12*Power(jMax,4));

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0) {
                continue;
            }

            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (a0_a0 + af_af - 2*(aMax_aMax + jMax*(jMax*t*t + v0 - vf)))/(2*aMax*jMax);
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = 0;
            profile.t[4] = t;
            profile.t[5] = 0;
            profile.t[6] = profile.t[4] - af/jMax;
            
            if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, jMax, 0, -jMax})) {
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

        if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
            add_profile(profile, Limits::NONE, jMax);
        }
        return;
    }

    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (-4*(2*Power(a0,3) + Power(af,3) + 3*jMax_jMax*pd - 3*a0*(af_af + jMax*(v0 - 2*vf)) - 3*af*jMax*vf))/(3.*jMax*(-a0_a0 + af_af + 2*jMax*(v0 - vf)));
    polynom[2] = (-3*Power(a0,4) + Power(af_af + 2*jMax*(v0 - vf),2) - 8*a0*(Power(af,3) + 3*jMax_jMax*pd - 3*af*jMax*vf) + 2*a0_a0*(5*af_af - 2*jMax*(v0 + 5*vf)))/(2*jMax_jMax*(-a0_a0 + af_af + 2*jMax*(v0 - vf)));
    polynom[3] = -(Power(a0,5) + 8*a0_a0*(Power(af,3) + 3*jMax_jMax*pd - 3*af*jMax*vf) + 8*jMax*v0*(Power(af,3) + 3*jMax_jMax*pd - 3*af*jMax*vf) + Power(a0,3)*(-6*af_af + 4*jMax*(v0 + 3*vf)) - 3*a0*(Power(af,4) + 4*af_af*jMax*(v0 - vf) - 4*jMax_jMax*(v0_v0 + 2*v0*vf - vf_vf)))/(3.*Power(jMax,3)*(-a0_a0 + af_af + 2*jMax*(v0 - vf)));
    polynom[4] = -(Power(a0,6) + Power(af,6) - 48*Power(af,3)*jMax_jMax*pd + 144*af*Power(jMax,3)*pd*vf - 6*Power(af,4)*jMax*(3*v0 + vf) + 16*Power(a0,3)*(Power(af,3) + 3*jMax_jMax*pd - 3*af*jMax*vf) + 48*a0*jMax*v0*(Power(af,3) + 3*jMax_jMax*pd - 3*af*jMax*vf) - 36*af_af*jMax_jMax*(v0_v0 - 2*v0*vf - vf_vf) - 72*Power(jMax,3)*(jMax*Power(pd,2) - (v0 - vf)*Power(v0 + vf,2)) + Power(a0,4)*(-9*af_af + 6*jMax*(v0 + 3*vf)) - 9*a0_a0*(Power(af,4) + 4*af_af*jMax*(v0 - vf) - 4*jMax_jMax*(v0_v0 + 2*v0*vf - vf_vf)))/(36.*Power(jMax,4)*(-a0_a0 + af_af + 2*jMax*(v0 - vf)));

    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        // Refine root
        if (std::abs(Roots::polyEval(polynom, t)) > 1e-9) {
            t = Roots::shrinkInterval(polynom, t - 1e-5, t + 1e-5, 1e-14);
        }

        profile.t[0] = t;
        profile.t[1] = 0;
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = 0;
        profile.t[4] = Sqrt(a0_a0/2 + af_af/2 + 2*a0*jMax*profile.t[0] + jMax*(jMax*t*t + v0 - vf))/Abs(jMax);
        profile.t[5] = 0;
        profile.t[6] = profile.t[4] + af/jMax;
        
        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, -jMax, 0, jMax})) {
            add_profile(profile, Limits::NONE, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        std::array<double, 7> polynom;
        polynom[0] = 1.0;
        polynom[1] = (6*a0)/jMax;
        polynom[2] = (53*a0_a0 + af_af + 2*jMax*(7*v0 + vf))/(4.*jMax_jMax);
        polynom[3] = (40*Power(a0,3) - Power(af,3) - 3*jMax_jMax*pd - 3*af*jMax*vf + 3*a0*(af_af + 13*jMax*v0 + 2*jMax*vf))/(3.*Power(jMax,3));
        polynom[4] = (51*Power(a0,4) - Power(af,4) + 4*af_af*jMax*(v0 - vf) + 2*a0_a0*(5*af_af + 58*jMax*v0 + 10*jMax*vf) - 8*a0*(Power(af,3) + 3*jMax_jMax*pd + 3*af*jMax*vf) + 4*jMax_jMax*(7*v0_v0 + 2*v0*vf - vf_vf))/(8.*Power(jMax,4));
        polynom[5] = (17*Power(a0,5) + 2*Power(a0,3)*(3*af_af + 34*jMax*v0 + 6*jMax*vf) - 8*a0_a0*(Power(af,3) + 3*jMax_jMax*pd + 3*af*jMax*vf) - 8*jMax*v0*(Power(af,3) + 3*jMax_jMax*pd + 3*af*jMax*vf) - 3*a0*(Power(af,4) + 4*af_af*jMax*(-v0 + vf) + 4*jMax_jMax*(-5*v0_v0 - 2*v0*vf + vf_vf)))/(12*Power(jMax,5));
        polynom[6] = -(-17*Power(a0,6) + Power(af,6) - 48*Power(af,3)*jMax_jMax*pd - 144*af*Power(jMax,3)*pd*vf + 6*Power(af,4)*jMax*(3*v0 + vf) - 3*Power(a0,4)*(3*af_af + 34*jMax*v0 + 6*jMax*vf) + 16*Power(a0,3)*(Power(af,3) + 3*jMax_jMax*pd + 3*af*jMax*vf) + 48*a0*jMax*v0*(Power(af,3) + 3*jMax_jMax*pd + 3*af*jMax*vf) - 36*af_af*jMax_jMax*(v0_v0 - 2*v0*vf - vf_vf) - 72*Power(jMax,3)*(jMax*Power(pd,2) + (v0 - vf)*Power(v0 + vf,2)) + 9*a0_a0*(Power(af,4) + 4*af_af*jMax*(-v0 + vf) + 4*jMax_jMax*(-5*v0_v0 - 2*v0*vf + vf_vf)))/(144.*Power(jMax,6));

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
        double tz_max = {1000.0};
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
            profile.t[4] = Sqrt(af_af/2 - a0_a0/2 - 2*a0*jMax*t - jMax*(jMax*t*t + v0 - vf))/(Abs(jMax));
            profile.t[5] = 0;
            profile.t[6] = profile.t[4] - af/jMax;

            if (profile.check(pf, vf, af, vMax, aMax, {jMax, 0, -jMax, 0, jMax, 0, -jMax})) {
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

    if (valid_profiles.empty()) {
        return false;
    }

    // auto min_profile = std::min_element(valid_profiles.begin(), valid_profiles.end(), [](const Profile& a, const Profile& b){ return a.t_sum[6] < b.t_sum[6]; }); 

    double t_brake = valid_profiles[0].t_brake.value_or(0.0);
    block = Block {valid_profiles[0].t_sum[6] + t_brake, valid_profiles[0]};

    if (valid_profiles.size() == 2) {
        
    } else if (valid_profiles.size() == 3) {
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
            }
        }

    } else if (valid_profiles.size() > 3) {
        
    }

    // if (valid_profiles.size() > 1) {
    //     std::cout << "---\n " << valid_profiles.size() << std::endl;
    //     for (auto p: valid_profiles) {
    //         std::cout << p.t_sum[6] + t_brake << " " << p.to_string() << std::endl;
    //     }
    // }

    return true;
}

} // namespace ruckig
