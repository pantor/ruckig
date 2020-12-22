#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>
#include <ruckig/wolfram.hpp>


namespace ruckig {

RuckigStep1::RuckigStep1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax): p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af) {

}

void RuckigStep1::add_profile(Profile profile, Profile::Limits limits, double jMax) {
    profile.limits = limits;
    profile.direction = (jMax > 0) ? Profile::Direction::UP : Profile::Direction::DOWN;

    if (valid_profiles.empty()) {
        fastest = profile;
    }

    valid_profiles.push_back(profile);
}

void RuckigStep1::time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    double a0_a0 = a0 * a0;
    double aMax_aMax = aMax * aMax;

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax_aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = aMax/jMax;
    profile.t[3] = (3*a0_a0*a0_a0 + 3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*a0_a0*(aMax_aMax - 2*jMax*v0) - 24*af*aMax*jMax*vf + 6*Power(af,2)*(aMax_aMax - 2*jMax*vf) - 12*jMax*(2*aMax*jMax*(p0 - pf) + aMax_aMax*(v0 + vf + 2*vMax) - jMax*(Power(v0,2) + Power(vf,2) - 2*Power(vMax,2))))/(24.*aMax*Power(jMax,2)*vMax);
    profile.t[4] = profile.t[2];
    profile.t[5] = (Power(af,2)/2 - aMax_aMax - jMax*vf + jMax*vMax)/(aMax*jMax);
    profile.t[6] = (af + aMax)/jMax;

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    if (profile.check(pf, vf, af, vMax, aMax)) {
        add_profile(profile, Limits::ACC0_ACC1_VEL, jMax);
    }
}

void RuckigStep1::time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    profile.t[0] = (-2*a0*jMax + Sqrt(2)*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*Abs(jMax))/(2*Power(jMax,2));
    profile.t[1] = 0;
    profile.t[2] = Sqrt(Power(a0,2)/2 + jMax*(-v0 + vMax))/Abs(jMax);
    profile.t[3] = (jMax*(3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 - 24*af*aMax*jMax*vf + 6*Power(af,2)*(Power(aMax,2) - 2*jMax*vf) - 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(vf + vMax) + jMax*(-Power(vf,2) + Power(vMax,2)))) + 6*Sqrt(2)*aMax*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*(Power(a0,2) - 2*jMax*(v0 + vMax))*Abs(jMax))/(24.*aMax*Power(jMax,3)*vMax);
    profile.t[4] = aMax/jMax;
    profile.t[5] = (Power(af,2)/2 - Power(aMax,2) - jMax*vf + jMax*vMax)/(aMax*jMax);
    profile.t[6] = (af + aMax)/jMax;

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    if (profile.check(pf, vf, af, vMax, aMax)) {
        add_profile(profile, Limits::ACC1_VEL, jMax);
    }
}

void RuckigStep1::time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Solution 2
    if (jMax > 0) {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (Power(a0,2)/2 - Power(aMax,2) - jMax*(v0 - vMax))/(aMax*jMax);
        profile.t[2] = aMax/jMax;
        profile.t[3] = (3*Power(a0,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) - 24*af*aMax*jMax*vf + 6*Sqrt(2)*Power(af,2)*aMax*Sqrt(Power(af,2) + 2*jMax*(-vf + vMax)) - 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(v0 + vMax) + jMax*(-Power(v0,2) + Power(vMax,2)) + Sqrt(2)*aMax*(vf + vMax)*Sqrt(Power(af,2) + 2*jMax*(-vf + vMax))))/(24.*aMax*Power(jMax,2)*vMax);
        profile.t[4] = Sqrt(Power(af,2)/2. + jMax*(-vf + vMax))/jMax;
        profile.t[5] = 0;
        profile.t[6] = (af + Sqrt(Power(af,2)/2. + jMax*(-vf + vMax)))/jMax;

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, af, vMax, aMax)) {
            add_profile(profile, Limits::ACC0_VEL, jMax);
        }
        return;
    }

    // Solution 1
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (Power(a0,2)/2 - Power(aMax,2) - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = aMax/jMax;
    profile.t[3] = (3*Power(a0,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) - 24*af*aMax*jMax*vf - 6*Sqrt(2)*Power(af,2)*aMax*Sqrt(Power(af,2) + 2*jMax*(-vf + vMax)) - 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(v0 + vMax) + jMax*(-Power(v0,2) + Power(vMax,2)) - Sqrt(2)*aMax*(vf + vMax)*Sqrt(Power(af,2) + 2*jMax*(-vf + vMax))))/(24.*aMax*Power(jMax,2)*vMax);
    profile.t[4] = -(Sqrt(Power(af,2)/2. + jMax*(-vf + vMax))/jMax);
    profile.t[5] = 0;
    profile.t[6] = (2*af - Sqrt(2)*Sqrt(Power(af,2) - 2*jMax*vf + 2*jMax*vMax))/(2.*jMax);

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    if (profile.check(pf, vf, af, vMax, aMax)) {
        add_profile(profile, Limits::ACC0_VEL, jMax);
    }
}

void RuckigStep1::time_up_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Solution 4
    if (jMax > 0) {
        profile.t[0] = ((-2*a0*jMax + Sqrt(2)*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*Abs(jMax))/(2*Power(jMax,2)));
        profile.t[1] = 0;
        profile.t[2] = Sqrt(Power(a0,2)/2 + jMax*(-v0 + vMax))/Abs(jMax);
        profile.t[3] = (jMax*(-4*Power(a0,3) + 4*Power(af,3) + 12*a0*jMax*v0 - 12*af*jMax*vf + 3*Sqrt(2)*Power(af,2)*Sqrt(Power(af,2) + 2*jMax*(-vf + vMax)) - 6*jMax*(2*jMax*(p0 - pf) + Sqrt(2)*(vf + vMax)*Sqrt(Power(af,2) - 2*jMax*vf + 2*jMax*vMax))) + 3*Sqrt(2)*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*(Power(a0,2) - 2*jMax*(v0 + vMax))*Abs(jMax))/(12.*Power(jMax,3)*vMax);
        profile.t[4] = Sqrt(Power(af,2)/2. + jMax*(-vf + vMax))/jMax;
        profile.t[5] = 0;
        profile.t[6] = (af + Sqrt(Power(af,2)/2. + jMax*(-vf + vMax)))/jMax;

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, af, vMax, aMax)) {
            add_profile(profile, Limits::VEL, jMax);
        }
        return;
    }
    
    // Solution 3
    profile.t[0] = ((-2*a0*jMax + Sqrt(2)*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*Abs(jMax))/(2*Power(jMax,2)));
    profile.t[1] = 0;
    profile.t[2] = Sqrt(Power(a0,2)/2 + jMax*(-v0 + vMax))/Abs(jMax);
    profile.t[3] = (jMax*(-4*Power(a0,3) + 4*Power(af,3) + 12*a0*jMax*v0 - 12*af*jMax*vf - 3*Sqrt(2)*Power(af,2)*Sqrt(Power(af,2) + 2*jMax*(-vf + vMax)) + 6*jMax*(-2*jMax*(p0 - pf) + Sqrt(2)*(vf + vMax)*Sqrt(Power(af,2) - 2*jMax*vf + 2*jMax*vMax))) + 3*Sqrt(2)*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*(Power(a0,2) - 2*jMax*(v0 + vMax))*Abs(jMax))/(12.*Power(jMax,3)*vMax);
    profile.t[4] = -(Sqrt(Power(af,2)/2. + jMax*(-vf + vMax))/jMax);
    profile.t[5] = 0;
    profile.t[6] = (2*af - Sqrt(2)*Sqrt(Power(af,2) - 2*jMax*vf + 2*jMax*vMax))/(2.*jMax);

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    if (profile.check(pf, vf, af, vMax, aMax)) {
        add_profile(profile, Limits::VEL, jMax);
    }
}

void RuckigStep1::time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(6)*Sqrt(Power(aMax,2)*(3*Power(a0,4) + 3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) - 24*af*aMax*jMax*vf + 6*Power(af,2)*(Power(aMax,2) - 2*jMax*vf) + 6*(Power(aMax,4) + 4*aMax*Power(jMax,2)*(-p0 + pf) - 2*Power(aMax,2)*jMax*(v0 + vf) + 2*Power(jMax,2)*(Power(v0,2) + Power(vf,2)))))*Abs(jMax);

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (6*Power(a0,2)*aMax*jMax - 18*Power(aMax,3)*jMax - 12*aMax*Power(jMax,2)*v0 + h1)/(12.*Power(aMax,2)*Power(jMax,2));
    profile.t[2] = aMax/jMax;
    profile.t[3] = 0;
    profile.t[4] = profile.t[2];
    profile.t[5] = (6*Power(af,2)*aMax*jMax - 18*Power(aMax,3)*jMax - 12*aMax*Power(jMax,2)*vf + h1)/(12.*Power(aMax,2)*Power(jMax,2));
    profile.t[6] = (af + aMax)/jMax;

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    if (profile.check(pf, vf, af, vMax, aMax)) {
        add_profile(profile, Limits::ACC0_ACC1, jMax);
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (3*Power(a0,4) - 3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(3*Power(aMax,2) - 2*jMax*v0) + 24*af*aMax*jMax*vf - 6*Power(af,2)*(Power(aMax,2) + 2*jMax*vf) - 12*(2*Power(aMax,4) + 2*aMax*Power(jMax,2)*(p0 - pf) + Power(aMax,2)*jMax*(3*v0 + vf) + Power(jMax,2)*(-Power(v0,2) + Power(vf,2))))/(24.*Power(aMax,3)*jMax);
        profile.t[2] = aMax/jMax;
        profile.t[3] = 0;
        profile.t[4] = profile.t[2];
        profile.t[5] = -(3*Power(a0,4) - 3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) + 24*af*aMax*jMax*vf - 6*Power(af,2)*(3*Power(aMax,2) + 2*jMax*vf) + 12*(2*Power(aMax,4) + 2*aMax*Power(jMax,2)*(-p0 + pf) - Power(aMax,2)*jMax*(v0 + 3*vf) + Power(jMax,2)*(Power(v0,2) - Power(vf,2))))/(24.*Power(aMax,3)*jMax);
        profile.t[6] = (-af + aMax)/jMax;

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
        if (profile.check(pf, vf, af, vMax, aMax)) {
            add_profile(profile, Limits::ACC0_ACC1, jMax);
        }
    }
}

void RuckigStep1::time_up_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (2*(2*a0 + aMax))/jMax;
    polynom[2] = (5*Power(a0,2) + 6*a0*aMax + Power(aMax,2) + 2*jMax*v0)/Power(jMax,2);
    polynom[3] = (2*(a0 + aMax)*(Power(a0,2) + a0*aMax + 2*jMax*v0))/Power(jMax,3);
    polynom[4] = (3*Power(a0,4) - 3*Power(af,4) + 8*Power(a0,3)*aMax - 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) + 2*jMax*v0) + 24*af*aMax*jMax*vf - 6*Power(af,2)*(Power(aMax,2) - 2*jMax*vf) + 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(v0 + vf) + jMax*(Power(v0,2) - Power(vf,2))))/(12.*Power(jMax,4));
    
    auto roots = Roots::solveQuart(polynom);
    for (double t: roots) {
        if (t <= 0.0) {
            continue;
        }

        profile.t[0] = t;
        profile.t[1] = 0;
        profile.t[2] = a0/jMax + t;
        profile.t[3] = 0;
        profile.t[4] = aMax/jMax;
        profile.t[5] = (Power(a0,2) + Power(af,2) - 2*Power(aMax,2) + 4*a0*jMax*t + 2*Power(jMax,2)*Power(t,2) + 2*jMax*v0 - 2*jMax*vf)/(2.*aMax*jMax);
        profile.t[6] = (af + aMax)/jMax;

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, af, vMax, aMax)) {
            add_profile(profile, Limits::ACC1, jMax);
        }
    }
}

void RuckigStep1::time_up_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (2*aMax)/jMax;
    polynom[2] = (-Power(af,2) + Power(aMax,2) + 2*jMax*vf)/Power(jMax,2);
    polynom[3] = (-2*aMax*(Power(af,2) - 2*jMax*vf))/Power(jMax,3);
    polynom[4] = (-3*Power(a0,4) + 3*Power(af,4) + 8*Power(a0,3)*aMax - 8*Power(af,3)*aMax - 24*a0*aMax*jMax*v0 - 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) + 24*af*aMax*jMax*vf - 6*Power(af,2)*(Power(aMax,2) + 2*jMax*vf) + 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(v0 + vf) + jMax*(-Power(v0,2) + Power(vf,2))))/(12.*Power(jMax,4));

    auto roots = Roots::solveQuart(polynom);
    for (double t: roots) {
        if (t <= 0.0) {
            continue;
        }

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (Power(a0,2) - Power(af,2) - 2*Power(aMax,2) + 2*Power(jMax,2)*Power(t,2) - 2*jMax*v0 + 2*jMax*vf)/(2.*aMax*jMax);
        profile.t[2] = aMax/jMax;
        profile.t[3] = 0;
        profile.t[4] = t;
        profile.t[5] = 0;
        profile.t[6] = af/jMax + t;

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, af, vMax, aMax)) {
            add_profile(profile, Limits::ACC0, jMax);
        }
    }
}

void RuckigStep1::time_up_none(Profile& profile, double vMax, double aMax, double jMax) {
    if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        profile.t[0] = std::cbrt((pf - p0)/(2*jMax));
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, af, vMax, aMax)) {
            add_profile(profile, Limits::NONE, jMax);
        }
        return;
    }

    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (-4*(2*Power(a0,3) + Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) - 3*a0*(Power(af,2) + jMax*(v0 - 2*vf)) - 3*af*jMax*vf))/(3.*jMax*(-Power(a0,2) + Power(af,2) + 2*jMax*(v0 - vf)));
    polynom[2] = (-3*Power(a0,4) + Power(Power(af,2) + 2*jMax*(v0 - vf),2) - 8*a0*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) - 3*af*jMax*vf) + 2*Power(a0,2)*(5*Power(af,2) - 2*jMax*(v0 + 5*vf)))/(2.*Power(jMax,2)*(-Power(a0,2) + Power(af,2) + 2*jMax*(v0 - vf)));
    polynom[3] = -(Power(a0,5) + 8*Power(a0,2)*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) - 3*af*jMax*vf) + 8*jMax*v0*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) - 3*af*jMax*vf) + Power(a0,3)*(-6*Power(af,2) + 4*jMax*(v0 + 3*vf)) - 3*a0*(Power(af,4) + 4*Power(af,2)*jMax*(v0 - vf) - 4*Power(jMax,2)*(Power(v0,2) + 2*v0*vf - Power(vf,2))))/(3.*Power(jMax,3)*(-Power(a0,2) + Power(af,2) + 2*jMax*(v0 - vf)));
    polynom[4] = -(Power(a0,6) + Power(af,6) + 48*Power(af,3)*Power(jMax,2)*(p0 - pf) - 144*af*Power(jMax,3)*(p0 - pf)*vf - 6*Power(af,4)*jMax*(3*v0 + vf) + 16*Power(a0,3)*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) - 3*af*jMax*vf) + 48*a0*jMax*v0*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) - 3*af*jMax*vf) - 36*Power(af,2)*Power(jMax,2)*(Power(v0,2) - 2*v0*vf - Power(vf,2)) - 72*Power(jMax,3)*(jMax*Power(p0 - pf,2) - (v0 - vf)*Power(v0 + vf,2)) + Power(a0,4)*(-9*Power(af,2) + 6*jMax*(v0 + 3*vf)) - 9*Power(a0,2)*(Power(af,4) + 4*Power(af,2)*jMax*(v0 - vf) - 4*Power(jMax,2)*(Power(v0,2) + 2*v0*vf - Power(vf,2))))/(36.*Power(jMax,4)*(-Power(a0,2) + Power(af,2) + 2*jMax*(v0 - vf)));

    auto roots = Roots::solveQuart(polynom);
    for (double t: roots) {
        if (t <= 0.0) {
            continue;
        }

        // Refine root
        if (std::abs(Roots::polyEval(polynom, t)) > 1e-9) {
            t = Roots::shrinkInterval(polynom, t - 1e-5, t + 1e-5, 1e-14);
        }

        profile.t[0] = t;
        profile.t[1] = 0;
        profile.t[2] = a0/jMax + t;
        profile.t[3] = 0;
        profile.t[4] = (2*a0*jMax + 2*Power(jMax,2)*profile.t[0] - 2*Power(jMax,2)*profile.t[2] + Sqrt(2)*Sqrt(Power(jMax,2)*(Power(a0,2) + Power(af,2) + 4*a0*jMax*profile.t[0] + 2*jMax*(jMax*Power(profile.t[0],2) + v0 - vf))))/(2.*Power(jMax,2));
        profile.t[5] = 0;
        profile.t[6] = (2*af*jMax + Sqrt(2)*Sqrt(Power(jMax,2)*(Power(a0,2) + Power(af,2) + 4*a0*jMax*profile.t[0] + 2*jMax*(jMax*Power(profile.t[0],2) + v0 - vf))))/(2.*Power(jMax,2));
        
        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];
            
        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, af, vMax, aMax)) {
            add_profile(profile, Limits::NONE, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        std::array<double, 7> polynom;
        polynom[0] = 1.0;
        polynom[1] = (6*a0)/jMax;
        polynom[2] = (53*Power(a0,2) + Power(af,2) + 2*jMax*(7*v0 + vf))/(4.*Power(jMax,2));
        polynom[3] = (40*Power(a0,3) - Power(af,3) + 3*Power(jMax,2)*(p0 - pf) - 3*af*jMax*vf + 3*a0*(Power(af,2) + 13*jMax*v0 + 2*jMax*vf))/(3.*Power(jMax,3));
        polynom[4] = (51*Power(a0,4) - Power(af,4) + 4*Power(af,2)*jMax*(v0 - vf) + 2*Power(a0,2)*(5*Power(af,2) + 58*jMax*v0 + 10*jMax*vf) - 8*a0*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) + 4*Power(jMax,2)*(7*Power(v0,2) + 2*v0*vf - Power(vf,2)))/(8.*Power(jMax,4));
        polynom[5] = (17*Power(a0,5) + 2*Power(a0,3)*(3*Power(af,2) + 34*jMax*v0 + 6*jMax*vf) - 8*Power(a0,2)*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) - 8*jMax*v0*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) - 3*a0*(Power(af,4) + 4*Power(af,2)*jMax*(-v0 + vf) + 4*Power(jMax,2)*(-5*Power(v0,2) - 2*v0*vf + Power(vf,2))))/(12.*Power(jMax,5));
        polynom[6] = -(-17*Power(a0,6) + Power(af,6) + 48*Power(af,3)*Power(jMax,2)*(p0 - pf) + 144*af*Power(jMax,3)*(p0 - pf)*vf + 6*Power(af,4)*jMax*(3*v0 + vf) - 3*Power(a0,4)*(3*Power(af,2) + 34*jMax*v0 + 6*jMax*vf) + 16*Power(a0,3)*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) + 48*a0*jMax*v0*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) - 36*Power(af,2)*Power(jMax,2)*(Power(v0,2) - 2*v0*vf - Power(vf,2)) - 72*Power(jMax,3)*(jMax*Power(p0 - pf,2) + (v0 - vf)*Power(v0 + vf,2)) + 9*Power(a0,2)*(Power(af,4) + 4*Power(af,2)*jMax*(-v0 + vf) + 4*Power(jMax,2)*(-5*Power(v0,2) - 2*v0*vf + Power(vf,2))))/(144.*Power(jMax,6));

        std::array<double, 6> deriv;
        deriv[0] = 6 * polynom[0];
        deriv[1] = 5 * polynom[1];
        deriv[2] = 4 * polynom[2];
        deriv[3] = 3 * polynom[3];
        deriv[4] = 2 * polynom[4];
        deriv[5] = polynom[5];

        auto dd_extremas = Roots::solveQuart(5 * deriv[0], 4 * deriv[1], 3 * deriv[2], 2 * deriv[3], deriv[4]);
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

            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + t;
            profile.t[3] = 0;
            profile.t[4] = -((29*Power(a0,7) + 12*Power(af,7) + 245*Power(a0,6)*jMax*t - Power(af,6)*jMax*t - 12*Power(af,5)*jMax*(jMax*Power(t,2) + 2*v0 - 7*vf) - 6*Power(af,4)*Power(jMax,2)*(6*p0 - 6*pf + t*(3*jMax*Power(t,2) + 3*v0 + vf)) + 3*Power(a0,5)*(3*Power(af,2) + 2*jMax*(117*jMax*Power(t,2) + 7*v0 + 3*vf)) + 36*Power(af,2)*Power(jMax,3)*(2*Power(jMax,2)*Power(t,5) - 2*pf*v0 + 5*t*Power(v0,2) + 2*p0*(v0 - 2*vf) + 4*pf*vf - 2*t*v0*vf - t*Power(vf,2) + jMax*Power(t,2)*(p0 - pf + 6*t*v0 - 2*t*vf)) - 72*af*Power(jMax,3)*vf*(2*Power(jMax,2)*Power(t,4) + 2*(v0 - vf)*vf + jMax*t*(2*p0 - 2*pf + 3*t*v0 + t*vf)) - 12*Power(af,3)*Power(jMax,2)*(4*Power(jMax,2)*Power(t,4) + 2*(5*v0 - 8*vf)*vf + jMax*t*(4*p0 - 4*pf + 6*t*v0 + 5*t*vf)) + Power(a0,4)*(-28*Power(af,3) + 105*Power(af,2)*jMax*t - 84*af*jMax*vf + 6*Power(jMax,2)*(14*p0 - 14*pf + 137*jMax*Power(t,3) - 11*t*v0 + 35*t*vf)) - 3*Power(a0,2)*jMax*(15*Power(af,4)*t + Power(af,3)*(84*jMax*Power(t,2) + 8*v0) - 12*Power(af,2)*jMax*t*(17*jMax*Power(t,2) + 11*v0 - 5*vf) + 12*af*jMax*(21*jMax*Power(t,2) + 2*v0)*vf - 12*Power(jMax,2)*(2*Power(jMax,2)*Power(t,5) + 2*p0*v0 - 2*pf*v0 - 31*t*Power(v0,2) + 22*t*v0*vf - 5*t*Power(vf,2) + jMax*Power(t,2)*(21*p0 - 21*pf - 46*t*v0 + 34*t*vf))) - Power(a0,3)*(21*Power(af,4) + 136*Power(af,3)*jMax*t - 12*Power(af,2)*jMax*(34*jMax*Power(t,2) + 5*v0 - 7*vf) + 408*af*Power(jMax,2)*t*vf - 12*Power(jMax,2)*(34*Power(jMax,2)*Power(t,4) - 15*Power(v0,2) + 10*v0*vf - 7*Power(vf,2) + jMax*t*(34*p0 - 34*pf - 87*t*v0 + 68*t*vf))) + 72*Power(jMax,4)*(2*Power(jMax,2)*Power(t,4)*(p0 - pf + t*(-v0 + vf)) - (v0 - vf)*(2*(-p0 + pf)*vf + t*(3*Power(v0,2) - 2*v0*vf - Power(vf,2))) + jMax*t*(Power(p0,2) - 2*p0*pf + Power(pf,2) + p0*t*(3*v0 + vf) - pf*t*(3*v0 + vf) - Power(t,2)*(5*Power(v0,2) - 6*v0*vf + Power(vf,2)))) - a0*(Power(af,6) + 24*Power(af,5)*jMax*t + 6*Power(af,4)*jMax*(9*jMax*Power(t,2) - 3*v0 + vf) + 24*Power(af,3)*Power(jMax,2)*(2*p0 - 2*pf + 8*jMax*Power(t,3) + 4*t*v0 + 5*t*vf) + 144*af*Power(jMax,3)*vf*(p0 - pf + t*(4*jMax*Power(t,2) + 2*v0 + vf)) - 36*Power(af,2)*Power(jMax,2)*(10*Power(jMax,2)*Power(t,4) + 3*Power(v0,2) + 2*v0*vf - Power(vf,2) + jMax*t*(2*p0 - 2*pf + 17*t*v0 - 6*t*vf)) - 72*Power(jMax,3)*(-3*Power(v0,3) + 3*Power(v0,2)*vf + v0*Power(vf,2) - Power(vf,3) + 2*Power(jMax,2)*Power(t,3)*(4*p0 - 4*pf - 6*t*v0 + 5*t*vf) + jMax*(Power(p0,2) - 2*p0*pf + Power(pf,2) + 2*p0*t*(2*v0 + vf) - 2*pf*t*(2*v0 + vf) + Power(t,2)*(-18*Power(v0,2) + 17*v0*vf - 3*Power(vf,2))))))/(jMax*(Power(a0,6) - 17*Power(af,6) + 48*Power(af,3)*Power(jMax,2)*(p0 - pf) + 6*Power(af,4)*jMax*(3*v0 - 17*vf) + 144*af*Power(jMax,3)*(p0 - pf)*vf + 3*Power(a0,4)*(3*Power(af,2) - 2*jMax*v0 + 6*jMax*vf) + 16*Power(a0,3)*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) - 48*a0*jMax*v0*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) + 36*Power(af,2)*Power(jMax,2)*(Power(v0,2) + 2*v0*vf - 5*Power(vf,2)) - 72*Power(jMax,3)*(jMax*Power(p0 - pf,2) + Power(v0 - vf,2)*(v0 + vf)) - 9*Power(a0,2)*(Power(af,4) + 4*Power(af,2)*jMax*(v0 + vf) + 4*Power(jMax,2)*(-Power(v0,2) + 2*v0*vf + Power(vf,2))))));
            profile.t[5] = 0;
            profile.t[6] = -((29*Power(a0,7) - 5*Power(af,7) - Power(af,6)*jMax*t + Power(a0,6)*(af + 245*jMax*t) - 6*Power(af,5)*jMax*(2*jMax*Power(t,2) + v0 + 3*vf) - 6*Power(af,4)*Power(jMax,2)*(-2*p0 + 2*pf + t*(3*jMax*Power(t,2) + 3*v0 + vf)) + 3*Power(a0,5)*(3*Power(af,2) + 2*jMax*(117*jMax*Power(t,2) + 7*v0 + 3*vf)) + 36*Power(af,2)*Power(jMax,3)*(2*Power(jMax,2)*Power(t,5) + 2*p0*v0 - 2*pf*v0 + 5*t*Power(v0,2) - 2*t*v0*vf - t*Power(vf,2) + jMax*Power(t,2)*(p0 - pf + 6*t*v0 - 2*t*vf)) - 12*Power(af,3)*Power(jMax,2)*(4*Power(jMax,2)*Power(t,4) - 3*Power(v0,2) + 4*v0*vf - Power(vf,2) + jMax*t*(4*p0 - 4*pf + 6*t*v0 + 5*t*vf)) + Power(a0,4)*(-19*Power(af,3) + 105*Power(af,2)*jMax*t - 6*af*jMax*(v0 + 11*vf) + 6*Power(jMax,2)*(14*p0 - 14*pf + 137*jMax*Power(t,3) - 11*t*v0 + 35*t*vf)) - 72*af*Power(jMax,3)*(Power(v0,3) + 2*Power(jMax,2)*Power(t,4)*vf - Power(v0,2)*vf + v0*Power(vf,2) - Power(vf,3) + jMax*(Power(p0,2) + Power(pf,2) - 2*pf*t*vf + Power(t,2)*vf*(3*v0 + vf) - 2*p0*(pf - t*vf))) - 3*Power(a0,2)*(3*Power(af,5) + 15*Power(af,4)*jMax*t - 12*Power(af,2)*Power(jMax,2)*t*(17*jMax*Power(t,2) + 11*v0 - 5*vf) + 4*Power(af,3)*jMax*(21*jMax*Power(t,2) + 5*v0 + 3*vf) + 12*af*Power(jMax,2)*(-Power(v0,2) + 4*v0*vf + vf*(21*jMax*Power(t,2) + vf)) - 12*Power(jMax,3)*(2*Power(jMax,2)*Power(t,5) + 2*p0*v0 - 2*pf*v0 - 31*t*Power(v0,2) + 22*t*v0*vf - 5*t*Power(vf,2) + jMax*Power(t,2)*(21*p0 - 21*pf - 46*t*v0 + 34*t*vf))) - Power(a0,3)*(5*Power(af,4) + 136*Power(af,3)*jMax*t - 12*Power(af,2)*jMax*(34*jMax*Power(t,2) + 5*v0 - 3*vf) + 24*af*Power(jMax,2)*(2*p0 - 2*pf + 17*t*vf) - 12*Power(jMax,2)*(34*Power(jMax,2)*Power(t,4) - 15*Power(v0,2) + 10*v0*vf - 7*Power(vf,2) + jMax*t*(34*p0 - 34*pf - 87*t*v0 + 68*t*vf))) + 72*Power(jMax,4)*(2*Power(jMax,2)*Power(t,4)*(p0 - pf + t*(-v0 + vf)) - (v0 - vf)*(2*(-p0 + pf)*vf + t*(3*Power(v0,2) - 2*v0*vf - Power(vf,2))) + jMax*t*(Power(p0,2) - 2*p0*pf + Power(pf,2) + p0*t*(3*v0 + vf) - pf*t*(3*v0 + vf) - Power(t,2)*(5*Power(v0,2) - 6*v0*vf + Power(vf,2)))) - a0*(Power(af,6) + 24*Power(af,5)*jMax*t + 6*Power(af,4)*jMax*(9*jMax*Power(t,2) + 5*v0 + vf) + 24*Power(af,3)*Power(jMax,2)*(2*p0 - 2*pf + 8*jMax*Power(t,3) + 4*t*v0 + 5*t*vf) + 144*af*Power(jMax,3)*(pf*(v0 - vf) + p0*(-v0 + vf) + t*vf*(4*jMax*Power(t,2) + 2*v0 + vf)) - 36*Power(af,2)*Power(jMax,2)*(10*Power(jMax,2)*Power(t,4) + 3*Power(v0,2) - 2*v0*vf - Power(vf,2) + jMax*t*(2*p0 - 2*pf + 17*t*v0 - 6*t*vf)) - 72*Power(jMax,3)*(-3*Power(v0,3) + 3*Power(v0,2)*vf + v0*Power(vf,2) - Power(vf,3) + 2*Power(jMax,2)*Power(t,3)*(4*p0 - 4*pf - 6*t*v0 + 5*t*vf) + jMax*(Power(p0,2) - 2*p0*pf + Power(pf,2) + 2*p0*t*(2*v0 + vf) - 2*pf*t*(2*v0 + vf) + Power(t,2)*(-18*Power(v0,2) + 17*v0*vf - 3*Power(vf,2))))))/(jMax*(Power(a0,6) - 17*Power(af,6) + 48*Power(af,3)*Power(jMax,2)*(p0 - pf) + 6*Power(af,4)*jMax*(3*v0 - 17*vf) + 144*af*Power(jMax,3)*(p0 - pf)*vf + 3*Power(a0,4)*(3*Power(af,2) - 2*jMax*v0 + 6*jMax*vf) + 16*Power(a0,3)*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) - 48*a0*jMax*v0*(Power(af,3) + 3*Power(jMax,2)*(-p0 + pf) + 3*af*jMax*vf) + 36*Power(af,2)*Power(jMax,2)*(Power(v0,2) + 2*v0*vf - 5*Power(vf,2)) - 72*Power(jMax,3)*(jMax*Power(p0 - pf,2) + Power(v0 - vf,2)*(v0 + vf)) - 9*Power(a0,2)*(Power(af,4) + 4*Power(af,2)*jMax*(v0 + vf) + 4*Power(jMax,2)*(-Power(v0,2) + 2*v0*vf + Power(vf,2))))));

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(pf, vf, af, vMax, aMax)) {
                add_profile(profile, Limits::NONE, jMax);
            }
        }
    }
}

void RuckigStep1::time_down_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, -vMax, -aMax, -jMax);
}

void RuckigStep1::time_down_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc1_vel(profile, -vMax, -aMax, -jMax);
}

void RuckigStep1::time_down_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_vel(profile, -vMax, -aMax, -jMax);
}

void RuckigStep1::time_down_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_vel(profile, -vMax, -aMax, -jMax);
}

void RuckigStep1::time_down_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, -vMax, -aMax, -jMax);
}

void RuckigStep1::time_down_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc1(profile, -vMax, -aMax, -jMax);
}

void RuckigStep1::time_down_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0(profile, -vMax, -aMax, -jMax);
}

void RuckigStep1::time_down_none(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_none(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep1::get_profile(const Profile& input, double vMax, double aMax, double jMax) {
    Profile profile = input;

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

    block = Block {fastest.t_sum[6] + fastest.t_brake.value_or(0.0)};

    if (valid_profiles.size() == 2) {
        
    } else if (valid_profiles.size() == 3) {
        auto& p_a = valid_profiles[1];
        auto& p_b = valid_profiles[2];
        if (p_a.direction == p_b.direction) {
            block.a = Block::Interval {p_a.t_sum[6], p_b.t_sum[6]};
        }
        
    } else if (valid_profiles.size() > 3) {
        
    }

    // if (valid_profiles.size() > 1) {
    //     std::cout << "---\n EVEN MORE " << valid_profiles.size() << std::endl;
    //     for (auto p: valid_profiles) {
    //         std::cout << p.t_sum[6] << " " << p.to_string() << std::endl;
    //     }
    // }

    return !valid_profiles.empty();
}

} // namespace ruckig
