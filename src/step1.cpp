#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>
#include <ruckig/wolfram.hpp>


namespace ruckig {

RuckigStep1::RuckigStep1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {

}

bool RuckigStep1::time_up_acc0_acc1_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
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
    return profile.check(pf, vf, af, vMax, aMax);
}

bool RuckigStep1::time_up_acc1_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    profile.t[0] = (-2*a0*jMax + Sqrt(2)*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*Abs(jMax))/(2*Power(jMax,2));
    profile.t[1] = 0;
    profile.t[2] = Sqrt(Power(a0,2)/2 + jMax*(-v0 + vMax))/Abs(jMax);
    profile.t[3] = (jMax*(3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 - 24*af*aMax*jMax*vf + 6*Power(af,2)*(Power(aMax,2) - 2*jMax*vf) - 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(vf + vMax) + jMax*(-Power(vf,2) + Power(vMax,2)))) + 6*Sqrt(2)*aMax*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*(Power(a0,2) - 2*jMax*(v0 + vMax))*Abs(jMax))/(24.*aMax*Power(jMax,3)*vMax);
    profile.t[4] = aMax/jMax;
    profile.t[5] = (Power(af,2)/2 - Power(aMax,2) - jMax*vf + jMax*vMax)/(aMax*jMax);
    profile.t[6] = (af + aMax)/jMax;

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(pf, vf, af, vMax, aMax);
}

bool RuckigStep1::time_up_acc0_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
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
        return profile.check(pf, vf, af, vMax, aMax);
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
    return profile.check(pf, vf, af, vMax, aMax);
}

bool RuckigStep1::time_up_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
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
        return profile.check(pf, vf, af, vMax, aMax);
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
    return profile.check(pf, vf, af, vMax, aMax);
}

bool RuckigStep1::time_up_acc0_acc1(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(6)*Sqrt(Power(aMax,2)*(3*Power(a0,4) + 3*Power(af,4) - 8*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) - 24*af*aMax*jMax*vf + 6*Power(af,2)*(Power(aMax,2) - 2*jMax*vf) + 6*(Power(aMax,4) + 4*aMax*Power(jMax,2)*(-p0 + pf) - 2*Power(aMax,2)*jMax*(v0 + vf) + 2*Power(jMax,2)*(Power(v0,2) + Power(vf,2)))))*Abs(jMax);

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (6*Power(a0,2)*aMax*jMax - 18*Power(aMax,3)*jMax - 12*aMax*Power(jMax,2)*v0 + h1)/(12.*Power(aMax,2)*Power(jMax,2));
    profile.t[2] = aMax/jMax;
    profile.t[3] = 0;
    profile.t[4] = profile.t[2];
    profile.t[5] = (6*Power(af,2)*aMax*jMax - 18*Power(aMax,3)*jMax - 12*aMax*Power(jMax,2)*vf + h1)/(12.*Power(aMax,2)*Power(jMax,2));
    profile.t[6] = (af + aMax)/jMax;

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(pf, vf, af, vMax, aMax);
}

bool RuckigStep1::time_up_acc1(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
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
            return true;
        }
    }
    return false;
}

bool RuckigStep1::time_up_acc0(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
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
            return true;
        }
    }
    return false;
}

bool RuckigStep1::time_up_none(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        profile.t[0] = std::cbrt((pf - p0)/(2*jMax));
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        return profile.check(pf, vf, af, vMax, aMax);
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
        if (std::abs(Roots::polyEval(polynom, t)) > 1e-7) {
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
            return true;
        }
    }
    return false;
}

bool RuckigStep1::time_down_acc0_acc1_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, af, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc1_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    return time_up_acc1_vel(profile, p0, v0, a0, pf, vf, af, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc0_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    return time_up_acc0_vel(profile, p0, v0, a0, pf, vf, af, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    return time_up_vel(profile, p0, v0, a0, pf, vf, af, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc0_acc1(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, p0, v0, a0, pf, vf, af, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc1(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    return time_up_acc1(profile, p0, v0, a0, pf, vf, af, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc0(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    return time_up_acc0(profile, p0, v0, a0, pf, vf, af, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_none(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    return time_up_none(profile, p0, v0, a0, pf, vf, af, -vMax, -aMax, -jMax);
}

bool RuckigStep1::get_profile(Profile& profile, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax) {
    // Test all cases to get ones that match
    if (pf > p0) {
        if (time_up_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1_VEL;

        } else if (time_down_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1_VEL;

        } else if (time_up_acc1_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1_VEL;

        } else if (time_down_acc1_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1_VEL;

        } else if (time_up_acc0_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_VEL;

        } else if (time_down_acc0_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_VEL;

        } else if (time_up_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_VEL;

        } else if (time_down_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_VEL;

        } else if (time_up_none(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_NONE;

        } else if (time_up_acc0(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0;

        } else if (time_up_acc1(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1;

        } else if (time_up_acc0_acc1(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1;

        } else if (time_down_none(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_NONE;

        } else if (time_down_acc0(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0;

        } else if (time_down_acc1(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1;

        } else if (time_down_acc0_acc1(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1;

        } else {
            return false;
        }

    } else {
        if (time_down_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1_VEL;

        } else if (time_up_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1_VEL;

        } else if (time_down_acc1_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1_VEL;

        } else if (time_up_acc1_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1_VEL;

        } else if (time_down_acc0_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_VEL;

        } else if (time_up_acc0_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_VEL;

        } else if (time_down_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_VEL;

        } else if (time_up_vel(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_VEL;

        } else if (time_down_none(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_NONE;

        } else if (time_down_acc0(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0;

        } else if (time_down_acc1(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1;

        } else if (time_down_acc0_acc1(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1;

        } else if (time_up_none(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_NONE;

        } else if (time_up_acc0(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0;

        } else if (time_up_acc1(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1;

        } else if (time_up_acc0_acc1(profile, p0, v0, a0, pf, vf, af, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1;

        } else {
            return false;
        }
    }
    return true;
}

} // namespace ruckig
