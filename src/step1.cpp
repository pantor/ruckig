#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>


namespace ruckig {

Step1::Step1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double vMin, double aMax, double jMax): p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af), vMax(vMax), vMin(vMin), aMax(aMax), jMax(jMax) {
    pd = pf - p0;

    v0_v0 = v0 * v0;
    vf_vf = vf * vf;

    a0_a0 = a0 * a0;
    af_af = af * af;

    a0_p3 = a0 * a0_a0;
    a0_p4 = a0_a0 * a0_a0;
    a0_p5 = a0_p3 * a0_a0;
    a0_p6 = a0_p4 * a0_a0;
    af_p3 = af * af_af;
    af_p4 = af_af * af_af;
    af_p5 = af_p3 * af_af;
    af_p6 = af_p4 * af_af;

    // max values needs to be invariant to plus minus sign change
    aMax_aMax = aMax * aMax;
    jMax_jMax = jMax * jMax;
}

void Step1::time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax_aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*(a0_p4 + af_p4) + 8*aMax*(af_p3 - a0_p3) + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 6*af_af*(aMax_aMax - 2*jMax*vf) - 12*jMax*(aMax_aMax*(v0 + vf + 2*vMax) - jMax*(2*aMax*pd + v0_v0 + vf_vf - 2*vMax*vMax)))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = profile.t[2];
    profile.t[5] = (af_af/2 - aMax_aMax - jMax*(vf - vMax))/(aMax*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<Teeth::UDDU, Limits::ACC0_ACC1_VEL>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(a0_a0/2 + jMax*(vMax - v0))/Abs(jMax);

    profile.t[0] = (-a0 + h1*jMax)/jMax;
    profile.t[1] = 0;
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*af_p4 + 8*aMax*(af_p3 - a0_p3) + 24*aMax*jMax*(a0*v0 - af*vf) + 6*af_af*(aMax_aMax - 2*jMax*vf) - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(vf + vMax) + jMax*(-vf_vf + vMax*vMax) - aMax*h1*(a0_a0 - 2*jMax*(v0 + vMax))))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = aMax/jMax;
    profile.t[5] = (af_af/2 - aMax_aMax + jMax*(vMax - vf))/(aMax*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<Teeth::UDDU, Limits::ACC1_VEL>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(af_af/2 + jMax*(vMax - vf))/Abs(jMax);

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax_aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*a0_p4 + 8*(af_p3 - a0_p3)*aMax + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 12*af_af*aMax*h1*jMax - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vMax) + jMax*(-v0_v0 + vMax*vMax) + 2*aMax*(vf + vMax)*h1*jMax))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<Teeth::UDDU, Limits::ACC0_VEL>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_up_vel(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt(af_af/2 + jMax*(vMax - vf))/Abs(jMax);
    const double h2 = Sqrt(a0_a0/2 + jMax*(vMax - v0))/Abs(jMax);

    // Solution 3/4
    profile.t[0] = h2 - a0/jMax;
    profile.t[1] = 0;
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (af_p3 - a0_p3)/(3*jMax_jMax*vMax) + (a0*v0 - af*vf + (af_af*h1 + a0_a0*h2)/2)/(jMax*vMax) - (v0/vMax + 1.0)*h2 - (vf/vMax + 1.0)*h1 + pd/vMax;
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<Teeth::UDDU, Limits::VEL>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt((a0_p4 + af_p4)/2 + 4./3*(af_p3 - a0_p3)*aMax + 4*aMax*jMax*(a0*v0 - af*vf) + a0_a0*(aMax_aMax - 2*jMax*v0) + af_af*(aMax_aMax - 2*jMax*vf) + (aMax_aMax*(aMax_aMax - 2*jMax*(v0 + vf)) + 2*jMax_jMax*(v0_v0 + vf_vf + 2*aMax*pd)));
    
    if (!std::isnan(h1)) {
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

            if (profile.check<Teeth::UDDU, Limits::ACC0_ACC1>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, jMax);
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

            if (profile.check<Teeth::UDDU, Limits::ACC0_ACC1>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, jMax);
            }
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

        if (profile.check<Teeth::UDUD, Limits::ACC0_ACC1>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, jMax);
        }
    }
}

void Step1::time_up_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = 2*(2*a0 + aMax)/jMax;
    polynom[2] = (5*a0_a0 + 6*a0*aMax + aMax_aMax + 2*jMax*v0)/jMax_jMax;
    polynom[3] = (2*(a0 + aMax)*(a0_a0 + a0*aMax + 2*jMax*v0))/(jMax_jMax*jMax);
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
            
        if (profile.check<Teeth::UDDU, Limits::ACC1>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = 2*(2*a0 - aMax)/jMax;
        polynom[2] = (5*a0_a0 - 6*a0*aMax - aMax_aMax + 2*jMax*v0)/jMax_jMax;
        polynom[3] = (2*(a0_p3 - a0*aMax*(2*a0 + aMax) + 2*(a0 - aMax)*jMax*v0))/(jMax_jMax*jMax);
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
            profile.t[5] = -((a0_a0 - af_af)/2 + aMax_aMax + jMax*(2*a0*t + jMax*t*t - (vf - v0)))/(aMax*jMax);
            profile.t[6] = profile.t[4] - af/jMax;

            if (profile.check<Teeth::UDUD, Limits::ACC1>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, jMax);
            }
        }
    }
}

void Step1::time_up_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (2*aMax)/jMax;
    polynom[2] = (-af_af + aMax_aMax + 2*jMax*vf)/jMax_jMax;
    polynom[3] = (-2*aMax*(af_af - 2*jMax*vf))/(jMax_jMax*jMax);
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
            
        if (profile.check<Teeth::UDDU, Limits::ACC0>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, jMax);
        }
    }

    // UDUD
    if (std::abs(af) > DBL_EPSILON) {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (-2*aMax)/jMax;
        polynom[2] = -((af_af + aMax_aMax + 2*jMax*vf)/jMax_jMax);
        polynom[3] = (2*aMax*(af_af + 2*jMax*vf))/(jMax*jMax_jMax);
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
            
            if (profile.check<Teeth::UDUD, Limits::ACC0>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, jMax);
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

        if (profile.check<Teeth::UDDU, Limits::NONE>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, jMax);
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
    polynom[3] = -(a0_p5 + 8*(a0_a0 + jMax*v0)*p3 + a0_p3*(-6*af_af + 4*jMax*(v0 + 3*vf)) - 3*a0*p2)/(3*jMax_jMax*jMax*p1);
    polynom[4] = -(a0_p6 + af_p6 - 48*af_p3*jMax_jMax*pd + 144*af*jMax_jMax*jMax*pd*vf - 6*af_p4*jMax*(3*v0 + vf) + 16*a0_p3*p3 + 48*a0*jMax*v0*p3 - 36*af_af*jMax_jMax*(v0_v0 - 2*v0*vf - vf_vf) - 72*jMax_jMax*jMax*(jMax*pd*pd + (vf - v0)*(v0 + vf)*(v0 + vf)) + a0_p4*(-9*af_af + 6*jMax*(v0 + 3*vf)) - 9*a0_a0*p2)/(36*jMax_jMax*jMax_jMax*p1);

    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        // Single Newton step (regarding polynom)
        if (std::abs(Roots::polyEval(polynom, t)) > Roots::tolerance) {
            t -= Roots::polyEval(polynom, t) / Roots::polyEval(Roots::polyDeri(polynom), t);
        }

        // Single Newton step (regarding pd)
        {
            double h1 = a0_a0 + af_af + 4*a0*jMax*t + 2*jMax*(jMax*t*t + v0 - vf);
            double h2 = a0_a0 + 4*a0*jMax*t + 2*jMax*(jMax*t*t + v0);
            double orig = -pd + (-3*Sqrt(2)*jMax_jMax*jMax*std::pow(h1,1.5) - 6*af*jMax_jMax*h1*Abs(jMax) + 6*jMax*h2*Sqrt(2*h1)*jMax_jMax + 2*(2*a0_p3 + af_p3 + 3*a0_a0*(af + 4*jMax*t) + 6*af*jMax*(jMax*t*t + v0) + 6*a0*jMax*(2*af*t + 3*jMax*t*t + v0) + 6*jMax_jMax*(jMax*Power(t,3) + 2*t*v0))*Power(Abs(jMax),3))/(12*jMax_jMax*Power(Abs(jMax),3));
            double deriv = (2*a0*(a0 + af))/jMax + 6*a0*t + 2*af*t + 3*jMax*t*t + 2*v0 - (3*jMax_jMax*(a0 + jMax*t)*Sqrt(2*h1))/(2*Power(Abs(jMax),3)) - (2*af*jMax*(a0 + jMax*t))/jMax_jMax + (2*(a0 + jMax*t)*h2)/(Sqrt(2*h1)*Abs(jMax)) + (2*(a0 + jMax*t)*Sqrt(2*h1))/Abs(jMax);
            t -= orig / deriv;
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
 
        if (profile.check<Teeth::UDDU, Limits::NONE>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, jMax);
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
        polynom[3] = (40*a0_p3 - af_p3 - 3*jMax_jMax*pd - 3*af*jMax*vf + 3*a0*(af_af + 13*jMax*v0 + 2*jMax*vf))/(3*jMax_jMax*jMax);
        polynom[4] = (51*a0_p4 - af_p4 - 4*af_af*jMax*(vf - v0) + 2*a0_a0*(5*af_af + 58*jMax*v0 + 10*jMax*vf) - 8*a0*p1 + 4*jMax_jMax*(7*v0_v0 + 2*v0*vf - vf_vf))/(8*jMax_jMax*jMax_jMax);
        polynom[5] = (17*a0_p5 + 2*a0_p3*p3 - 8*(a0_a0 + jMax*v0)*p1 - 3*a0*p2)/(12*jMax_jMax*jMax_jMax*jMax);
        polynom[6] = -(-17*a0_p6 + af_p6 - 48*af_p3*jMax_jMax*pd - 144*af*jMax_jMax*jMax*pd*vf + 6*af_p4*jMax*(3*v0 + vf) - 3*a0_p4*p3 + 16*a0_p3*p1 + 48*a0*jMax*v0*p1 - 36*af_af*jMax_jMax*(v0_v0 - 2*v0*vf - vf_vf) - 72*jMax_jMax*jMax*(jMax*pd*pd - (vf - v0)*(v0 + vf)*(v0 + vf)) + 9*a0_a0*p2)/(144*jMax_jMax*jMax_jMax*jMax_jMax);

        std::array<double, 6> deriv = Roots::polyMonicDeri(polynom);
        std::array<double, 5> dderiv = Roots::polyMonicDeri(deriv);

        auto dd_extremas = Roots::solveQuartMonic(dderiv);
        std::set<std::pair<double, double>> dd_tz_intervals;

        double tz_min {0.0};
        double tz_max = {1e5};
        double dd_tz_current {tz_min};

        for (double tz: dd_extremas) {
            if (tz <= 0.0 || tz >= tz_max) {
                continue;
            }

            // if (std::abs(Roots::polyEval(dderiv, tz)) > 1e-14) {
            //     tz -= Roots::polyEval(dderiv, tz) / Roots::polyEval(Roots::polyDeri(dderiv), tz);
            // }

            if (Roots::polyEval(deriv, dd_tz_current) * Roots::polyEval(deriv, tz) < 0) {
                dd_tz_intervals.insert({dd_tz_current, tz});
            }
            dd_tz_current = tz;
        }
        if (Roots::polyEval(deriv, dd_tz_current) * Roots::polyEval(deriv, tz_max) < 0) {
            dd_tz_intervals.insert({dd_tz_current, tz_max});
        }

        std::set<double> roots;
        double tz_current {tz_min};

        for (auto interval: dd_tz_intervals) {
            double tz = Roots::shrinkInterval(deriv, interval.first, interval.second);

            if (tz <= 0.0) {
                continue;
            }

            if (Roots::polyEval(polynom, tz_current) * Roots::polyEval(polynom, tz) < 0) {
                roots.insert(Roots::shrinkInterval(polynom, tz_current, tz));
            }
            tz_current = tz;
        }
        if (Roots::polyEval(polynom, tz_current) * Roots::polyEval(polynom, tz_max) < 0) {
            roots.insert(Roots::shrinkInterval(polynom, tz_current, tz_max));
        }

        for (double t: roots) {
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

            if (profile.check<Teeth::UDUD, Limits::NONE>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, jMax);
            }
        }
    }
}

void Step1::time_down_acc0_acc1_vel(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, vMin, -aMax, -jMax);
}

void Step1::time_down_acc1_vel(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc1_vel(profile, vMin, -aMax, -jMax);
}

void Step1::time_down_acc0_vel(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc0_vel(profile, vMin, -aMax, -jMax);
}

void Step1::time_down_vel(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_vel(profile, vMin, -aMax, -jMax);
}

void Step1::time_down_acc0_acc1(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, vMin, -aMax, -jMax);
}

void Step1::time_down_acc1(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc1(profile, vMin, -aMax, -jMax);
}

void Step1::time_down_acc0(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc0(profile, vMin, -aMax, -jMax);
}

void Step1::time_down_none(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_none(profile, vMin, -aMax, -jMax);
}

bool Step1::calculate_block() {
    if (valid_profile_counter == 1) {
        block = Block(valid_profiles[0]);
        return true;
    }

    // Findex index of fastest profile
    auto idx_min_it = std::min_element(valid_profiles.cbegin(), valid_profiles.cbegin() + valid_profile_counter, [](auto& a, auto& b) { return a.t_sum[6] + a.t_brake.value_or(0.0) < b.t_sum[6] + b.t_brake.value_or(0.0); });
    size_t idx_min = std::distance(valid_profiles.cbegin(), idx_min_it);
    
    const double t_brake = valid_profiles[idx_min].t_brake.value_or(0.0); // t_brake is equal for every idx
    size_t count_udud = std::count_if(valid_profiles.cbegin(), valid_profiles.cbegin() + valid_profile_counter, [](auto& a) { return a.teeth == Profile::Teeth::UDUD; });
    // if (valid_profile_counter - count_udud > 3) {
    //     std::cout << valid_profile_counter - count_udud << std::endl;
    // }
    

    // std::array<Profile, 5> up_profiles, down_profiles;
    // size_t up_profile_counter {0}, down_profile_counter {0};
    // for (size_t i = 0; i < valid_profile_counter; ++i) {
    //     if (valid_profiles[i].direction == Profile::Direction::UP) {
    //         up_profiles[up_profile_counter] = valid_profiles[i];
    //         ++up_profile_counter;
    //     } else {
    //         down_profiles[down_profile_counter] = valid_profiles[i];
    //         ++down_profile_counter;
    //     }
    // }

    // std::stable_sort(up_profiles.begin(), up_profiles.begin() + up_profile_counter, [](auto& p1, auto& p2) { return p1.t_sum[6] < p2.t_sum[6]; });
    // std::stable_sort(down_profiles.begin(), down_profiles.begin() + down_profile_counter, [](auto& p1, auto& p2) { return p1.t_sum[6] < p2.t_sum[6]; });

    // size_t up_idx {0}, down_idx {0};


    /* if (down_profile_counter == 0 || (up_profile_counter > 0 && up_profiles[0].t_sum[6] < down_profiles[0].t_sum[6])) {
        size_t up_idx {0};
        if (up_profile_counter % 2 == 0) {
            ++up_idx;
        }

        block = Block(up_profiles[up_idx]);

        size_t inter {0};

        if (up_profile_counter >= 3) {
            block.a = Block::Interval(up_profiles[up_idx + 1].t_sum[6] + t_brake, up_profiles[up_idx + 2].t_sum[6] + t_brake, up_profiles[up_idx + 2]);
            ++inter;
        } else if (up_profile_counter >= 5) {
            block.b = Block::Interval(up_profiles[up_idx + 3].t_sum[6] + t_brake, up_profiles[up_idx + 4].t_sum[6] + t_brake, up_profiles[up_idx + 4]);
            ++inter;
        }

        size_t down_idx {0};
        if (down_profile_counter % 2 == 1) {
            std::cout << "---\n " << valid_profile_counter << " " << up_profile_counter << " " << down_profile_counter << std::endl;
            for (size_t i = 0; i < valid_profile_counter; ++i) {
                auto& p = valid_profiles[i];
                std::cout << p.t_sum[6] << " " << p.to_string() << std::endl;
            }

            ++down_idx;
        }

        if (down_profile_counter >= 2) {
            if (inter) {
                block.a = Block::Interval(down_profiles[down_idx + 1].t_sum[6] + t_brake, down_profiles[down_idx + 2].t_sum[6] + t_brake, down_profiles[down_idx + 2]);
            } else {
                block.b = Block::Interval(down_profiles[down_idx + 1].t_sum[6] + t_brake, down_profiles[down_idx + 2].t_sum[6] + t_brake, down_profiles[down_idx + 2]);
            }
        } else if (down_profile_counter >= 4) {
            block.b = Block::Interval(down_profiles[down_idx + 3].t_sum[6] + t_brake, down_profiles[down_idx + 4].t_sum[6] + t_brake, down_profiles[down_idx + 4]);
        }

        return true;
    } */

    // if (valid_profile_counter > 3) {
    //     std::cout << "---\n " << valid_profile_counter << std::endl;
    //     for (size_t i = 0; i < valid_profile_counter; ++i) {
    //         auto& p = valid_profiles[i];
    //         std::cout << p.t_sum[6] << " " << p.to_string() << std::endl;
    //     }
    // }

    if (valid_profile_counter == 2) {
        size_t idx_else_1 = (idx_min + 1) % 2;
        
        if (count_udud == 1) { // Ignore UDUD
            block = Block(valid_profiles[idx_min]);
        } else {
            block = Block(valid_profiles[idx_else_1]);
        }
        
    } else if (valid_profile_counter == 3) {
        size_t idx_else_1 = (idx_min + 1) % 3;
        size_t idx_else_2 = (idx_min + 2) % 3;

        block = Block(valid_profiles[idx_min]);
        add_interval(block.a, idx_else_1, idx_else_2, t_brake);

    } else if (valid_profile_counter == 4) {
        size_t idx_else_1 = (idx_min + 1) % 4;
        size_t idx_else_2 = (idx_min + 2) % 4;
        size_t idx_else_3 = (idx_min + 3) % 4;

        if (count_udud == 2) {
            if (valid_profiles[idx_min].teeth == valid_profiles[idx_else_1].teeth) {
                block = Block(valid_profiles[idx_else_1]);
                add_interval(block.a, idx_else_2, idx_else_3, t_brake);
            
            } else if (valid_profiles[idx_min].teeth == valid_profiles[idx_else_2].teeth) {
                block = Block(valid_profiles[idx_else_2]);
                add_interval(block.a, idx_else_1, idx_else_3, t_brake);
            
            } else {
                block = Block(valid_profiles[idx_else_3]);
                add_interval(block.a, idx_else_1, idx_else_2, t_brake);
            }
        
        } else if (count_udud == 1) { // Ignore UDUD
            block = Block(valid_profiles[idx_min]);
            if (valid_profiles[idx_else_1].teeth == Profile::Teeth::UDUD) {
                add_interval(block.a, idx_else_2, idx_else_3, t_brake);
        
            } else if (valid_profiles[idx_else_2].teeth == Profile::Teeth::UDUD) {
                add_interval(block.a, idx_else_1, idx_else_3, t_brake);
        
            } else {
                add_interval(block.a, idx_else_1, idx_else_2, t_brake);
            }
        
        } else {
            if (valid_profiles[idx_min].direction == valid_profiles[idx_else_1].direction) {
                block = Block(valid_profiles[idx_else_1]);
                add_interval(block.a, idx_else_2, idx_else_3, t_brake);
            
            } else if (valid_profiles[idx_min].direction == valid_profiles[idx_else_2].direction) {
                block = Block(valid_profiles[idx_else_2]);
                add_interval(block.a, idx_else_1, idx_else_3, t_brake);
            
            } else {
                block = Block(valid_profiles[idx_else_3]);
                add_interval(block.a, idx_else_1, idx_else_2, t_brake);
            }
        }
        
    } else if (valid_profile_counter == 5) {
        size_t idx_else_1 = (idx_min + 1) % 5;
        size_t idx_else_2 = (idx_min + 2) % 5;
        size_t idx_else_3 = (idx_min + 3) % 5;
        size_t idx_else_4 = (idx_min + 4) % 5;

        block = Block(valid_profiles[idx_min]);
        if (valid_profiles[idx_else_1].direction == valid_profiles[idx_else_2].direction) {
            add_interval(block.a, idx_else_1, idx_else_2, t_brake);
            add_interval(block.b, idx_else_3, idx_else_4, t_brake);
        } else {
            add_interval(block.a, idx_else_1, idx_else_4, t_brake);
            add_interval(block.b, idx_else_2, idx_else_3, t_brake);
        }

    } else if (valid_profile_counter == 6) {
        size_t idx_else_1 = (idx_min + 1) % 6;
        size_t idx_else_2 = (idx_min + 2) % 6;
        size_t idx_else_3 = (idx_min + 3) % 6;
        size_t idx_else_4 = (idx_min + 4) % 6;
        size_t idx_else_5 = (idx_min + 5) % 6;

        block = Block(valid_profiles[idx_else_1]);
        if (valid_profiles[idx_else_2].direction == valid_profiles[idx_else_3].direction) {
            add_interval(block.a, idx_else_2, idx_else_3, t_brake);
            add_interval(block.b, idx_else_4, idx_else_5, t_brake);
        } else {
            add_interval(block.a, idx_else_2, idx_else_5, t_brake);
            add_interval(block.b, idx_else_3, idx_else_4, t_brake);
        }

    } else {
        return false;
    }

    return true;
}

bool Step1::get_profile(const Profile& input) {
    Profile profile = input;
    profile.a[0] = a0;
    profile.v[0] = v0;
    profile.p[0] = p0;
    valid_profile_counter = 0;

    if (pf > p0) {
        time_up_acc0_acc1_vel(profile, vMax, aMax, jMax);
        time_down_acc0_acc1_vel(profile, vMin, aMax, jMax);
        time_up_acc1_vel(profile, vMax, aMax, jMax);
        time_down_acc1_vel(profile, vMin, aMax, jMax);
        time_up_acc0_vel(profile, vMax, aMax, jMax);
        time_down_acc0_vel(profile, vMin, aMax, jMax);
        time_up_vel(profile, vMax, aMax, jMax);
        time_down_vel(profile, vMin, aMax, jMax);
        time_up_none(profile, vMax, aMax, jMax);
        time_up_acc0(profile, vMax, aMax, jMax);
        time_up_acc1(profile, vMax, aMax, jMax);
        time_up_acc0_acc1(profile, vMax, aMax, jMax);
        time_down_none(profile, vMin, aMax, jMax);
        time_down_acc0(profile, vMin, aMax, jMax);
        time_down_acc1(profile, vMin, aMax, jMax);
        time_down_acc0_acc1(profile, vMin, aMax, jMax);

    } else {
        time_down_acc0_acc1_vel(profile, vMin, aMax, jMax);
        time_up_acc0_acc1_vel(profile, vMax, aMax, jMax);
        time_down_acc1_vel(profile, vMin, aMax, jMax);
        time_up_acc1_vel(profile, vMax, aMax, jMax);
        time_down_acc0_vel(profile, vMin, aMax, jMax);
        time_up_acc0_vel(profile, vMax, aMax, jMax);
        time_down_vel(profile, vMin, aMax, jMax);
        time_up_vel(profile, vMax, aMax, jMax);
        time_down_none(profile, vMin, aMax, jMax);
        time_down_acc0(profile, vMax, aMax, jMax);
        time_down_acc1(profile, vMax, aMax, jMax);
        time_down_acc0_acc1(profile, vMax, aMax, jMax);
        time_up_none(profile, vMin, aMax, jMax);
        time_up_acc0(profile, vMin, aMax, jMax);
        time_up_acc1(profile, vMin, aMax, jMax);
        time_up_acc0_acc1(profile, vMin, aMax, jMax);
    }

    return valid_profile_counter && calculate_block();
}

} // namespace ruckig
