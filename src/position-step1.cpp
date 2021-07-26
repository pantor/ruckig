#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>


namespace ruckig {

PositionStep1::PositionStep1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double vMin, double aMax, double aMin, double jMax): p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af), _vMax(vMax), _vMin(vMin), _aMax(aMax), _aMin(aMin), _jMax(jMax) {
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
    jMax_jMax = jMax * jMax;
}

void PositionStep1::time_acc0_acc1_vel(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax*aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = aMax/jMax;
    profile.t[3] = (3*(a0_p4*aMin - af_p4*aMax) + 8*aMax*aMin*(af_p3 - a0_p3 + 3*jMax*(a0*v0 - af*vf)) + 6*a0_a0*aMin*(aMax*aMax - 2*jMax*v0) - 6*af_af*aMax*(aMin*aMin - 2*jMax*vf) - 12*jMax*(aMax*aMin*(aMax*(v0 + vMax) - aMin*(vf + vMax) - 2*jMax*pd) + (aMin - aMax)*jMax*vMax*vMax + jMax*(aMax*vf_vf - aMin*v0_v0)))/(24*aMax*aMin*jMax_jMax*vMax);
    profile.t[4] = -aMin/jMax;
    profile.t[5] = -(af_af/2 - aMin*aMin - jMax*(vf - vMax))/(aMin*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1_VEL>(jMax, vMax, vMin, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void PositionStep1::time_acc1_vel(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }

    const double h1 = Sqrt(a0_a0/(2*jMax_jMax) + (vMax - v0)/jMax);

    profile.t[0] = -a0/jMax + h1;
    profile.t[1] = 0;
    profile.t[2] = h1;
    profile.t[3] = -(3*af_p4 - 8*aMin*(af_p3 - a0_p3) - 24*aMin*jMax*(a0*v0 - af*vf) + 6*af_af*(aMin*aMin - 2*jMax*vf) - 12*jMax*(2*aMin*jMax*pd + aMin*aMin*(vf + vMax) + jMax*(vMax*vMax - vf_vf) + aMin*h1*(a0_a0 - 2*jMax*(v0 + vMax))))/(24*aMin*jMax_jMax*vMax);
    profile.t[4] = -aMin/jMax;
    profile.t[5] = -(af_af/2 - aMin*aMin + jMax*(vMax - vf))/(aMin*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (std::abs(profile.t[0]) < DBL_EPSILON && profile.t[0] < 0.0) {
        profile.t[0] = 0;
    }

    if (profile.check<JerkSigns::UDDU, Limits::ACC1_VEL>(jMax, vMax, vMin, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void PositionStep1::time_acc0_vel(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }

    const double h1 = Sqrt(af_af/(2*jMax_jMax) + (vMax - vf)/jMax);

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax*aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = aMax/jMax;
    profile.t[3] = (3*a0_p4 + 8*(af_p3 - a0_p3)*aMax + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax*aMax - 2*jMax*v0) - 12*jMax*(-2*aMax*jMax*pd + aMax*aMax*(v0 + vMax) + jMax*(vMax*vMax - v0_v0) + (2*(vf + vMax)*jMax - af_af)*aMax*h1))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::ACC0_VEL>(jMax, vMax, vMin, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void PositionStep1::time_vel(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }

    const double h1 = Sqrt(af_af/(2*jMax_jMax) + (vMax - vf)/jMax);
    const double h2 = Sqrt(a0_a0/(2*jMax_jMax) + (vMax - v0)/jMax);

    // Solution 3/4
    profile.t[0] = h2 - a0/jMax;
    profile.t[1] = 0;
    profile.t[2] = h2;
    profile.t[3] = (af_p3 - a0_p3)/(3*jMax_jMax*vMax) + (a0*v0 - af*vf + (af_af*h1 + a0_a0*h2)/2)/(jMax*vMax) - (v0/vMax + 1.0)*h2 - (vf/vMax + 1.0)*h1 + pd/vMax;
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = h1 + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::VEL>(jMax, vMax, vMin, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void PositionStep1::time_acc0_acc1(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    const double h1 = Sqrt(3*(aMax - aMin)*(3*(af_p4*aMax - a0_p4*aMin) + 8*(a0_p3 - af_p3)*aMax*aMin + 24*aMax*aMin*jMax*(af*vf - a0*v0) + 6*af_af*aMax*(aMin*aMin - 2*jMax*vf) - 6*a0_a0*aMin*(aMax*aMax - 2*jMax*v0) + 3*(aMax*aMax*aMax*aMin*aMin - 4*aMin*jMax_jMax*v0_v0 - aMax*aMax*aMin*(aMin*aMin - 4*jMax*v0) + 4*aMax*jMax*(-2*aMin*jMax*pd - aMin*aMin*vf + jMax*vf_vf))))*Abs(jMax)/(3*(aMax - aMin)*jMax);

    if (!std::isnan(h1)) {
        // UDDU: Solution 2
        {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (a0_a0 + aMax*aMin - 2*(aMax*aMax + jMax*v0) - h1)/(2*aMax*jMax);
            profile.t[2] = aMax/jMax;
            profile.t[3] = 0;
            profile.t[4] = -aMin/jMax;
            profile.t[5] = -(af_af + aMax*aMin - 2*(aMin*aMin + jMax*vf) - h1)/(2*aMin*jMax);
            profile.t[6] = profile.t[4] + af/jMax;

            if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1>(jMax, vMax, vMin, aMax, aMin)) {
                add_profile(profile, jMax);
            }
        }

        // UDDU: Solution 1
        {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (a0_a0 + aMax*aMin - 2*(aMax*aMax + jMax*v0) + h1)/(2*aMax*jMax);
            profile.t[2] = aMax/jMax;
            profile.t[3] = 0;
            profile.t[4] = -aMin/jMax;
            profile.t[5] = -(af_af + aMax*aMin - 2*(aMin*aMin + jMax*vf) + h1)/(2*aMin*jMax);
            profile.t[6] = profile.t[4] + af/jMax;

            if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1>(jMax, vMax, vMin, aMax, aMin)) {
                add_profile(profile, jMax);
            }
        }
    }
}

void PositionStep1::time_acc1(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (2*(2*a0 - aMin))/jMax;
    polynom[2] = (5*a0_a0 - 6*a0*aMin + aMin*aMin + 2*jMax*v0)/jMax_jMax;
    polynom[3] = (2*(a0 - aMin)*(a0_a0 - a0*aMin + 2*jMax*v0))/(jMax_jMax*jMax);
    polynom[4] = (3*(a0_p4 - af_p4) - 8*(a0_p3 - af_p3)*aMin - 24*aMin*jMax*(a0*v0 + af*vf) + 6*a0_a0*(aMin*aMin + 2*jMax*v0) - 6*af_af*(aMin*aMin - 2*jMax*vf) + 12*jMax*(2*aMin*jMax*pd + aMin*aMin*(v0 + vf) + jMax*(v0_v0 - vf_vf)))/(12*jMax_jMax*jMax_jMax);

    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        // Double Newton step (regarding pd)
        if (t > DBL_EPSILON) {
            double h1 = jMax*t*t + v0;
            double orig = -pd + (3*(a0_p4 - af_p4) + 8*af_p3*aMin + 8*a0_p3*(-aMin + 3*jMax*t) + 24*a0*jMax*(-aMin + 2*jMax*t)*(-aMin*t + h1) + 6*a0_a0*(aMin*aMin + 2*jMax*(-4*aMin*t + 5*h1 - 4*v0)) - 24*af*aMin*jMax*vf - 6*af_af*(aMin*aMin - 2*jMax*vf) + 12*jMax*(-2*aMin*jMax*t*(h1 + v0) + aMin*aMin*(h1 + vf) + jMax*(h1*h1 - vf_vf)))/(-24*aMin*jMax_jMax);
            double deriv = -((a0 - aMin + jMax*t)*(a0_a0 - aMin*jMax*t + a0*(-aMin + 4*jMax*t) + 2*jMax*h1))/(aMin*jMax);

            t -= std::min(orig / deriv, t);

            h1 = jMax*t*t + v0;
            orig = -pd + (3*(a0_p4 - af_p4) + 8*af_p3*aMin + 8*a0_p3*(-aMin + 3*jMax*t) + 24*a0*jMax*(-aMin + 2*jMax*t)*(-aMin*t + h1) + 6*a0_a0*(aMin*aMin + 2*jMax*(-4*aMin*t + 5*h1 - 4*v0)) - 24*af*aMin*jMax*vf - 6*af_af*(aMin*aMin - 2*jMax*vf) + 12*jMax*(-2*aMin*jMax*t*(h1 + v0) + aMin*aMin*(h1 + vf) + jMax*(h1*h1 - vf_vf)))/(-24*aMin*jMax_jMax);
            if (std::abs(orig) > 1e-9) {
                deriv = -((a0 - aMin + jMax*t)*(a0_a0 - aMin*jMax*t + a0*(-aMin + 4*jMax*t) + 2*jMax*h1))/(aMin*jMax);

                t -= orig / deriv;
            }
        }

        profile.t[0] = t;
        profile.t[1] = 0;
        profile.t[2] = (a0 - aMin)/jMax + t;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = -((a0_a0 + af_af)/2 - aMin*aMin + jMax*t*(2*a0 + jMax*t) - jMax*(vf - v0))/(aMin*jMax);
        profile.t[6] = (af - aMin)/jMax;

        if (profile.check<JerkSigns::UDDU, Limits::ACC1>(jMax, vMax, vMin, aMax, aMin)) {
            add_profile(profile, jMax);
        }
    }
}

void PositionStep1::time_acc0(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    // Combined UDDU and UDUD
    // UDUD Strategy t7 == 0 => is equal to UDDU
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (-2*aMax)/jMax;
    polynom[2] = (-af_af + aMax*aMax + 2*jMax*vf)/jMax_jMax;
    polynom[3] = 0;
    polynom[4] = (3*(af_p4 - a0_p4) + 8*(a0_p3 - af_p3)*aMax + 24*aMax*jMax*(af*vf - a0*v0) - 6*a0_a0*(aMax*aMax - 2*jMax*v0) + 6*af_af*(aMax*aMax - 2*jMax*vf) + 12*jMax*(jMax*(vf_vf - v0_v0 - 2*aMax*pd) - aMax*aMax*(vf - v0)))/(12*jMax_jMax*jMax_jMax);

    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        // Single Newton step (regarding pd)
        if (t > DBL_EPSILON) {
            double orig = (3*(af_p4 - a0_p4) + 8*(a0_p3 - af_p3)*aMax + 24*aMax*jMax*(af*vf - a0*v0) - 6*a0_a0*(aMax*aMax - 2*jMax*v0) + 6*af_af*(aMax*aMax - 2*jMax*(jMax*t*t + vf)) + 12*jMax*(-2*aMax*jMax*(pd + jMax*t*t*t) + aMax*aMax*(jMax*t*t + v0 - vf) + jMax*(jMax*t*t*(jMax*t*t + 2*vf) + vf_vf - v0_v0)))/(24*aMax*jMax_jMax);
            double deriv = t*(-af_af + aMax*aMax - 3*aMax*jMax*t + 2*jMax*(jMax*t*t + vf))/aMax;

            t -= orig / deriv;
        }

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (a0_a0 - af_af + 2*jMax*((-2*aMax + jMax*t)*t + vf - v0))/(2*aMax*jMax);
        profile.t[2] = t;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = (af - aMax)/jMax + t;

        if (profile.check<JerkSigns::UDDU, Limits::ACC0>(jMax, vMax, vMin, aMax, aMin)) {
            add_profile(profile, jMax);
        }
    }
}

void PositionStep1::time_none(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    // UDDU
    {
        if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
            profile.t[0] = std::cbrt(pd/(2*jMax));
            profile.t[1] = 0;
            profile.t[2] = 2*profile.t[0];
            profile.t[3] = 0;
            profile.t[4] = 0;
            profile.t[5] = 0;
            profile.t[6] = profile.t[0];

            if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, vMax, vMin, aMax, aMin)) {
                add_profile(profile, jMax);
            }
            return;
        }

        if (std::abs(a0 - af) < DBL_EPSILON && std::abs(v0 + vf) < DBL_EPSILON && std::abs(p0 - pf) < DBL_EPSILON) {
            const double h1 = std::sqrt(a0_a0 - 2*jMax*v0);

            // Solution 3
            {
                profile.t[0] = -(a0 + h1)/jMax;
                profile.t[1] = 0;
                profile.t[2] = profile.t[0];
                profile.t[3] = 0;
                profile.t[4] = 0;
                profile.t[5] = 0;
                profile.t[6] = 0;

                if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, vMax, vMin, aMax, aMin)) {
                    add_profile(profile, jMax);
                }
            }

            // Solution 4
            {
                profile.t[0] = -(a0 - h1)/jMax;
                profile.t[1] = 0;
                profile.t[2] = profile.t[0];
                profile.t[3] = 0;
                profile.t[4] = 0;
                profile.t[5] = 0;
                profile.t[6] = 0;

                if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, vMax, vMin, aMax, aMin)) {
                    add_profile(profile, jMax);
                }
            }

            return;
        }
    }

    // UDDU / UDUD modern, this one is in particular prone to numerical issues
    {
        const double h2 = (a0_a0 - af_af + 2*jMax*(vf - v0))/(2*jMax_jMax);
        
        // UDUD Strategy: t7 == 0 (equals UDDU)
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = 0;
        polynom[2] = (-2*(a0_a0 + af_af - 2*jMax*(v0 + vf)))/jMax_jMax;
        polynom[3] = (4*(a0_p3 - af_p3 + 3*jMax*(af*vf - a0*v0 - jMax*pd)))/(3*jMax*jMax_jMax);
        polynom[4] = -h2*h2;

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0) {
                continue;
            }

            // Single Newton-step (regarding pd)
            {
                const double h1 = (a0_a0 - af_af)/(2*jMax) + (vf - v0);
                const double orig = (-h1*h1 + 4*h1*t*(af + jMax*t))/(4*jMax*t) + (4*a0_p3 + 2*af_p3 - 6*a0_a0*(af + 2*jMax*t) + 12*(af - a0)*jMax*v0 + 3*jMax_jMax*(-4*pd + jMax*t*t*t + 8*t*v0))/(12*jMax_jMax);
                const double deriv = h1 + 2*v0 - a0_a0/jMax + h1*h1/(4*jMax*t*t) + (3*jMax*t*t)/4;
                
                t -= orig / deriv;
            }

            const double h0 = ((a0_a0 - af_af)/jMax + 2*(vf - v0))/(4*jMax*t);
            profile.t[0] = h0 + t/2 - a0/jMax;
            profile.t[1] = 0;
            profile.t[2] = t;
            profile.t[3] = 0;
            profile.t[4] = 0;
            profile.t[5] = 0;
            profile.t[6] = -h0 + t/2 + af/jMax;

            if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, vMax, vMin, aMax, aMin)) {
                add_profile(profile, jMax);
            }
        }
    }
}


void PositionStep1::time_acc1_vel_two_step(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    profile.t[0] = 0;
    profile.t[1] = 0;
    profile.t[2] = a0/jMax;
    profile.t[3] = -(3*af_p4 - 8*aMin*(af_p3 - a0_p3) - 24*aMin*jMax*(a0*v0 - af*vf) + 6*af_af*(aMin*aMin - 2*jMax*vf) - 12*jMax*(2*aMin*jMax*pd + aMin*aMin*(vf + vMax) + jMax*(vMax*vMax - vf_vf) + aMin*a0*(a0_a0 - 2*jMax*(v0 + vMax))/jMax))/(24*aMin*jMax_jMax*vMax);
    profile.t[4] = -aMin/jMax;
    profile.t[5] = -(af_af/2 - aMin*aMin + jMax*(vMax - vf))/(aMin*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::ACC1_VEL>(jMax, vMax, vMin, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void PositionStep1::time_acc0_two_step(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    // Two step
    {
        profile.t[0] = 0;
        profile.t[1] = (af_af - a0_a0 + 2*jMax*(vf - v0))/(2*a0*jMax);
        profile.t[2] = (a0 - af)/jMax;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::ACC0>(jMax, vMax, vMin, aMax, aMin)) {
            add_profile(profile, jMax);
            return;
        }
    }

    // Three step
    {
        const double h0 = 3*(af_af - a0_a0 + 2*jMax*(v0 + vf));
        const double h1 = Sqrt(2*(2*Power2(a0_p3 + 2*af_p3 + 6*jMax_jMax*pd + 6*(af - a0)*jMax*vf - 3*a0*af_af) + h0*(a0_p4 - 6*a0_a0*(af_af + 2*jMax*vf) + 8*a0*(af_p3 + 3*jMax_jMax*pd + 3*af*jMax*vf) - 3*(af_p4 + 4*af_af*jMax*vf + 4*jMax_jMax*(vf_vf - v0_v0))))) * Abs(jMax) / jMax;
        profile.t[0] = (4*af_p3 + 2*a0_p3 - 6*a0*af_af + 12*jMax_jMax*pd + 12*(af - a0)*jMax*vf + h1)/(2*jMax*h0);
        profile.t[1] = -h1/(jMax*h0);
        profile.t[2] = (-4*a0_p3 - 2*af_p3 + 6*a0_a0*af + 12*jMax_jMax*pd - 12*(af - a0)*jMax*v0 + h1)/(2*jMax*h0);
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::ACC0>(jMax, vMax, vMin, aMax, aMin)) {
            add_profile(profile, jMax);
            return;
        }
    }
}

void PositionStep1::time_vel_two_step(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    // Four step
    {
        const double h1 = Sqrt(af_af/(2*jMax_jMax) + (vMax - vf)/jMax);

        // Solution 3/4
        profile.t[0] = -a0/jMax;
        profile.t[1] = 0;
        profile.t[2] = 0;
        profile.t[3] = (af_p3 - a0_p3)/(3*jMax_jMax*vMax) + (a0*v0 - af*vf + (af_af*h1)/2)/(jMax*vMax) - (vf/vMax + 1.0)*h1 + pd/vMax;
        profile.t[4] = h1;
        profile.t[5] = 0;
        profile.t[6] = h1 + af/jMax;

        if (profile.check<JerkSigns::UDDU, Limits::VEL>(jMax, vMax, vMin, aMax, aMin)) {
            add_profile(profile, jMax);
            return;
        }
    }


    // Four step
    {
        const double h1 = Sqrt(af_af/(2*jMax_jMax) + (vMax - vf)/jMax);

        profile.t[0] = 0;
        profile.t[1] = 0;
        profile.t[2] = a0/jMax;
        profile.t[3] = (af_p3 - a0_p3)/(3*jMax_jMax*vMax) + (a0*v0 - af*vf + (af_af*h1 + a0_p3/jMax)/2)/(jMax*vMax) - (v0/vMax + 1.0)*a0/jMax - (vf/vMax + 1.0)*h1 + pd/vMax;
        profile.t[4] = h1;
        profile.t[5] = 0;
        profile.t[6] = h1 + af/jMax;

        if (profile.check<JerkSigns::UDDU, Limits::VEL>(jMax, vMax, vMin, aMax, aMin)) {
            add_profile(profile, jMax);
            return;
        }
    }
}

void PositionStep1::time_none_two_step(Profile& profile, double vMax, double vMin, double aMax, double aMin, double jMax) {
    // Two step
    {
        const double h0 = Sqrt((a0_a0 + af_af)/2 + jMax*(vf - v0)) * Abs(jMax) / jMax;
        profile.t[0] = (h0 - a0)/jMax;
        profile.t[1] = 0;
        profile.t[2] = (h0 - af)/jMax;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, vMax, vMin, aMax, aMin)) {
            add_profile(profile, jMax);
            return;
        }
    }

    // Single step
    {
        profile.t[0] = (af - a0)/jMax;
        profile.t[1] = 0;
        profile.t[2] = 0;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, vMax, vMin, aMax, aMin)) {
            add_profile(profile, jMax);
            return;
        }
    }
}


bool PositionStep1::get_profile(const Profile& input, Block& block) {
    Profile profile = input;
    profile.set_boundary(p0, v0, a0, pf, vf, af);
    valid_profile_counter = 0;

    if (std::abs(pf - p0) < DBL_EPSILON && std::abs(v0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        if (pf >= p0) {
            time_none(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        } else {
            time_none(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        }

    } else {
        time_acc0_acc1_vel(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        time_acc0_acc1_vel(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        time_acc1_vel(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        time_acc1_vel(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        time_acc0_vel(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        time_acc0_vel(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        time_vel(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        time_vel(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        time_none(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        time_acc0(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        time_acc1(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        time_acc0_acc1(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        time_none(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        time_acc0(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        time_acc1(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        time_acc0_acc1(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);

        if (valid_profile_counter == 0 || valid_profile_counter == 2 || valid_profile_counter == 4) {
            time_none_two_step(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        }
        if (valid_profile_counter == 0 || valid_profile_counter == 2 || valid_profile_counter == 4) {
            time_none_two_step(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        }
        if (valid_profile_counter == 0 || valid_profile_counter == 2) {
            time_acc0_two_step(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        }
        if (valid_profile_counter == 0 || valid_profile_counter == 2) {
            time_acc0_two_step(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        }
        if (valid_profile_counter == 0 || valid_profile_counter == 2) {
            time_vel_two_step(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        }
        if (valid_profile_counter == 0 || valid_profile_counter == 2) {
            time_vel_two_step(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        }
        if (valid_profile_counter == 0 || valid_profile_counter == 2) {
            time_acc1_vel_two_step(profile, _vMax, _vMin, _aMax, _aMin, _jMax);
        }
        if (valid_profile_counter == 0 || valid_profile_counter == 2) {
            time_acc1_vel_two_step(profile, _vMin, _vMax, _aMin, _aMax, -_jMax);
        }
    }

    return Block::calculate_block(block, valid_profiles, valid_profile_counter);
}

} // namespace ruckig
