#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>


namespace ruckig {

Step2::Step2(double tf, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double vMin, double aMax, double jMax): tf(tf), p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af), vMax(vMax), vMin(vMin), aMax(aMax), jMax(jMax)  {
    pd = pf - p0;
    tf_tf = tf * tf;
    tf_p3 = tf_tf * tf;
    tf_p4 = tf_tf * tf_tf;

    vd = vf - v0;
    vd_vd = vd * vd;
    v0_v0 = v0 * v0;
    vf_vf = vf * vf;

    ad = af - a0;
    ad_ad = ad * ad;
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

bool Step2::time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    if (tf < 2*aMax/jMax) {
        return false;
    }

    // Profile UDDU
    {
        const double h0b = aMax_aMax + jMax*(vd - aMax*tf);
        const double h0a = af_af + 2*(af*aMax + aMax_aMax - aMax*jMax*tf - jMax*vd);
        const double h1 = Sqrt(-a0_p4 - af_p4 + 4./3*aMax*(a0_p3 - af_p3) + 2*a0*h0a*(a0 - 2*aMax) + 4*af*h0b*(af + 2*aMax) + 4*(aMax_aMax*aMax_aMax - 2*aMax_aMax*aMax*jMax*tf + aMax_aMax*jMax_jMax*tf_tf - jMax_jMax*vd_vd + 2*aMax*jMax_jMax*(-2*pd + tf*(v0 + vf))));
        const double h2 = 2*aMax*(ad + 3*aMax - jMax*tf) + h1;
        const double h3 = 4*aMax*jMax;
        const double h4 = a0_a0 - af_af + 2*jMax*vd;

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = -(h2 - h4)/h3;
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = tf - (ad + 4*aMax)/jMax + 2*h2/h3;
        profile.t[4] = profile.t[2];
        profile.t[5] = -(h2 + h4)/h3;
        profile.t[6] = profile.t[4] + af/jMax;

        if (profile.check<Teeth::UDDU, Limits::ACC0_ACC1_VEL>(tf, pf, vf, af, jMax, vMax, aMax)) {
            return true;
        }
    }

    // Profile UDUD
    {
        const double h1 = 12*aMax*jMax*(a0_a0 + af_af - 2*(a0 + af)*aMax + 2*(aMax_aMax - aMax*jMax*tf + jMax*vd));
        const double h2 = 3*(a0_p4 + af_p4) - 4*(a0_p3 + af_p3)*aMax;
        const double h3 = -4*af_p3*aMax + 24*(a0 + af)*aMax_aMax*aMax - 6*(af_af + a0_a0)*(aMax_aMax - 2*jMax*vd) + 6*a0_a0*(af_af - 2*af*aMax - 2*aMax*jMax*tf) - 12*aMax_aMax*(2*aMax_aMax - 2*aMax*jMax*tf + jMax*vd) - 24*af*aMax*jMax*vd + 12*jMax_jMax*(2*aMax*(-pd + tf*v0) + vd_vd);

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (h2 + h3)/h1;
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = -(a0_a0 + af_af - 2*aMax*(a0 + af + jMax*tf) + 4*aMax_aMax + 2*jMax*vd)/(2*aMax*jMax);
        profile.t[4] = profile.t[2];
        profile.t[5] = tf - (profile.t[0] + profile.t[1] + profile.t[2] + profile.t[3] + 2*profile.t[4] - af/jMax);
        profile.t[6] = profile.t[4] - af/jMax;

        if (profile.check<Teeth::UDUD, Limits::ACC0_ACC1_VEL>(tf, pf, vf, af, jMax, vMax, aMax)) {
            return true;
        }
    }
    
    return false;
}

bool Step2::time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    if (tf < aMax/jMax) {
        return false;
    }
    
    // Profile UDDU
    {
        const double ph1 = a0_a0 + af_af + aMax*(a0 + 2*af) + aMax_aMax - 2*jMax*(vd + aMax*tf);
        const double ph2 = -2*aMax*jMax*(-pd + tf*v0) - aMax_aMax*vd + jMax*vd_vd;
        const double ph3 = af_af + 2*af*aMax + aMax_aMax - 2*jMax*(vd + aMax*tf);

        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (2*(2*a0 + aMax))/jMax;
        polynom[2] = (4*a0_a0 + ph1 + 3*a0*aMax)/jMax_jMax;
        polynom[3] = (2*a0*ph1)/(jMax_jMax*jMax);
        polynom[4] = (3*(a0_p4 + af_p4) + 4*(a0_p3 + 2*af_p3)*aMax + 6*af_af*(aMax_aMax - 2*jMax*vd) + 12*jMax*ph2 - 24*af*aMax*jMax*vd + 6*a0_a0*ph3)/(12*jMax_jMax*jMax_jMax);
        auto roots = Roots::solveQuartMonic(polynom);

        for (double t: roots) {
            if (t < 0.0 || t > tf - aMax/jMax) {
                continue;
            }

            const double h1 = ((a0_a0 + af_af)/2 + jMax*(v0 - vf + 2*a0*t + jMax*t*t))/aMax;
            
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = tf - (h1 + aMax + a0 + af)/jMax - 2*t;
            profile.t[4] = aMax/jMax;
            profile.t[5] = (h1 - aMax)/jMax;
            profile.t[6] = profile.t[4] + af/jMax;
            
            if (profile.check<Teeth::UDDU, Limits::ACC1_VEL>(tf, pf, vf, af, jMax, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        const double ph1 = a0_a0 - af_af + (2*af - a0)*aMax - aMax_aMax - 2*jMax*(vd - aMax*tf);
        const double ph2 = aMax_aMax + 2*jMax*vd;
        const double ph3 = af_af + ph2 - 2*aMax*(af + jMax*tf);
        const double ph4 = 2*aMax*jMax*(-pd + tf*v0) + aMax_aMax*vd + jMax*vd_vd;

        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (4*a0 - 2*aMax)/jMax;
        polynom[2] = (4*a0_a0 - 3*a0*aMax + ph1)/jMax_jMax;
        polynom[3] = (2*a0*ph1)/(jMax_jMax*jMax);
        polynom[4] = (3*(a0_p4 + af_p4) - 4*(a0_p3 + 2*af_p3)*aMax - 24*af*aMax*jMax*vd + 12*jMax*ph4 - 6*a0_a0*ph3 + 6*af_af*ph2)/(12*jMax_jMax*jMax_jMax);

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0 || t > tf - aMax/jMax) {
                continue;
            }

            const double h1 = ((a0_a0 - af_af)/2 + jMax_jMax*t*t - jMax*(vd - 2*a0*t))/aMax;

            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + t;
            profile.t[3] = tf + (h1 + ad - aMax)/jMax - 2*t;
            profile.t[4] = aMax/jMax;
            profile.t[5] = -(h1 + aMax)/jMax;
            profile.t[6] = profile.t[4] - af/jMax;
            
            if (profile.check<Teeth::UDUD, Limits::ACC1_VEL>(tf, pf, vf, af, jMax, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool Step2::time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    if (tf < aMax/jMax) {
        return false;
    }

    const double ph1 = 12*jMax*(-aMax_aMax*vd - jMax*vd_vd + 2*aMax*jMax*(-pd + tf*vf));
       
    // Profile UDDU
    {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (2*aMax)/jMax;
        polynom[2] = (a0_a0 - af_af + 2*ad*aMax + aMax_aMax + 2*jMax*(vd - aMax*tf))/jMax_jMax;
        polynom[3] = 0;
        polynom[4] = -(-3*(a0_p4 + af_p4) + 4*(af_p3 + 2*a0_p3)*aMax - 12*a0*aMax*(af_af - 2*jMax*vd) + 6*a0_a0*(af_af - aMax_aMax - 2*jMax*vd) + 6*af_af*(aMax_aMax - 2*aMax*jMax*tf + 2*jMax*vd) + ph1)/(12*jMax_jMax*jMax_jMax);

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0 || t > tf - aMax/jMax) {
                continue;
            }

            const double h1 = ((a0_a0 - af_af)/2 + jMax*(jMax*t*t + vd))/aMax;

            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (h1 - aMax)/jMax;
            profile.t[2] = aMax/jMax;
            profile.t[3] = tf - (h1 + ad + aMax)/jMax - 2*t;
            profile.t[4] = t;
            profile.t[5] = 0;
            profile.t[6] = af/jMax + t;
            
            if (profile.check<Teeth::UDDU, Limits::ACC0_VEL>(tf, pf, vf, af, jMax, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (-2*aMax)/jMax;
        polynom[2] = -(a0_a0 + af_af - 2*(a0 + af)*aMax + aMax_aMax + 2*jMax*(vd - aMax*tf))/jMax_jMax;
        polynom[3] = 0;
        polynom[4] = (3*(a0_p4 + af_p4) - 4*(af_p3 + 2*a0_p3)*aMax + 6*a0_a0*(af_af + aMax_aMax + 2*jMax*vd) - 12*a0*aMax*(af_af + 2*jMax*vd) + 6*af_af*(aMax_aMax - 2*aMax*jMax*tf + 2*jMax*vd) - ph1)/(12*jMax_jMax*jMax_jMax);

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0 || t > tf - aMax/jMax) {
                continue;
            }

            const double h1 = ((a0_a0 + af_af)/2 + jMax*(vd - jMax*t*t))/aMax;

            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (h1 - aMax)/jMax;
            profile.t[2] = aMax/jMax;
            profile.t[3] = tf - (h1 - a0 - af + aMax)/jMax - 2*t;
            profile.t[4] = t;
            profile.t[5] = 0;
            profile.t[6] = -(af/jMax) + t;
            
            if (profile.check<Teeth::UDUD, Limits::ACC0_VEL>(tf, pf, vf, af, jMax, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool Step2::time_up_vel(Profile& profile, double vMax, double aMax, double jMax) {
    const double g1 = (-pd + tf*v0);

    // Profile UDDU
    {
        const double p1 = af_af - 2*jMax*(-2*af*tf + jMax*tf_tf + 3*vd);
        const double ph1 = af_p3 - 3*jMax_jMax*g1 - 3*af*jMax*vd;
        const double ph2 = af_p4 + 8*af_p3*jMax*tf - 24*jMax_jMax*jMax*tf*g1 + 12*jMax*(-af_af*vd + 3*jMax*vd_vd + 2*af*jMax*(-pd + 2*tf*v0 - tf*vf));
        const double ph3 = a0*(af - jMax*tf);

        // Find root of 5th order polynom
        std::array<double, 6> polynom;
        polynom[0] = 1.0;
        polynom[1] = (15*a0_a0 + af_af + 4*af*jMax*tf - 16*ph3 - 2*jMax*(jMax*tf_tf + 3*vd))/(4*jMax*(-ad + jMax*tf));
        polynom[2] = (29*a0_p3 - 2*af_p3 - 33*a0*ph3 + 6*jMax_jMax*g1 + 6*af*jMax*vd + 6*a0*p1)/(6*jMax_jMax*(-ad + jMax*tf));
        polynom[3] = (61*a0_p4 + ph2 - 76*a0_a0*ph3 - 16*a0*ph1 + 30*a0_a0*p1)/(24*jMax_jMax*jMax*(-ad + jMax*tf));
        polynom[4] = (a0*(7*a0_p4 + ph2 - 10*a0_a0*ph3 - 4*a0*ph1 + 6*a0_a0*p1))/(12*jMax_jMax*jMax_jMax*(-ad + jMax*tf));
        polynom[5] = (7*a0_p6 + af_p6 - 12*a0_p4*ph3 + 48*af_p3*jMax_jMax*g1 - 8*a0_p3*ph1 - 72*jMax_jMax*jMax*(jMax*g1*g1 + vd_vd*vd) - 6*af_p4*jMax*vd - 144*af*jMax_jMax*jMax*g1*vd + 36*af_af*jMax_jMax*vd_vd + 9*a0_p4*p1 + 3*a0_a0*ph2)/(144*jMax_jMax*jMax_jMax*jMax*(-ad + jMax*tf));

        std::array<double, 5> deriv = Roots::polyMonicDeri(polynom);

        // Solve 4th order derivative analytically
        auto extremas = Roots::solveQuartMonic(deriv);
        
        std::set<double> roots;
        double tz_min {0.0};
        double tz_max = std::min<double>(tf, (tf - a0/jMax) / 2);
        double tz_current {tz_min};

        for (double tz: extremas) {
            if (tz <= 0.0 || tz >= tz_max) {
                continue;
            }

            const double res = std::abs(Roots::polyEval(Roots::polyDeri(deriv), tz)) * Roots::tolerance;
            const double val_new = Roots::polyEval(polynom, tz);
            if (std::abs(val_new) < res) {
                roots.insert(tz);
                
            } else if (Roots::polyEval(polynom, tz_current) * val_new < 0) {
                roots.insert(Roots::shrinkInterval(polynom, tz_current, tz));
            }
            tz_current = tz;
        }
        if (Roots::polyEval(polynom, tz_current) * Roots::polyEval(polynom, tz_max) < 0) {
            roots.insert(Roots::shrinkInterval(polynom, tz_current, tz_max));
        }

        for (double t: roots) {
            // Single Newton step (regarding pd)
            {
                const double h2 = Sqrt((a0_a0 + af_af)/2 + jMax*(2*a0*t + jMax*t*t - vd))/Abs(jMax);
                const double orig = -pd - (2*a0_p3 + 4*af_p3 + 24*a0*jMax*t*(af + jMax*(t - tf) + jMax*h2) + 6*a0_a0*(af + jMax*(2*t - tf) + jMax*h2) + 6*af_af*jMax*h2 + 12*af*jMax*(jMax*t*t - vd) + 12*jMax_jMax*(jMax*t*t*(t - tf) - tf*v0 - h2*(vd - jMax*t*t)))/(12*jMax_jMax);                
                const double deriv = -(a0 + jMax*t)*(3*h2 + (a0 + 2*af)/jMax + (3*t - 2*tf));
                t -= orig / deriv;
            }

            if (t < 0.0 || t > tf) {
                continue;
            }

            const double h1 = Sqrt((a0_a0 + af_af)/2 + jMax*(2*a0*t + jMax*t*t - vd))/Abs(jMax);
            
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = tf - 2*(t + h1) - (a0 + af)/jMax;
            profile.t[4] = h1;
            profile.t[5] = 0;
            profile.t[6] = profile.t[4] + af/jMax;

            if (profile.check<Teeth::UDDU, Limits::VEL>(tf, pf, vf, af, jMax, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        const double ph1 = af_af - 2*jMax*(2*af*tf + jMax*tf_tf - 3*vd);
        const double ph2 = af_p3 - 3*jMax_jMax*g1 + 3*af*jMax*vd;
        const double ph3 = 2*jMax*tf*g1 + 3*vd_vd;
        const double ph4 = af_p4 - 8*af_p3*jMax*tf + 12*jMax*(jMax*ph3 + af_af*vd + 2*af*jMax*(g1 - tf*vd));
        const double ph5 = af + jMax*tf;

        // Find root of 6th order polynom
        std::array<double, 7> polynom;
        polynom[0] = 1.0;
        polynom[1] = (5*a0 - ph5)/jMax;
        polynom[2] = (39*a0_a0 - ph1 - 16*a0*ph5)/(4*jMax_jMax);
        polynom[3] = (55*a0_p3 - 33*a0_a0*ph5 - 6*a0*ph1 + 2*ph2)/(6*jMax_jMax*jMax);
        polynom[4] = (101*a0_p4 + ph4 - 76*a0_p3*ph5 - 30*a0_a0*ph1 + 16*a0*ph2)/(24*jMax_jMax*jMax_jMax);
        polynom[5] = (a0*(11*a0_p4 + ph4 - 10*a0_p3*ph5 - 6*a0_a0*ph1 + 4*a0*ph2))/(12*jMax_jMax*jMax_jMax*jMax);
        polynom[6] = (11*a0_p6 - af_p6 - 12*a0_p5*ph5 - 48*af_p3*jMax_jMax*g1 - 9*a0_p4*ph1 + 72*jMax_jMax*jMax*(jMax*g1*g1 - vd_vd*vd) - 6*af_p4*jMax*vd - 144*af*jMax_jMax*jMax*g1*vd - 36*af_af*jMax_jMax*vd_vd + 8*a0_p3*ph2 + 3*a0_a0*ph4)/(144*jMax_jMax*jMax_jMax*jMax_jMax);

        std::array<double, 6> deriv = Roots::polyMonicDeri(polynom);
        std::array<double, 5> dderiv = Roots::polyMonicDeri(deriv);

        auto dd_extremas = Roots::solveQuartMonic(dderiv);
        std::set<std::pair<double, double>> dd_tz_intervals;

        double tz_min {0.0};
        double tz_max = std::min<double>(tf, (tf - a0/jMax) / 2);
        double dd_tz_current {tz_min};

        for (double tz: dd_extremas) {
            if (tz <= 0.0 || tz >= tz_max) {
                continue;
            }

            if (std::abs(Roots::polyEval(dderiv, tz)) > Roots::tolerance) {
                tz -= Roots::polyEval(dderiv, tz) / Roots::polyEval(Roots::polyDeri(dderiv), tz);
            }

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
            const double tz = Roots::shrinkInterval(deriv, interval.first, interval.second);

            if (tz <= 0.0 || tz >= tz_max) {
                continue;
            }

            const double res = std::abs(Roots::polyEval(dderiv, tz)) * Roots::tolerance;
            const double p_val = Roots::polyEval(polynom, tz);

            if (std::abs(p_val) < res) {
                roots.insert(tz);
            } else if (Roots::polyEval(polynom, tz_current) * p_val < 0) {
                roots.insert(Roots::shrinkInterval(polynom, tz_current, tz));
            }
            tz_current = tz;
        }
        if (Roots::polyEval(polynom, tz_current) * Roots::polyEval(polynom, tz_max) < 0) {
            roots.insert(Roots::shrinkInterval(polynom, tz_current, tz_max));
        }

        for (double t: roots) {
            // Single Newton step (regarding pd)
            {
                const double h2 = (af_af - a0_a0)/2 - jMax*(2*a0*t + jMax*t*t - vd);
                const double orig = -pd + (af_p3 - a0_p3 - 12*a0*jMax_jMax*t*(t - tf) + 3*a0_a0*jMax*(-2*t + tf) - 6*af*h2)/(6*jMax_jMax) + (jMax*t*t*(-t + tf) + tf*v0) + std::pow(h2,1.5)/(jMax*Abs(jMax)); 
                const double deriv = -(a0 + jMax*t)*(3*jMax*Sqrt(h2)/Abs(jMax) - 2*af + (a0 + 3*jMax*t - 2*jMax*tf))/jMax;

                t -= orig / deriv;
            }

            const double h1 = Sqrt((af_af - a0_a0)/2 - jMax*(2*a0*t + jMax*t*t - vd))/Abs(jMax);
            
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = tf - 2*(t + h1) + ad/jMax;
            profile.t[4] = h1;
            profile.t[5] = 0;
            profile.t[6] = profile.t[4] - af/jMax;

            if (profile.check<Teeth::UDUD, Limits::VEL>(tf, pf, vf, af, jMax, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool Step2::time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    if (std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        const double h1 = (v0 + vf)/aMax - (vd_vd + 4*aMax*pd)/(2*aMax_aMax*tf);

        profile.t[0] = tf/2 + h1;
        profile.t[1] = -(tf + 4*h1 - vd/aMax)/2;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = -(tf + 4*h1 + vd/aMax)/2;
        profile.t[6] = profile.t[0];
        const double jf = aMax/profile.t[0];

        return profile.check<Teeth::UDDU, Limits::ACC0_ACC1>(tf, pf, vf, af, jf, vMax, aMax, jMax);
    }

    const double h1 = -2*ad_ad*(3*(a0_p3 - af_p3) - 4*(a0_a0 + af_af)*aMax + 3*af*a0_a0 + 12*(ad + 2*aMax)*aMax_aMax - a0*(3*af_af + 16*af*aMax));
    const double h2 = aMax_aMax*tf_tf - vd_vd - 2*aMax*(2*pd - tf*(v0 + vf));
    const double h3 = 2*aMax_aMax*aMax*tf + (af_af + 2*af*aMax)*(aMax*tf - vd) + (a0_a0 - 2*a0*aMax)*(aMax*tf + vd);
    const double h4 = Sqrt(h3*h3 - h1*h2/3);
    const double jf = (h3 + h4)/(2*h2);

    profile.t[0] = (-a0 + aMax)/jf;
    profile.t[1] = (a0_a0 - af_af + 2*(a0 - af)*aMax - 8*aMax_aMax + 6*aMax*tf + 6*vd)/(12*aMax);
    profile.t[2] = profile.t[0] + a0/jf;
    profile.t[3] = 0;
    profile.t[4] = profile.t[2];
    profile.t[5] = tf - (profile.t[0] + profile.t[1] + profile.t[3] + 3*profile.t[2] + af/jf);
    profile.t[6] = profile.t[4] + af/jf;
    
    return profile.check<Teeth::UDDU, Limits::ACC0_ACC1>(tf, pf, vf, af, jf, vMax, aMax, jMax);
}

bool Step2::time_up_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    // a3 != 0
    
    // Case UDDU, Solution 2
    {
        const double h0a = a0_p3 - af_p3 + 3*a0_a0*aMax + 3*a0*aMax_aMax + 3*aMax_aMax*jMax*tf - 3*af*aMax*(aMax - 2*jMax*tf) - 3*af_af*(aMax - jMax*tf) - 3*jMax_jMax*(-2*pd + aMax*tf_tf + 2*tf*vf);
        const double h0b = a0_a0 + af_af + 2*(a0 + af)*aMax + 2*(aMax_aMax - jMax*(aMax*tf + vd));
        const double h0c = a0_p4 + 3*af_p4 + 4*(a0_p3 + 2*af_p3)*aMax + 6*a0_a0*aMax_aMax + 6*af_af*(aMax_aMax - 2*jMax*vd) + 12*jMax*(-2*aMax*jMax*(-pd + tf*v0) - aMax_aMax*vd + jMax*vd_vd) - 24*af*aMax*jMax*vd - 4*a0*(af_p3 + 3*af*aMax*(aMax - 2*jMax*tf) + 3*af_af*(aMax - jMax*tf) + 3*jMax*(-(aMax_aMax*tf) + jMax*(-2*pd + aMax*tf_tf + 2*tf*vf)));
        const double h1 = Abs(jMax)/jMax*Sqrt(4*h0a*h0a - 6*h0b*h0c);
        const double h2 = 6*jMax*h0b;

        profile.t[0] = 0;
        profile.t[1] = 0;
        profile.t[2] = (2*h0a + h1)/h2;
        profile.t[3] = -(a0_a0 + af_af + 2*(a0 + af)*aMax + 2*(aMax_aMax - aMax*jMax*tf - jMax*vd))/(2*jMax*(a0 + aMax - jMax*profile.t[2]));
        profile.t[4] = (aMax + a0)/jMax - profile.t[2];
        profile.t[5] = tf - (profile.t[2] + profile.t[3] + profile.t[4] + (af + aMax)/jMax);
        profile.t[6] = (af + aMax)/jMax;

        if (profile.check<Teeth::UDDU, Limits::ACC1>(tf, pf, vf, af, jMax, vMax, aMax)) {
            return true;
        }
    }

    // Case UDUD, Solution 1
    {
        const double h0a = -a0_p3 + af_p3 + 3*a0_a0*aMax - 3*a0*aMax_aMax + 3*af*aMax*(aMax - 2*jMax*tf) - 3*af_af*(aMax - jMax*tf) + 3*jMax*(aMax_aMax*tf + jMax*(-2*pd - aMax*tf_tf + 2*tf*vf));
        const double h0b = a0_a0 - af_af + 2*ad*aMax + 2*jMax*(aMax*tf - vd);
        const double h0c = a0_p4 + 3*af_p4 - 4*(a0_p3 + 2*af_p3)*aMax + 6*a0_a0*aMax_aMax - 24*af*aMax*jMax*vd + 12*jMax*(2*aMax*jMax*(-pd + tf*v0) + jMax*vd_vd + aMax_aMax*vd) + 6*af_af*(aMax_aMax + 2*jMax*vd) - 4*a0*(af_p3 + 3*af*aMax*(aMax - 2*jMax*tf) - 3*af_af*(aMax - jMax*tf) + 3*jMax*(aMax_aMax*tf + jMax*(-2*pd - aMax*tf_tf + 2*tf*vf)));
        const double h1 = Abs(jMax)/jMax*Sqrt(4*h0a*h0a - 6*h0b*h0c);
        const double h2 = 6*jMax*h0b;

        profile.t[0] = 0;
        profile.t[1] = 0;
        profile.t[2] = -(2*h0a + h1)/h2;
        profile.t[3] = 2*h1/h2;
        profile.t[4] = (aMax - a0)/jMax + profile.t[2];
        profile.t[5] = tf - (profile.t[2] + profile.t[3] + profile.t[4] + (-af + aMax)/jMax);
        profile.t[6] = (-af + aMax)/jMax;

        if (profile.check<Teeth::UDUD, Limits::ACC1>(tf, pf, vf, af, jMax, vMax, aMax)) {
            return true;
        }
    }
    return false;
}

bool Step2::time_up_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    // a3 != 0

    const double h0a = a0_p3 + 2*af_p3 - 6*(af_af + aMax_aMax)*aMax - 6*(a0 + af)*aMax*jMax*tf + 9*aMax_aMax*(af + jMax*tf) + 3*a0*aMax*(-2*af + 3*aMax) + 3*a0_a0*(af - 2*aMax + jMax*tf) - 6*jMax_jMax*(-pd + tf*v0) + 6*(af - aMax)*jMax*vd - 3*aMax*jMax_jMax*tf_tf;
    const double h0b = a0_a0 + af_af + 2*(aMax_aMax - (a0 + af)*aMax + jMax*(vd - aMax*tf));
    const double h1 = Abs(jMax)/jMax*Sqrt(4*h0a*h0a - 18*h0b*h0b*h0b);
    const double h2 = 6*jMax*h0b;

    // Solution 1, UDDU?
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = ad/jMax - 2 * profile.t[0] - (2*h0a - h1)/h2 + tf;
    profile.t[2] = -(2*h0a + h1)/h2;
    profile.t[3] = (2*h0a - h1)/h2;
    profile.t[4] = tf - (profile.t[0] + profile.t[1] + profile.t[2] + profile.t[3]);
    profile.t[5] = 0;
    profile.t[6] = 0;

    if (profile.check<Teeth::UDDU, Limits::ACC0>(tf, pf, vf, af, jMax, vMax, aMax)) {
        return true;
    }

    return false;
}

bool Step2::time_up_none(Profile& profile, double vMax, double aMax, double jMax) {
    if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        const double h1 = Sqrt(tf_tf*vf_vf + Power(4*pd - tf*vf,2));
        const double jf = 4*(4*pd - 2*tf*vf + h1)/tf_p3;
        
        profile.t[0] = tf/4;
        profile.t[1] = 0;
        profile.t[2] = 2*profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        if (profile.check<Teeth::UDDU, Limits::NONE>(tf, pf, vf, af, jf, vMax, aMax, jMax)) {
            return true;
        }
    }

    if (std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        // Solution 1
        {
            const double h1 = Sqrt(tf_tf*vd_vd + 4*Power(2*pd - tf*(v0 + vf),2));
            const double h2 = Sqrt(16*pd*(pd - tf*(v0 + vf)) + tf_tf*(5*v0_v0 + 6*v0*vf + 5*vf_vf));
            const double jf = 4*(4*pd - 2*tf*(v0 + vf) - h2)/tf_p3;
        
            profile.t[0] = (tf*(v0 + 3*vf) - 4*pd)/(4*vd);
            profile.t[1] = 0;
            profile.t[2] = tf/2;
            profile.t[3] = 0;
            profile.t[4] = 0;
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];
            
            if (profile.check<Teeth::UDDU, Limits::NONE>(tf, pf, vf, af, jf, vMax, aMax, jMax)) {
                return true;
            }
        }
    }

    // Profiles with a3 != 0, Solution UDDU
    {
        // First acc, then constant
        {
            const double ph1 = af + jMax*tf;

            std::array<double, 5> polynom;
            polynom[0] = 1.0;
            polynom[1] = (-2*(ad + jMax*tf))/jMax;
            polynom[2] = (2*a0_a0 + 2*af_af + 2*af*jMax*tf - 4*a0*ph1 + jMax*(jMax*tf_tf + 2*vd))/jMax_jMax;
            polynom[3] = (2*(a0_p3 - af_p3 - 3*af_af*jMax*tf - 3*a0_a0*ph1 + 3*a0*ph1*ph1 - 6*jMax_jMax*(-pd + tf*vf)))/(3*jMax_jMax*jMax);
            polynom[4] = (a0_p4 + af_p4 + 4*af_p3*jMax*tf + 6*af_af*jMax_jMax*tf_tf - 4*a0_p3*ph1 + 6*a0_a0*ph1*ph1 + 24*af*jMax_jMax*(-pd + tf*v0) - 4*a0*(af_p3 + 3*af_af*jMax*tf + 6*jMax_jMax*(-pd + tf*vf)) + 12*jMax_jMax*(vd_vd + jMax*tf*(-2*pd + tf*(v0 + vf))))/(12*jMax_jMax*jMax_jMax);
            auto roots = Roots::solveQuartMonic(polynom);

            for (double t: roots) {
                if (t < 0.0 || t > tf) {
                    continue;
                }

                // Single Newton step (regarding pd)
                {
                    const double h1 = (a0_a0 + af_af - 2*af*jMax*t - 2*a0*(af + jMax*(-t + tf)) + 2*jMax*(jMax*t*(t - tf) + vd))/(2*jMax*(-ad + 2*jMax*t - jMax*tf));
                    const double h2 = (-a0_a0 - af_af + 2*jMax_jMax*t*t + af*jMax*tf - 2*jMax_jMax*t*tf + jMax_jMax*tf_tf + a0*(-ad + 2*af + jMax*tf) + ad*(af - 2*jMax*t + jMax*tf) - 2*jMax*vd)/Power(ad + jMax*(-2*t + tf),2);
                    const double orig = (-a0_p3 + af_p3 + 3*(a0_a0 + af_af)*jMax*(h1 - t) + 3*(af - a0)*jMax_jMax*Power(h1 - t,2) + 3*a0*af*(a0 - af) - 3*a0*(2*af*jMax*(h1 - t)) + 3*jMax_jMax*(a0*tf_tf - 2*pd + 2*tf*v0 + h1*h1*jMax*(tf - 2*t) + jMax*tf*(2*h1*t - t*t - (h1 - t)*tf)))/(6*jMax_jMax);
                    const double deriv = ((a0 - af + 2*jMax*t - jMax*tf)*(a0 - af + jMax*tf)*(h2 - 1))/(2*jMax) + h1*(a0 - af + jMax*tf + (-a0 + af - 2*jMax*t + jMax*tf)*h2 - jMax*h1);
                    
                    t -= orig / deriv;
                }

                profile.t[0] = t;
                profile.t[1] = 0;
                profile.t[2] = (a0_a0 + af_af - 2*af*jMax*t - 2*a0*(af + jMax*(-t + tf)) + 2*jMax*(jMax*t*(t - tf) + vd))/(2*jMax*(-ad + 2*jMax*t - jMax*tf));
                profile.t[3] = ad/jMax + (tf - 2*t);
                profile.t[4] = tf - (t + profile.t[2] + profile.t[3]);
                profile.t[5] = 0;
                profile.t[6] = 0;
                
                if (profile.check<Teeth::UDDU, Limits::NONE>(tf, pf, vf, af, jMax, vMax, aMax)) {
                    return true;
                }
            }
        }

        // First constant, then acc
        {
            const double ph1 = ad_ad + 2*(af + a0)*jMax*tf - jMax*(jMax*tf_tf + 4*vd);
            const double ph2 = jMax*tf_tf*(-pd + tf*v0) - vd*(-2*pd - tf*v0 + 3*tf*vf);
            const double ph3 = 5*af_af - 8*af*jMax*tf + 2*jMax*(2*jMax*tf_tf - vd);
            const double ph4 = jMax_jMax*tf_p4 - 2*vd_vd + 8*jMax*tf*(-pd + tf*vf);
            const double ph5 = (5*af_p4 - 8*af_p3*jMax*tf - 12*af_af*jMax*(jMax*tf_tf + vd) + 24*af*jMax_jMax*(-2*pd + jMax*tf_p3 + 2*tf*vf) - 6*jMax_jMax*ph4);

            std::array<double, 5> polynom;
            polynom[0] = 1.0;
            polynom[1] = -(4*a0_p3 - 4*af_p3 - 12*a0_a0*(af - jMax*tf) + 6*a0*(2*af_af - 2*af*jMax*tf + jMax*(jMax*tf_tf - 2*vd)) + 6*af*jMax*(3*jMax*tf_tf + 2*vd) - 6*jMax_jMax*(-4*pd + jMax*tf_p3 - 2*tf*v0 + 6*tf*vf))/(3*jMax*ph1);
            polynom[2] = -(-a0_p4 - af_p4 + 4*a0_p3*(af - jMax*tf) + a0_a0*(-6*af_af + 8*af*jMax*tf - 4*jMax*(jMax*tf_tf - vd)) + 2*af_af*jMax*(jMax*tf_tf + 2*vd) - 4*af*jMax_jMax*(-3*pd + jMax*tf_p3 + 2*tf*v0 + tf*vf) + jMax_jMax*(jMax_jMax*tf_p4 - 8*vd_vd + 4*jMax*tf*(-3*pd + tf*v0 + 2*tf*vf)) + 2*a0*(2*af_p3 - 2*af_af*jMax*tf + af*jMax*(-3*jMax*tf_tf + 4*v0 - 4*vf) + jMax_jMax*(-6*pd + jMax*tf_p3 - 4*tf*v0 + 10*tf*vf)))/(jMax_jMax*ph1);
            polynom[3] = -(a0_p5 - af_p5 + af_p4*jMax*tf - 5*a0_p4*(af - jMax*tf) + 2*a0_p3*ph3 + 4*af_p3*jMax*(jMax*tf_tf + vd) - 12*af_af*jMax_jMax*(-2*pd + tf*(v0 + vf)) + 12*af*jMax_jMax*(-vd_vd + jMax*tf*(-2*pd + 3*tf*v0 - tf*vf)) - 2*a0_a0*(5*af_p3 - 9*af_af*jMax*tf - 6*af*jMax*vd + 6*jMax_jMax*(-2*pd - tf*v0 + 3*tf*vf)) - 12*jMax_jMax*jMax*ph2 + a0*ph5)/(3*jMax_jMax*jMax*ph1);
            polynom[4] = -(-a0_p6 - af_p6 + 6*a0_p5*(af - jMax*tf) - 48*af_p3*jMax_jMax*(-pd + tf*v0) + 72*jMax_jMax*jMax*(jMax*Power(-pd + tf*v0,2) + vd_vd*vd) - 3*a0_p4*ph3 + 144*af*jMax_jMax*jMax*(-pd + tf*v0)*vd - 36*af_af*jMax_jMax*vd_vd + 6*af_p4*jMax*vd + 4*a0_p3*(5*af_p3 - 9*af_af*jMax*tf - 6*af*jMax*vd + 6*jMax_jMax*(-2*pd - tf*v0 + 3*tf*vf)) - 3*a0_a0*ph5 + 6*a0*(af_p5 - af_p4*jMax*tf - 4*af_p3*jMax*(jMax*tf_tf + vd) + 12*af_af*jMax_jMax*(-2*pd + tf*(v0 + vf)) - 12*af*jMax_jMax*(-vd_vd + jMax*tf*(-2*pd + 3*tf*v0 - tf*vf)) + 12*jMax_jMax*jMax*ph2))/(18*jMax_jMax*jMax_jMax*ph1);
            auto roots = Roots::solveQuartMonic(polynom);

            for (double t: roots) {
                if (t < 0.0 || t > tf) {
                    continue;
                }

                // Single Newton step (regarding pd)
                {
                    const double h1 = a0_a0 + af_af + 2*af*jMax*t - 2*a0*(af + jMax*(t - tf)) + 2*jMax*(jMax*t*(t - tf) - vd);
                    const double h2 = (a0 - af + jMax*(-2*t + tf));
                    const double orig = (af_p3 - a0_p3 + 3*af*jMax*t*(af + jMax*t) + 3*a0_a0*(af + jMax*t) - 3*a0*(af_af + 2*af*jMax*t + jMax_jMax*(t*t - tf_tf)) + 3*jMax_jMax*(-2*pd + jMax*t*(t - tf)*tf + 2*tf*v0))/(6*jMax_jMax) - (jMax*std::pow(h1,1.5))/(2*Sqrt(2)*jMax_jMax*Abs(jMax)) + ((a0 - af - jMax*t)*h1)/(2*jMax_jMax);
                    const double deriv = (3*Sqrt(2)*jMax*h2*Sqrt(h1)/Abs(jMax) - 2*(3*a0_a0 - 6*a0*af + 3*af_af + af*jMax*(8*t - 2*tf) + 4*a0*jMax*(-2*t + tf) + 2*jMax*(jMax*t*(3*t - 2*tf) - vd)) + 2*(a0 - af - jMax*tf)*h2)/(4*jMax);

                    t -= orig / deriv;
                }

                const double h1 = Sqrt(2*(a0_a0 + af_af + 2*af*jMax*t - 2*a0*(af + jMax*(t - tf)) + 2*jMax*(jMax*t*(t - tf) - vd)))/Abs(jMax);

                // Solution 2 with aPlat
                profile.t[0] = 0;
                profile.t[1] = 0;
                profile.t[2] = t;
                profile.t[3] = tf - 2*t - ad/jMax - h1;
                profile.t[4] = h1/2;
                profile.t[5] = 0;
                profile.t[6] = tf - (t + profile.t[3] + profile.t[4]);

                if (profile.check<Teeth::UDDU, Limits::NONE>(tf, pf, vf, af, jMax, vMax, aMax)) {
                    return true;
                }
            }       
        }
    }

    // Profiles with a3 != 0, Solution UDUD
    {
        // First constant, then acc
        {
            const double ph1 = -ad + jMax*tf;
            const double ph2 = jMax*tf_tf*(-pd + tf*v0) - vd*(-2*pd - tf*v0 + 3*tf*vf);
            const double ph3 = 5*af_af - 8*af*jMax*tf + 2*jMax*(2*jMax*tf_tf - vd);
            const double ph4 = jMax_jMax*tf_p4 - 2*vd_vd + 8*jMax*tf*(-pd + tf*vf);
            const double ph5 = (5*af_p4 - 8*af_p3*jMax*tf - 12*af_af*jMax*(jMax*tf_tf + vd) + 24*af*jMax_jMax*(-2*pd + jMax*tf_p3 + 2*tf*vf) - 6*jMax_jMax*ph4);
            const double ph6 = -vd_vd + jMax*tf*(-2*pd + 3*tf*v0 - tf*vf);
            const double ph7 = jMax_jMax*ph1*ph1;

            std::array<double, 5> polynom;
            polynom[0] = 1.0;
            polynom[1] = (4*af*tf - 2*jMax*tf_tf - 4*vd)/ph1;
            polynom[2] = (-2*(a0_p4 + af_p4) + 8*af_p3*jMax*tf + 6*af_af*jMax_jMax*tf_tf + 8*a0_p3*(af - jMax*tf) - 12*a0_a0*Power(af - jMax*tf,2) - 12*af*jMax_jMax*(-pd + jMax*tf_p3 - 2*tf*v0 + 3*tf*vf) + 2*a0*(4*af_p3 - 12*af_af*jMax*tf + 9*af*jMax_jMax*tf_tf - 3*jMax_jMax*(2*pd + jMax*tf_p3 - 2*tf*vf)) + 3*jMax_jMax*(jMax_jMax*tf_p4 + 4*vd_vd - 4*jMax*tf*(pd + tf*v0 - 2*tf*vf)))/(3*ph7);
            polynom[3] = (-a0_p5 + af_p5 - af_p4*jMax*tf + 5*a0_p4*(af - jMax*tf) - 2*a0_p3*ph3 - 4*af_p3*jMax*(jMax*tf_tf + vd) + 12*af_af*jMax_jMax*(-2*pd + tf*(v0 + vf)) - 12*af*jMax_jMax*ph6 + 2*a0_a0*(5*af_p3 - 9*af_af*jMax*tf - 6*af*jMax*vd + 6*jMax_jMax*(-2*pd - tf*v0 + 3*tf*vf)) + 12*jMax_jMax*jMax*ph2 + a0*(-5*af_p4 + 8*af_p3*jMax*tf + 12*af_af*jMax*(jMax*tf_tf + vd) - 24*af*jMax_jMax*(-2*pd + jMax*tf_p3 + 2*tf*vf) + 6*jMax_jMax*ph4))/(3*jMax*ph7);
            polynom[4] = -(a0_p6 + af_p6 - 6*a0_p5*(af - jMax*tf) + 48*af_p3*jMax_jMax*(-pd + tf*v0) - 72*jMax_jMax*jMax*(jMax*Power(-pd + tf*v0,2) + vd_vd*vd) + 3*a0_p4*ph3 - 6*af_p4*jMax*vd - 144*af*jMax_jMax*jMax*(-pd + tf*v0)*vd + 36*af_af*jMax_jMax*vd_vd - 4*a0_p3*(5*af_p3 - 9*af_af*jMax*tf - 6*af*jMax*vd + 6*jMax_jMax*(-2*pd - tf*v0 + 3*tf*vf)) + 3*a0_a0*ph5 - 6*a0*(af_p5 - af_p4*jMax*tf - 4*af_p3*jMax*(jMax*tf_tf + vd) + 12*af_af*jMax_jMax*(-2*pd + tf*(v0 + vf)) - 12*af*jMax_jMax*ph6 + 12*jMax_jMax*jMax*ph2))/(18*jMax_jMax*ph7);
            auto roots = Roots::solveQuartMonic(polynom);

            for (double t: roots) {
                if (t < 0.0 || t > tf) {
                    continue;
                }

                double h1 = Sqrt(((a0_a0 + af_af)/2 - af*(a0 + jMax*t) + a0*jMax*(t + tf) + jMax*(jMax*t*tf - vd)))/Abs(jMax);

                profile.t[0] = t;
                profile.t[1] = tf - ad/jMax - 2*h1;
                profile.t[2] = h1;
                profile.t[3] = 0;
                profile.t[4] = ad/jMax + h1 - t;
                profile.t[5] = 0;
                profile.t[6] = 0;
 
                if (profile.check<Teeth::UDUD, Limits::NONE>(tf, pf, vf, af, jMax, vMax, aMax)) {
                    return true;
                }
            }
        }
    }

    return false;
}

bool Step2::time_down_acc0_acc1_vel(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, vMin, -aMax, -jMax);
}

bool Step2::time_down_acc1_vel(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc1_vel(profile, vMin, -aMax, -jMax);
}

bool Step2::time_down_acc0_vel(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc0_vel(profile, vMin, -aMax, -jMax);
}

bool Step2::time_down_vel(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_vel(profile, vMin, -aMax, -jMax);
}

bool Step2::time_down_acc0_acc1(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, vMin, -aMax, -jMax);
}

bool Step2::time_down_acc1(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc1(profile, vMin, -aMax, -jMax);
}

bool Step2::time_down_acc0(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_acc0(profile, vMin, -aMax, -jMax);
}

bool Step2::time_down_none(Profile& profile, double vMin, double aMax, double jMax) {
    return time_up_none(profile, vMin, -aMax, -jMax);
}

bool Step2::get_profile(Profile& profile) {
    profile.a[0] = a0;
    profile.v[0] = v0;
    profile.p[0] = p0;

    // Test all cases to get ones that match
    if (pf > p0) {
        return time_up_acc0_acc1_vel(profile, vMax, aMax, jMax)
            || time_down_acc0_acc1_vel(profile, vMin, aMax, jMax)
            || time_up_acc0_vel(profile, vMax, aMax, jMax)
            || time_down_acc0_vel(profile, vMin, aMax, jMax)
            || time_up_acc1_vel(profile, vMax, aMax, jMax)
            || time_down_acc1_vel(profile, vMin, aMax, jMax)
            || time_up_vel(profile, vMax, aMax, jMax)
            || time_down_vel(profile, vMin, aMax, jMax)
            || time_up_none(profile, vMax, aMax, jMax)
            || time_up_acc0(profile, vMax, aMax, jMax)
            || time_up_acc1(profile, vMax, aMax, jMax)
            || time_up_acc0_acc1(profile, vMax, aMax, jMax)
            || time_down_acc0(profile, vMin, aMax, jMax)
            || time_down_acc1(profile, vMin, aMax, jMax)
            || time_down_acc0_acc1(profile, vMin, aMax, jMax)
            || time_down_none(profile, vMin, aMax, jMax);

    } else {
        return time_down_acc0_acc1_vel(profile, vMin, aMax, jMax)
            || time_up_acc0_acc1_vel(profile, vMax, aMax, jMax)
            || time_down_acc0_vel(profile, vMin, aMax, jMax)
            || time_up_acc0_vel(profile, vMax, aMax, jMax)
            || time_down_acc1_vel(profile, vMin, aMax, jMax)
            || time_up_acc1_vel(profile, vMax, aMax, jMax)
            || time_down_vel(profile, vMin, aMax, jMax)
            || time_up_vel(profile, vMax, aMax, jMax)
            || time_down_none(profile, vMin, aMax, jMax)
            || time_down_acc0(profile, vMin, aMax, jMax)
            || time_down_acc1(profile, vMin, aMax, jMax)
            || time_down_acc0_acc1(profile, vMin, aMax, jMax)
            || time_up_acc0(profile, vMax, aMax, jMax)
            || time_up_acc1(profile, vMax, aMax, jMax)
            || time_up_acc0_acc1(profile, vMax, aMax, jMax)
            || time_up_none(profile, vMax, aMax, jMax);
    }
}

} // namespace ruckig
