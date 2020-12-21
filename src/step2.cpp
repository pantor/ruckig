#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>
#include <ruckig/wolfram.hpp>


namespace ruckig {

RuckigStep2::RuckigStep2(double tf, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax): tf(tf), p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af) {
    
}

bool RuckigStep2::time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Profile UDDU
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = -(-3*Power(a0,2)*aMax*jMax + 3*Power(af,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 6*af*Power(aMax,2)*jMax + 18*Power(aMax,3)*jMax - 6*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,2)*v0 - 6*aMax*Power(jMax,2)*vf + Sqrt(3)*Sqrt(Power(aMax,2)*Power(jMax,2)*(-3*Power(a0,4) - 3*Power(af,4) + 4*Power(a0,3)*aMax - 4*Power(af,3)*aMax + 6*Power(a0,2)*(Power(af,2) + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))) - 12*a0*aMax*(Power(af,2) + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))) + 12*Power(af,2)*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)) + 24*af*aMax*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)) + 12*(Power(aMax,4) - 2*Power(aMax,3)*jMax*tf + Power(aMax,2)*Power(jMax,2)*Power(tf,2) - Power(jMax,2)*Power(v0 - vf,2) + 2*aMax*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf))))))/(12.*Power(aMax,2)*Power(jMax,2));
        profile.t[2] = aMax/jMax;
        profile.t[3] = (-6*Power(aMax,3)*jMax + Sqrt(3)*Sqrt(Power(aMax,2)*Power(jMax,2)*(-3*Power(a0,4) - 3*Power(af,4) + 4*Power(a0,3)*aMax - 4*Power(af,3)*aMax + 6*Power(a0,2)*(Power(af,2) + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))) - 12*a0*aMax*(Power(af,2) + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))) + 12*Power(af,2)*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)) + 24*af*aMax*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)) + 12*(Power(aMax,4) - 2*Power(aMax,3)*jMax*tf + Power(aMax,2)*Power(jMax,2)*Power(tf,2) - Power(jMax,2)*Power(v0 - vf,2) + 2*aMax*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf))))))/(6.*Power(aMax,2)*Power(jMax,2));
        profile.t[4] = profile.t[2];
        profile.t[5] = -(3*Power(a0,2)*aMax*jMax - 3*Power(af,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 6*af*Power(aMax,2)*jMax + 18*Power(aMax,3)*jMax - 6*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,2)*v0 + 6*aMax*Power(jMax,2)*vf + Sqrt(3)*Sqrt(Power(aMax,2)*Power(jMax,2)*(-3*Power(a0,4) - 3*Power(af,4) + 4*Power(a0,3)*aMax - 4*Power(af,3)*aMax + 6*Power(a0,2)*(Power(af,2) + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))) - 12*a0*aMax*(Power(af,2) + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))) + 12*Power(af,2)*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)) + 24*af*aMax*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)) + 12*(Power(aMax,4) - 2*Power(aMax,3)*jMax*tf + Power(aMax,2)*Power(jMax,2)*Power(tf,2) - Power(jMax,2)*Power(v0 - vf,2) + 2*aMax*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf))))))/(12.*Power(aMax,2)*Power(jMax,2));
        profile.t[6] = (af + aMax)/jMax;

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(tf, pf, vf, af, vMax, aMax)) {
            return true;
        }
    }

    // Profile UDUD
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (3*Power(a0,4) + 3*Power(af,4) - 4*Power(a0,3)*aMax - 8*Power(af,3)*aMax + 24*a0*Power(aMax,3) + 24*af*aMax*(Power(aMax,2) + jMax*(v0 - vf)) - 6*Power(af,2)*(Power(aMax,2) + 2*jMax*(v0 - vf)) + 6*Power(a0,2)*(Power(af,2) - 2*af*aMax - Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf) - 12*(2*Power(aMax,4) - 2*Power(aMax,3)*jMax*tf - 2*aMax*Power(jMax,2)*(p0 - pf + tf*v0) - Power(jMax,2)*Power(v0 - vf,2) + Power(aMax,2)*jMax*(-v0 + vf)))/(12.*aMax*jMax*(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf - jMax*v0 + jMax*vf)));
        profile.t[2] = aMax/jMax;
        profile.t[3] = -(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + 4*Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf)/(2.*aMax*jMax);
        profile.t[4] = profile.t[2];
        profile.t[5] = (3*Power(a0,4) + 3*Power(af,4) - 8*Power(a0,3)*aMax - 4*Power(af,3)*aMax + 24*af*Power(aMax,3) - 6*Power(af,2)*(Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*(v0 - vf)) + 6*Power(a0,2)*(Power(af,2) - Power(aMax,2) - 2*jMax*v0 + 2*jMax*vf) - 12*a0*aMax*(Power(af,2) - 2*(Power(aMax,2) + jMax*v0 - jMax*vf)) - 12*(2*Power(aMax,4) - 2*Power(aMax,3)*jMax*tf - Power(jMax,2)*Power(v0 - vf,2) + Power(aMax,2)*jMax*(-v0 + vf) + 2*aMax*Power(jMax,2)*(p0 - pf + tf*vf)))/(12.*aMax*jMax*(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf - jMax*v0 + jMax*vf)));
        profile.t[6] = (-af + aMax)/jMax;

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
        if (profile.check(tf, pf, vf, af, vMax, aMax)) {
            return true;
        }
    }
    
    return false;
}

bool RuckigStep2::time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Profile UDDU
    {
        std::array<double, 5> polynom;

        polynom[0] = 1.0;
        polynom[1] = (2*(2*a0 + aMax))/jMax;
        polynom[2] = (5*Power(a0,2) + Power(af,2) + 4*a0*aMax + 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf)/Power(jMax,2);
        polynom[3] = (2*a0*(Power(a0,2) + Power(af,2) + a0*aMax + 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf))/Power(jMax,3);
        polynom[4] = (3*Power(a0,4) + 3*Power(af,4) + 4*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 6*Power(af,2)*(Power(aMax,2) + 2*jMax*(v0 - vf)) + 12*jMax*(-2*aMax*jMax*(p0 - pf + tf*v0) + Power(aMax,2)*(v0 - vf) + jMax*Power(v0 - vf,2)) + 24*af*aMax*jMax*(v0 - vf) + 6*Power(a0,2)*(Power(af,2) + 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf))/(12.*Power(jMax,4));

        auto roots = Roots::solveQuart(polynom);
        for (double t: roots) {
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + t;
            profile.t[3] = -(Power(a0,2) + Power(af,2) + 2*a0*(aMax + 2*jMax*t) + 2*(af*aMax + Power(aMax,2) + aMax*jMax*(2*t - tf) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[4] = aMax/jMax;
            profile.t[5] = (Power(a0,2) + Power(af,2) - 2*Power(aMax,2) + 4*a0*jMax*t + 2*jMax*(jMax*Power(t,2) + v0 - vf))/(2.*aMax*jMax);
            profile.t[6] = (af + aMax)/jMax;
            
            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (4*a0 - 2*aMax)/jMax;
        polynom[2] = -((-5*Power(a0,2) + Power(af,2) + 4*a0*aMax - 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf)/Power(jMax,2));
        polynom[3] = (2*a0*(Power(a0,2) - Power(af,2) - a0*aMax + 2*af*aMax - Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf))/Power(jMax,3);
        polynom[4] = (3*Power(a0,4) + 3*Power(af,4) - 4*Power(a0,3)*aMax - 8*Power(af,3)*aMax + 24*af*aMax*jMax*(v0 - vf) - 6*Power(a0,2)*(Power(af,2) - 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf) + 12*jMax*(2*aMax*jMax*(p0 - pf + tf*v0) + jMax*Power(v0 - vf,2) + Power(aMax,2)*(-v0 + vf)) + 6*Power(af,2)*(Power(aMax,2) + 2*jMax*(-v0 + vf)))/(12.*Power(jMax,4));

        auto roots = Roots::solveQuart(polynom);
        for (double t: roots) {
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + t;
            profile.t[3] = (Power(a0,2) - Power(af,2) - 2*a0*aMax + 2*af*aMax - 2*Power(aMax,2) + 4*a0*jMax*t - 4*aMax*jMax*t + 2*Power(jMax,2)*Power(t,2) + 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf)/(2.*aMax*jMax);
            profile.t[4] = aMax/jMax;
            profile.t[5] = -(Power(a0,2) - Power(af,2) + 4*a0*jMax*t + 2*(Power(aMax,2) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[6] = (-af + aMax)/jMax;
            
            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool RuckigStep2::time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {   
    // Profile UDDU
    {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (2*aMax)/jMax;
        polynom[2] = (Power(a0,2) - Power(af,2) - 2*a0*aMax + 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf)/Power(jMax,2);
        polynom[3] = 0;
        polynom[4] = -(-3*Power(a0,4) - 3*Power(af,4) + 8*Power(a0,3)*aMax + 4*Power(af,3)*aMax - 12*a0*aMax*(Power(af,2) + 2*jMax*(v0 - vf)) + 6*Power(a0,2)*(Power(af,2) - Power(aMax,2) + 2*jMax*v0 - 2*jMax*vf) + 6*Power(af,2)*(Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(-v0 + vf)) + 12*jMax*(Power(aMax,2)*(v0 - vf) - jMax*Power(v0 - vf,2) + 2*aMax*jMax*(p0 - pf + tf*vf)))/(12.*Power(jMax,4));

        auto roots = Roots::solveQuart(polynom);
        for (double t: roots) {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (Power(a0,2) - Power(af,2) - 2*Power(aMax,2) + 2*jMax*(jMax*Power(t,2) - v0 + vf))/(2.*aMax*jMax);
            profile.t[2] = aMax/jMax;
            profile.t[3] = (-Power(a0,2) + Power(af,2) + 2*a0*aMax - 2*(af*aMax + Power(aMax,2) + aMax*jMax*(2*t - tf) + jMax*(jMax*Power(t,2) - v0 + vf)))/(2.*aMax*jMax);
            profile.t[4] = t;
            profile.t[5] = 0;
            profile.t[6] = af/jMax + t;
            
            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (-2*aMax)/jMax;
        polynom[2] = -((Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf)/Power(jMax,2));
        polynom[3] = 0;
        polynom[4] = (3*Power(a0,4) + 3*Power(af,4) - 8*Power(a0,3)*aMax - 4*Power(af,3)*aMax + 6*Power(a0,2)*(Power(af,2) + Power(aMax,2) - 2*jMax*v0 + 2*jMax*vf) - 12*a0*aMax*(Power(af,2) + 2*jMax*(-v0 + vf)) + 6*Power(af,2)*(Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(-v0 + vf)) - 12*jMax*(Power(aMax,2)*(v0 - vf) - jMax*Power(v0 - vf,2) + 2*aMax*jMax*(p0 - pf + tf*vf)))/(12.*Power(jMax,4));

        auto roots = Roots::solveQuart(polynom);
        for (double t: roots) {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (Power(a0,2) + Power(af,2) - 2*(Power(aMax,2) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[2] = aMax/jMax;
            profile.t[3] = -(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*(af*aMax - Power(aMax,2) + aMax*jMax*(-2*t + tf) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[4] = t;
            profile.t[5] = 0;
            profile.t[6] = -(af/jMax) + t;
            
            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool RuckigStep2::time_up_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Profile UDDU
    {
        // Find root of 5th order polynom
        std::array<double, 6> polynom;
        polynom[0] = 1.0;
        polynom[1] = (15*Power(a0,2) + 16*a0*jMax*tf - 2*Power(jMax,2)*Power(tf,2) + 6*jMax*v0 - 6*jMax*vf)/(4*a0*jMax + 4*Power(jMax,2)*tf);
        polynom[2] = (29*Power(a0,3) + 33*Power(a0,2)*jMax*tf + 6*Power(jMax,2)*(p0 - pf + tf*v0) - 12*a0*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf))/(6.*Power(jMax,2)*(a0 + jMax*tf));
        polynom[3] = (61*Power(a0,4) + 76*Power(a0,3)*jMax*tf + 48*a0*Power(jMax,2)*(p0 - pf + tf*v0) - 24*Power(jMax,3)*tf*(p0 - pf + tf*v0) + 36*Power(jMax,2)*Power(v0 - vf,2) - 60*Power(a0,2)*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf))/(24.*Power(jMax,3)*(a0 + jMax*tf));
        polynom[4] = (a0*(7*Power(a0,4) + 10*Power(a0,3)*jMax*tf + 12*a0*Power(jMax,2)*(p0 - pf + tf*v0) - 24*Power(jMax,3)*tf*(p0 - pf + tf*v0) + 36*Power(jMax,2)*Power(v0 - vf,2) - 12*Power(a0,2)*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf)))/(12.*Power(jMax,4)*(a0 + jMax*tf));
        polynom[5] = (7*Power(a0,6) + 12*Power(a0,5)*jMax*tf + 24*Power(a0,3)*Power(jMax,2)*(p0 - pf + tf*v0) - 36*Power(a0,2)*Power(jMax,2)*(2*jMax*tf*(p0 - pf + tf*v0) - 3*Power(v0 - vf,2)) - 72*Power(jMax,3)*(jMax*Power(p0 - pf + tf*v0,2) - Power(v0 - vf,3)) - 18*Power(a0,4)*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf))/(144.*Power(jMax,5)*(a0 + jMax*tf));

        std::array<double, 5> deriv;
        deriv[0] = 5 * polynom[0];
        deriv[1] = 4 * polynom[1];
        deriv[2] = 3 * polynom[2];
        deriv[3] = 2 * polynom[3];
        deriv[4] = polynom[4];

        // Solve 4th order derivative analytically
        auto extremas = Roots::solveQuart(deriv[0], deriv[1], deriv[2], deriv[3], deriv[4]);
        std::set<std::tuple<double, double>> tz_intervals;

        double tz_min {0.0};
        double tz_max = std::min<double>(tf, (tf - a0/jMax) / 2);
        double tz_current {tz_min};

        for (double tz: extremas) {
            if (tz <= 0.0 || tz >= tz_max) {
                continue;
            }

            // Check that polynom(lower) and polynom(upper) have different signs (should only happen at first and last boundary)
            double val_current = Roots::polyEval(polynom, tz_current);
            double val_new = Roots::polyEval(polynom, tz);
            if (std::abs(val_new) < 1e-15) {
                // if (val_current * val_new < 0) {
                //     tz_intervals.insert({tz_current, tz});
                //     tz += 1e-14;
                // } else {
                    tz_intervals.insert({tz - 1e-12, tz + 1e-12});
                    tz += 1e-14;
                // }
                
            } else if (val_current * val_new < 0) {
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
            double t = Roots::shrinkInterval(polynom, lower, upper, 2e-16);

            double vPlat = Power(a0,2)/(2.*jMax) + 2*a0*t + jMax*Power(t,2) + v0;
            double h1 = 2*Sqrt((vPlat - vf)/jMax);

            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + t;
            // profile.t[3] = -((23*Power(a0,7) + Power(a0,6)*jMax*(200*t + 39*tf) + 6*Power(a0,5)*jMax*(2*jMax*(45*Power(t,2) + 22*t*tf - 5*Power(tf,2)) + 13*(v0 - vf)) + 6*Power(a0,4)*Power(jMax,2)*(20*p0 - 20*pf + 2*jMax*t*(41*Power(t,2) + 52*t*tf - 23*Power(tf,2)) + 24*t*(v0 - vf) + 5*tf*(-3*v0 + 7*vf)) + 12*Power(a0,3)*Power(jMax,2)*(2*Power(jMax,2)*Power(t,2)*(6*Power(t,2) + 22*t*tf - 11*Power(tf,2)) - 19*Power(v0 - vf,2) + 2*jMax*(20*p0*t - 20*pf*t - 6*p0*tf + 6*pf*tf - 25*Power(t,2)*v0 - 26*t*tf*v0 + Power(tf,2)*v0 + 25*Power(t,2)*vf + 46*t*tf*vf - 7*Power(tf,2)*vf)) + 36*Power(a0,2)*Power(jMax,3)*(2*Power(jMax,2)*Power(t,3)*(2*t - tf)*tf - 2*jMax*(pf*(9*Power(t,2) + 2*t*tf - Power(tf,2)) + p0*(-9*Power(t,2) - 2*t*tf + Power(tf,2)) + 4*t*(3*t - tf)*(t + 2*tf)*v0 + (-12*Power(t,3) - 29*Power(t,2)*tf + 6*t*Power(tf,2) + Power(tf,3))*vf) - (v0 - vf)*(4*p0 - 4*pf + 40*t*(v0 - vf) + tf*(v0 + 3*vf))) + 72*a0*Power(jMax,3)*(-3*Power(v0 - vf,3) - 4*Power(jMax,2)*t*(pf*(Power(t,2) + 2*t*tf - Power(tf,2)) + p0*(-Power(t,2) - 2*t*tf + Power(tf,2)) + t*(Power(t,2) + 4*t*tf - 2*Power(tf,2))*v0 + (-Power(t,3) - 5*Power(t,2)*tf + Power(tf,3))*vf) + jMax*(Power(p0,2) + Power(pf,2) + 4*(-5*Power(t,2) + Power(tf,2))*Power(v0,2) + 2*(20*Power(t,2) - Power(tf,2))*v0*vf - (20*Power(t,2) + Power(tf,2))*Power(vf,2) + pf*tf*(-6*v0 + 4*vf) - 2*p0*(pf - 3*tf*v0 + 2*tf*vf))) + 72*Power(jMax,4)*(2*Power(jMax,2)*Power(t,2)*(2*t - tf)*tf*(p0 - pf - t*v0 + (t + tf)*vf) + Power(v0 - vf,2)*(4*p0 - 4*pf + tf*(v0 + 3*vf)) + jMax*(Power(p0,2)*(4*t - 3*tf) + Power(pf,2)*(4*t - 3*tf) - 2*pf*(Power(t,2) + 2*t*tf - Power(tf,2))*v0 + 2*t*(-3*Power(t,2) + Power(tf,2))*Power(v0,2) + 2*pf*(Power(t,2) - 2*t*tf + 2*Power(tf,2))*vf + 2*(6*Power(t,3) + Power(t,2)*tf - Power(tf,3))*v0*vf - (6*Power(t,3) + 2*Power(t,2)*tf - 2*t*Power(tf,2) + Power(tf,3))*Power(vf,2) - 2*p0*(pf*(4*t - 3*tf) + Power(t,2)*(-v0 + vf) - 2*t*tf*(v0 + vf) + Power(tf,2)*(v0 + 2*vf)))))/(jMax*(-Power(a0,6) + 6*Power(a0,4)*jMax*(v0 - vf) - 36*Power(a0,2)*Power(jMax,2)*Power(v0 - vf,2) + 48*Power(a0,3)*Power(jMax,2)*(p0 - pf + tf*vf) - 144*a0*Power(jMax,3)*(v0 - vf)*(p0 - pf + tf*vf) + 72*Power(jMax,3)*(Power(v0 - vf,3) + jMax*Power(p0 - pf + tf*vf,2)))));
            // profile.t[4] = (12*Power(a0,7) + Power(a0,6)*jMax*(101*t + 19*tf) + 6*Power(a0,5)*jMax*(45*jMax*Power(t,2) + 22*jMax*t*tf - 5*jMax*Power(tf,2) + 6*v0 - 6*vf) + 6*Power(a0,4)*Power(jMax,2)*(6*p0 - 6*pf + jMax*t*(41*Power(t,2) + 52*t*tf - 23*Power(tf,2)) - 7*tf*v0 + 11*t*(v0 - vf) + 13*tf*vf) + 12*Power(a0,3)*Power(jMax,2)*(Power(jMax,2)*Power(t,2)*(6*Power(t,2) + 22*t*tf - 11*Power(tf,2)) - 8*Power(v0 - vf,2) + jMax*(16*p0*t - 16*pf*t - 4*p0*tf + 4*pf*tf - 25*Power(t,2)*v0 - 26*t*tf*v0 + Power(tf,2)*v0 + 25*Power(t,2)*vf + 42*t*tf*vf - 5*Power(tf,2)*vf)) + 144*a0*Power(jMax,3)*(-Power(v0 - vf,3) + jMax*(v0 - vf)*(p0*(t + tf) - pf*(t + tf) - 5*Power(t,2)*v0 + Power(tf,2)*v0 + 5*Power(t,2)*vf + t*tf*vf) + Power(jMax,2)*t*(p0*(Power(t,2) + 2*t*tf - Power(tf,2)) + pf*(-Power(t,2) - 2*t*tf + Power(tf,2)) - t*(Power(t,2) + 4*t*tf - 2*Power(tf,2))*v0 + (Power(t,3) + 5*Power(t,2)*tf - Power(tf,3))*vf)) + 36*Power(a0,2)*Power(jMax,3)*(Power(jMax,2)*Power(t,3)*(2*t - tf)*tf - (19*t + tf)*Power(v0 - vf,2) - jMax*(pf*(9*Power(t,2) + 2*t*tf - Power(tf,2)) + p0*(-9*Power(t,2) - 2*t*tf + Power(tf,2)) + 4*t*(3*t - tf)*(t + 2*tf)*v0 + (-12*Power(t,3) - 29*Power(t,2)*tf + 6*t*Power(tf,2) + Power(tf,3))*vf)) + 72*Power(jMax,4)*(Power(jMax,2)*Power(t,2)*(2*t - tf)*tf*(p0 - pf - t*v0 + (t + tf)*vf) + Power(v0 - vf,2)*(2*p0 - 2*pf - t*v0 + tf*v0 + (t + tf)*vf) + jMax*(Power(p0,2)*(t - tf) + Power(pf,2)*(t - tf) + pf*(-Power(t,2) - 2*t*tf + Power(tf,2))*v0 + t*(-3*Power(t,2) + Power(tf,2))*Power(v0,2) + pf*(Power(t,2) + Power(tf,2))*vf + (6*Power(t,3) + Power(t,2)*tf - Power(tf,3))*v0*vf - Power(t,2)*(3*t + tf)*Power(vf,2) - p0*(2*pf*(t - tf) - 2*t*tf*v0 + Power(t,2)*(-v0 + vf) + Power(tf,2)*(v0 + vf)))))/(jMax*(-Power(a0,6) + 6*Power(a0,4)*jMax*(v0 - vf) - 36*Power(a0,2)*Power(jMax,2)*Power(v0 - vf,2) + 48*Power(a0,3)*Power(jMax,2)*(p0 - pf + tf*vf) - 144*a0*Power(jMax,3)*(v0 - vf)*(p0 - pf + tf*vf) + 72*Power(jMax,3)*(Power(v0 - vf,3) + jMax*Power(p0 - pf + tf*vf,2))));
            // profile.t[5] = 0.0;
            // profile.t[6] = (12*Power(a0,7) + Power(a0,6)*jMax*(101*t + 19*tf) + 6*Power(a0,5)*jMax*(45*jMax*Power(t,2) + 22*jMax*t*tf - 5*jMax*Power(tf,2) + 6*v0 - 6*vf) + 6*Power(a0,4)*Power(jMax,2)*(6*p0 - 6*pf + jMax*t*(41*Power(t,2) + 52*t*tf - 23*Power(tf,2)) - 7*tf*v0 + 11*t*(v0 - vf) + 13*tf*vf) + 12*Power(a0,3)*Power(jMax,2)*(Power(jMax,2)*Power(t,2)*(6*Power(t,2) + 22*t*tf - 11*Power(tf,2)) - 8*Power(v0 - vf,2) + jMax*(16*p0*t - 16*pf*t - 4*p0*tf + 4*pf*tf - 25*Power(t,2)*v0 - 26*t*tf*v0 + Power(tf,2)*v0 + 25*Power(t,2)*vf + 42*t*tf*vf - 5*Power(tf,2)*vf)) + 144*a0*Power(jMax,3)*(-Power(v0 - vf,3) + jMax*(v0 - vf)*(p0*(t + tf) - pf*(t + tf) - 5*Power(t,2)*v0 + Power(tf,2)*v0 + 5*Power(t,2)*vf + t*tf*vf) + Power(jMax,2)*t*(p0*(Power(t,2) + 2*t*tf - Power(tf,2)) + pf*(-Power(t,2) - 2*t*tf + Power(tf,2)) - t*(Power(t,2) + 4*t*tf - 2*Power(tf,2))*v0 + (Power(t,3) + 5*Power(t,2)*tf - Power(tf,3))*vf)) + 36*Power(a0,2)*Power(jMax,3)*(Power(jMax,2)*Power(t,3)*(2*t - tf)*tf - (19*t + tf)*Power(v0 - vf,2) - jMax*(pf*(9*Power(t,2) + 2*t*tf - Power(tf,2)) + p0*(-9*Power(t,2) - 2*t*tf + Power(tf,2)) + 4*t*(3*t - tf)*(t + 2*tf)*v0 + (-12*Power(t,3) - 29*Power(t,2)*tf + 6*t*Power(tf,2) + Power(tf,3))*vf)) + 72*Power(jMax,4)*(Power(jMax,2)*Power(t,2)*(2*t - tf)*tf*(p0 - pf - t*v0 + (t + tf)*vf) + Power(v0 - vf,2)*(2*p0 - 2*pf - t*v0 + tf*v0 + (t + tf)*vf) + jMax*(Power(p0,2)*(t - tf) + Power(pf,2)*(t - tf) + pf*(-Power(t,2) - 2*t*tf + Power(tf,2))*v0 + t*(-3*Power(t,2) + Power(tf,2))*Power(v0,2) + pf*(Power(t,2) + Power(tf,2))*vf + (6*Power(t,3) + Power(t,2)*tf - Power(tf,3))*v0*vf - Power(t,2)*(3*t + tf)*Power(vf,2) - p0*(2*pf*(t - tf) - 2*t*tf*v0 + Power(t,2)*(-v0 + vf) + Power(tf,2)*(v0 + vf)))))/(jMax*(-Power(a0,6) + 6*Power(a0,4)*jMax*(v0 - vf) - 36*Power(a0,2)*Power(jMax,2)*Power(v0 - vf,2) + 48*Power(a0,3)*Power(jMax,2)*(p0 - pf + tf*vf) - 144*a0*Power(jMax,3)*(v0 - vf)*(p0 - pf + tf*vf) + 72*Power(jMax,3)*(Power(v0 - vf,3) + jMax*Power(p0 - pf + tf*vf,2))));

            profile.t[3] = -(4*Power(a0,3) + 12*a0*jMax*(v0 + jMax*t*(3*t + h1)) + 3*Power(a0,2)*jMax*(8*t + h1) + 6*Power(jMax,2)*(2*p0 - 2*pf + 2*jMax*Power(t,3) + 4*t*v0 + jMax*Power(t,2)*h1 + v0*h1 + h1*vf))/(6.*jMax*(Power(a0,2) + 4*a0*jMax*t + 2*jMax*(jMax*Power(t,2) + v0)));
            profile.t[4] = Sqrt((-vf + vPlat)/jMax);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        // Find root of 6th order polynom
        std::array<double, 7> polynom;
        polynom[0] = 144*Power(jMax,6);
        polynom[1] = -144*Power(jMax,5)*(-5*a0 + jMax*tf);
        polynom[2] = 36*Power(jMax,4)*(39*Power(a0,2) - 16*a0*jMax*tf + 2*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf));
        polynom[3] = 24*Power(jMax,3)*(55*Power(a0,3) - 33*Power(a0,2)*jMax*tf - 6*Power(jMax,2)*(p0 - pf + tf*v0) + 12*a0*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf));
        polynom[4] = 6*Power(jMax,2)*(101*Power(a0,4) - 76*Power(a0,3)*jMax*tf - 48*a0*Power(jMax,2)*(p0 - pf + tf*v0) + 12*Power(jMax,2)*(2*jMax*tf*(p0 - pf + tf*v0) + 3*Power(v0 - vf,2)) + 60*Power(a0,2)*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf));
        polynom[5] = 12*a0*jMax*(11*Power(a0,4) - 10*Power(a0,3)*jMax*tf - 12*a0*Power(jMax,2)*(p0 - pf + tf*v0) + 12*Power(jMax,2)*(2*jMax*tf*(p0 - pf + tf*v0) + 3*Power(v0 - vf,2)) + 12*Power(a0,2)*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf));
        polynom[6] = 11*Power(a0,6) - 12*Power(a0,5)*jMax*tf - 24*Power(a0,3)*Power(jMax,2)*(p0 - pf + tf*v0) + 36*Power(a0,2)*Power(jMax,2)*(2*jMax*tf*(p0 - pf + tf*v0) + 3*Power(v0 - vf,2)) + 72*Power(jMax,3)*(jMax*Power(p0 - pf + tf*v0,2) + Power(v0 - vf,3)) + 18*Power(a0,4)*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf);
    
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
        double tz_max = std::min<double>(tf, (tf - a0/jMax) / 2);
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

            if (tz <= 0.0 || tz >= tz_max) {
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
            double tz = Roots::shrinkInterval(polynom, lower, upper, 1e-14);

            double vPlat = Power(a0,2)/(2.*jMax) + 2*a0*tz + jMax*Power(tz,2) + v0;
            // std::cout << "BE CAREFUL " << vPlat << std::endl;

            profile.t[0] = tz;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + tz;
            profile.t[3] = tf - (2*tz + a0/jMax + 2 * Sqrt((vf - vPlat)/jMax));
            profile.t[4] = Sqrt((vf - vPlat)/jMax); 
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool RuckigStep2::time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    if (std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        profile.t[0] = (Power(aMax,2)*Power(tf,2) - Power(v0 - vf,2) + 2*aMax*(2*p0 - 2*pf + tf*(v0 + vf)))/(2.*Power(aMax,2)*tf);
        profile.t[1] = -(Power(aMax,2)*Power(tf,2) - 2*Power(v0 - vf,2) + aMax*(8*p0 - 8*pf + 5*tf*v0 + 3*tf*vf))/(2.*Power(aMax,2)*tf);
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = -(Power(aMax,2)*Power(tf,2) - 2*Power(v0 - vf,2) + aMax*(8*p0 - 8*pf + 3*tf*v0 + 5*tf*vf))/(2.*Power(aMax,2)*tf);
        profile.t[6] = profile.t[0];
        jMax = aMax/profile.t[0];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        return profile.check(tf, pf, vf, af, vMax, aMax);
    }

    double h3 = Sqrt(aMax*(6*Power(aMax,5)*Power(tf,2) + 3*Power(af,4)*(2*p0 - 2*pf + tf*(aMax*tf + 2*v0)) + 3*Power(a0,4)*(2*p0 - 2*pf + tf*(aMax*tf + 2*vf)) + 6*Power(af,2)*aMax*(Power(aMax,2)*Power(tf,2) + 2*Power(v0 - vf,2) + aMax*(-4*p0 + 4*pf + tf*v0 - 5*tf*vf)) + 4*Power(af,3)*(2*Power(aMax,2)*Power(tf,2) + Power(v0 - vf,2) + 2*aMax*(p0 - pf + 2*tf*v0 - tf*vf)) - 4*Power(a0,3)*(2*Power(aMax,2)*Power(tf,2) + Power(v0 - vf,2) + 2*aMax*(p0 - pf - tf*v0 + 2*tf*vf)) - 12*af*Power(aMax,2)*(-Power(v0 - vf,2) + aMax*(4*p0 - 4*pf + tf*v0 + 3*tf*vf)) - 6*Power(a0,2)*(Power(af,2)*(2*p0 - 2*pf + tf*(v0 + vf)) + 2*af*aMax*(2*p0 - 2*pf + tf*(v0 + vf)) - aMax*(Power(aMax,2)*Power(tf,2) + 2*Power(v0 - vf,2) + aMax*(-4*p0 + 4*pf - 5*tf*v0 + tf*vf))) + 12*a0*aMax*(Power(af,2)*(2*p0 - 2*pf + tf*(v0 + vf)) + 2*af*aMax*(2*p0 - 2*pf + tf*(v0 + vf)) + aMax*(-Power(v0 - vf,2) + aMax*(4*p0 - 4*pf + 3*tf*v0 + tf*vf)))));
    double h1 = (3*Power(a0,2)*aMax*tf + 3*Power(af,2)*aMax*tf - 6*a0*Power(aMax,2)*tf + 6*af*Power(aMax,2)*tf + 6*Power(aMax,3)*tf - 3*Power(a0,2)*v0 + 3*Power(af,2)*v0 + 6*a0*aMax*v0 + 6*af*aMax*v0 + 3*Power(a0,2)*vf - 3*Power(af,2)*vf - 6*a0*aMax*vf - 6*af*aMax*vf - Sqrt(6)*h3);
    double h5 = (3*Power(a0,3) - 3*Power(af,3) + Power(a0,2)*(3*af - 4*aMax) - 4*Power(af,2)*aMax + 12*af*Power(aMax,2) + 24*Power(aMax,3) - a0*(3*Power(af,2) + 16*af*aMax + 12*Power(aMax,2)));
    double h2 = ((a0 - af)*h5);
    double h4 = (2.*(a0 - af)*aMax*h5);

    profile.t[0] = (2*(a0 - aMax)*h1)/h2;
    profile.t[1] = (-6*Power(a0,2)*Power(af,2)*aMax*tf + 6*Power(af,4)*aMax*tf - 4*Power(a0,3)*Power(aMax,2)*tf - 12*Power(a0,2)*af*Power(aMax,2)*tf + 16*Power(af,3)*Power(aMax,2)*tf + 18*Power(a0,2)*Power(aMax,3)*tf + 30*Power(af,2)*Power(aMax,3)*tf - 36*a0*Power(aMax,4)*tf + 36*af*Power(aMax,4)*tf + 48*Power(aMax,5)*tf + 4*Power(a0,3)*aMax*v0 - 12*a0*Power(af,2)*aMax*v0 + 8*Power(af,3)*aMax*v0 - 24*Power(a0,2)*Power(aMax,2)*v0 - 24*a0*af*Power(aMax,2)*v0 + 48*Power(af,2)*Power(aMax,2)*v0 + 24*a0*Power(aMax,3)*v0 + 72*af*Power(aMax,3)*v0 - 4*Power(a0,3)*aMax*vf + 12*a0*Power(af,2)*aMax*vf - 8*Power(af,3)*aMax*vf + 24*Power(a0,2)*Power(aMax,2)*vf + 24*a0*af*Power(aMax,2)*vf - 48*Power(af,2)*Power(aMax,2)*vf - 24*a0*Power(aMax,3)*vf - 72*af*Power(aMax,3)*vf + Sqrt(6)*Power(a0,2)*h3 - Sqrt(6)*Power(af,2)*h3 + 2*Sqrt(6)*a0*aMax*h3 - 2*Sqrt(6)*af*aMax*h3 - 8*Sqrt(6)*Power(aMax,2)*h3)/h4;
    profile.t[2] = (-2*aMax*h1)/h2;
    profile.t[3] = 0;
    profile.t[4] = profile.t[2];
    profile.t[5] = (6*Power(a0,4)*aMax*tf - 6*Power(a0,2)*Power(af,2)*aMax*tf - 16*Power(a0,3)*Power(aMax,2)*tf + 12*a0*Power(af,2)*Power(aMax,2)*tf + 4*Power(af,3)*Power(aMax,2)*tf + 30*Power(a0,2)*Power(aMax,3)*tf + 18*Power(af,2)*Power(aMax,3)*tf - 36*a0*Power(aMax,4)*tf + 36*af*Power(aMax,4)*tf + 48*Power(aMax,5)*tf + 8*Power(a0,3)*aMax*v0 - 12*Power(a0,2)*af*aMax*v0 + 4*Power(af,3)*aMax*v0 - 48*Power(a0,2)*Power(aMax,2)*v0 + 24*a0*af*Power(aMax,2)*v0 + 24*Power(af,2)*Power(aMax,2)*v0 + 72*a0*Power(aMax,3)*v0 + 24*af*Power(aMax,3)*v0 - 8*Power(a0,3)*aMax*vf + 12*Power(a0,2)*af*aMax*vf - 4*Power(af,3)*aMax*vf + 48*Power(a0,2)*Power(aMax,2)*vf - 24*a0*af*Power(aMax,2)*vf - 24*Power(af,2)*Power(aMax,2)*vf - 72*a0*Power(aMax,3)*vf - 24*af*Power(aMax,3)*vf - Sqrt(6)*Power(a0,2)*h3 + Sqrt(6)*Power(af,2)*h3 + 2*Sqrt(6)*a0*aMax*h3 - 2*Sqrt(6)*af*aMax*h3 - 8*Sqrt(6)*Power(aMax,2)*h3)/h4;
    profile.t[6] = (-2*(af + aMax)*h1)/h2;
    jMax = -(-12*Power(a0,2)*aMax*tf - 12*Power(af,2)*aMax*tf + 24*a0*Power(aMax,2)*tf - 24*af*Power(aMax,2)*tf - 24*Power(aMax,3)*tf + 12*Power(a0,2)*v0 - 12*Power(af,2)*v0 - 24*a0*aMax*v0 - 24*af*aMax*v0 - 12*Power(a0,2)*vf + 12*Power(af,2)*vf + 24*a0*aMax*vf + 24*af*aMax*vf - 4*Sqrt(9*Power(2*Power(aMax,3)*tf + Power(af,2)*(aMax*tf + v0 - vf) + 2*af*aMax*(aMax*tf + v0 - vf) + Power(a0,2)*(aMax*tf - v0 + vf) - 2*a0*aMax*(aMax*tf - v0 + vf),2) + 3*(a0 - af)*h5*(Power(aMax,2)*Power(tf,2) - Power(v0 - vf,2) + 2*aMax*(2*p0 - 2*pf + tf*(v0 + vf)))))/(24.*(Power(aMax,2)*Power(tf,2) - Power(v0 - vf,2) + 2*aMax*(2*p0 - 2*pf + tf*(v0 + vf))));

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(tf, pf, vf, af, vMax, aMax);
}

bool RuckigStep2::time_up_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    // a3 != 0
    
    // Case UDDU, Solution 2
    {
        double h1 = Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) - Power(af,3) + 3*Power(a0,2)*aMax + 3*a0*Power(aMax,2) + 3*Power(aMax,2)*jMax*tf - 3*af*aMax*(aMax - 2*jMax*tf) - 3*Power(af,2)*(aMax - jMax*tf) - 3*Power(jMax,2)*(2*p0 - 2*pf + aMax*Power(tf,2) + 2*tf*vf),2) - 3*(Power(a0,2) + Power(af,2) + 2*a0*aMax + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*v0 - jMax*vf))*(Power(a0,4) + 3*Power(af,4) + 4*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 6*Power(af,2)*(Power(aMax,2) + 2*jMax*(v0 - vf)) + 12*jMax*(-2*aMax*jMax*(p0 - pf + tf*v0) + Power(aMax,2)*(v0 - vf) + jMax*Power(v0 - vf,2)) + 24*af*aMax*jMax*(v0 - vf) - 4*a0*(Power(af,3) + 3*af*aMax*(aMax - 2*jMax*tf) + 3*Power(af,2)*(aMax - jMax*tf) + 3*jMax*(-(Power(aMax,2)*tf) + jMax*(2*p0 - 2*pf + aMax*Power(tf,2) + 2*tf*vf))))));
        double h2 = (6.*Power(jMax,2)*(Power(a0,2) + Power(af,2) + 2*a0*aMax + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*v0 - jMax*vf)));

        profile.t[0] = 0;
        profile.t[1] = 0;
        profile.t[2] = -(-2*Power(a0,3)*jMax + 2*Power(af,3)*jMax - 6*Power(a0,2)*aMax*jMax + 6*Power(af,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 6*af*Power(aMax,2)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf - 6*Power(af,2)*Power(jMax,2)*tf - 12*af*aMax*Power(jMax,2)*tf - 6*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*Power(jMax,3)*tf*vf - h1)/h2;
        profile.t[3] = -(4*Power(a0,3)*jMax + 6*a0*Power(af,2)*jMax + 2*Power(af,3)*jMax + 12*Power(a0,2)*aMax*jMax + 12*a0*af*aMax*jMax + 12*Power(af,2)*aMax*jMax + 18*a0*Power(aMax,2)*jMax + 18*af*Power(aMax,2)*jMax + 12*Power(aMax,3)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf - 6*Power(af,2)*Power(jMax,2)*tf - 12*a0*aMax*Power(jMax,2)*tf - 12*af*aMax*Power(jMax,2)*tf - 18*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*a0*Power(jMax,2)*v0 + 12*aMax*Power(jMax,2)*v0 - 12*a0*Power(jMax,2)*vf - 12*aMax*Power(jMax,2)*vf + 12*Power(jMax,3)*tf*vf + h1)/h2;
        profile.t[4] = -(-4*Power(a0,3)*jMax - 6*a0*Power(af,2)*jMax - 2*Power(af,3)*jMax - 12*Power(a0,2)*aMax*jMax - 12*a0*af*aMax*jMax - 12*Power(af,2)*aMax*jMax - 18*a0*Power(aMax,2)*jMax - 18*af*Power(aMax,2)*jMax - 12*Power(aMax,3)*jMax - 12*Power(jMax,3)*p0 + 12*Power(jMax,3)*pf + 6*Power(af,2)*Power(jMax,2)*tf + 12*a0*aMax*Power(jMax,2)*tf + 12*af*aMax*Power(jMax,2)*tf + 18*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) - 12*a0*Power(jMax,2)*v0 - 12*aMax*Power(jMax,2)*v0 + 12*a0*Power(jMax,2)*vf + 12*aMax*Power(jMax,2)*vf - 12*Power(jMax,3)*tf*vf + h1)/h2;
        profile.t[5] = -(2*Power(a0,3)*jMax + 6*Power(a0,2)*af*jMax + 4*Power(af,3)*jMax + 12*Power(a0,2)*aMax*jMax + 12*a0*af*aMax*jMax + 12*Power(af,2)*aMax*jMax + 18*a0*Power(aMax,2)*jMax + 18*af*Power(aMax,2)*jMax + 12*Power(aMax,3)*jMax - 12*Power(jMax,3)*p0 + 12*Power(jMax,3)*pf - 6*Power(a0,2)*Power(jMax,2)*tf - 12*a0*aMax*Power(jMax,2)*tf - 12*af*aMax*Power(jMax,2)*tf - 18*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*af*Power(jMax,2)*v0 + 12*aMax*Power(jMax,2)*v0 - 12*Power(jMax,3)*tf*v0 - 12*af*Power(jMax,2)*vf - 12*aMax*Power(jMax,2)*vf - h1)/h2;
        profile.t[6] = (af + aMax)/jMax;

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(tf, pf, vf, af, vMax, aMax)) {
            return true;
        }
    }

    // Case UDUD, Solution 1
    {
        double h1 = Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(-Power(a0,3) + Power(af,3) + 3*Power(a0,2)*aMax - 3*a0*Power(aMax,2) + 3*af*aMax*(aMax - 2*jMax*tf) - 3*Power(af,2)*(aMax - jMax*tf) + 3*jMax*(Power(aMax,2)*tf + jMax*(2*p0 - 2*pf - aMax*Power(tf,2) + 2*tf*vf)),2) - 3*(Power(a0,2) - Power(af,2) - 2*a0*aMax + 2*af*aMax + 2*jMax*(aMax*tf + v0 - vf))*(Power(a0,4) + 3*Power(af,4) - 4*Power(a0,3)*aMax - 8*Power(af,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 24*af*aMax*jMax*(v0 - vf) + 12*jMax*(2*aMax*jMax*(p0 - pf + tf*v0) + jMax*Power(v0 - vf,2) + Power(aMax,2)*(-v0 + vf)) + 6*Power(af,2)*(Power(aMax,2) + 2*jMax*(-v0 + vf)) - 4*a0*(Power(af,3) + 3*af*aMax*(aMax - 2*jMax*tf) - 3*Power(af,2)*(aMax - jMax*tf) + 3*jMax*(Power(aMax,2)*tf + jMax*(2*p0 - 2*pf - aMax*Power(tf,2) + 2*tf*vf))))));
        double h2 = (6.*Power(jMax,2)*(Power(a0,2) - Power(af,2) - 2*a0*aMax + 2*af*aMax + 2*jMax*(aMax*tf + v0 - vf)));

        profile.t[0] = 0;
        profile.t[1] = 0;
        profile.t[2] = -(-2*Power(a0,3)*jMax + 2*Power(af,3)*jMax + 6*Power(a0,2)*aMax*jMax - 6*Power(af,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 6*af*Power(aMax,2)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf + 6*Power(af,2)*Power(jMax,2)*tf - 12*af*aMax*Power(jMax,2)*tf + 6*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) + 12*Power(jMax,3)*tf*vf + h1)/h2;
        profile.t[3] = (h1)/(3.*Power(jMax,2)*(Power(a0,2) - Power(af,2) - 2*a0*aMax + 2*af*aMax + 2*jMax*(aMax*tf + v0 - vf)));
        profile.t[4] = -(4*Power(a0,3)*jMax - 6*a0*Power(af,2)*jMax + 2*Power(af,3)*jMax - 12*Power(a0,2)*aMax*jMax + 12*a0*af*aMax*jMax + 6*a0*Power(aMax,2)*jMax - 6*af*Power(aMax,2)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf + 6*Power(af,2)*Power(jMax,2)*tf + 12*a0*aMax*Power(jMax,2)*tf - 12*af*aMax*Power(jMax,2)*tf - 6*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) + 12*a0*Power(jMax,2)*v0 - 12*aMax*Power(jMax,2)*v0 - 12*a0*Power(jMax,2)*vf + 12*aMax*Power(jMax,2)*vf + 12*Power(jMax,3)*tf*vf + h1)/h2;
        profile.t[5] = (Power(a0,3) - Power(af,3) + 3*Power(a0,2)*(af - 2*aMax + jMax*tf) + 3*Power(af,2)*(2*aMax + jMax*tf) - 3*a0*(Power(af,2) - 2*(Power(aMax,2) + jMax*(v0 - vf))) - 6*af*(Power(aMax,2) + jMax*(-v0 + vf)) + 6*jMax*(-(aMax*(aMax*tf + 2*v0 - 2*vf)) + jMax*(2*p0 - 2*pf + tf*(v0 + vf))))/(3.*jMax*(Power(a0,2) - Power(af,2) - 2*a0*aMax + 2*af*aMax + 2*jMax*(aMax*tf + v0 - vf)));
        profile.t[6] = (-af + aMax)/jMax;

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
        if (profile.check(tf, pf, vf, af, vMax, aMax)) {
            return true;
        }
    }
    return false;
}

bool RuckigStep2::time_up_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    // a3 != 0

    double h1 = Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) + 2*Power(af,3) - 6*Power(af,2)*aMax + 9*af*Power(aMax,2) - 6*Power(aMax,3) - 6*af*aMax*jMax*tf + 9*Power(aMax,2)*jMax*tf + 3*a0*aMax*(-2*af + 3*aMax - 2*jMax*tf) + 3*Power(a0,2)*(af - 2*aMax + jMax*tf) - 6*Power(jMax,2)*(p0 - pf + tf*v0) + 6*af*jMax*(-v0 + vf) - 3*aMax*jMax*(jMax*Power(tf,2) - 2*v0 + 2*vf),2) - 9*Power(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf - jMax*v0 + jMax*vf),3)));
    double h2 = (6.*Power(jMax,2)*(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf - jMax*v0 + jMax*vf)));

    // Solution 1
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = -(-4*Power(a0,3)*jMax - 6*a0*Power(af,2)*jMax - 2*Power(af,3)*jMax + 12*Power(a0,2)*aMax*jMax + 12*a0*af*aMax*jMax + 12*Power(af,2)*aMax*jMax - 18*a0*Power(aMax,2)*jMax - 18*af*Power(aMax,2)*jMax + 12*Power(aMax,3)*jMax - 12*Power(jMax,3)*p0 + 12*Power(jMax,3)*pf - 6*Power(af,2)*Power(jMax,2)*tf + 12*a0*aMax*Power(jMax,2)*tf + 12*af*aMax*Power(jMax,2)*tf - 18*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*a0*Power(jMax,2)*v0 - 12*aMax*Power(jMax,2)*v0 - 12*a0*Power(jMax,2)*vf + 12*aMax*Power(jMax,2)*vf - 12*Power(jMax,3)*tf*vf - h1)/h2;
    profile.t[2] = -(2*Power(a0,3)*jMax + 6*Power(a0,2)*af*jMax + 4*Power(af,3)*jMax - 12*Power(a0,2)*aMax*jMax - 12*a0*af*aMax*jMax - 12*Power(af,2)*aMax*jMax + 18*a0*Power(aMax,2)*jMax + 18*af*Power(aMax,2)*jMax - 12*Power(aMax,3)*jMax - 12*Power(jMax,3)*p0 + 12*Power(jMax,3)*pf + 6*Power(a0,2)*Power(jMax,2)*tf - 12*a0*aMax*Power(jMax,2)*tf - 12*af*aMax*Power(jMax,2)*tf + 18*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) - 12*af*Power(jMax,2)*v0 + 12*aMax*Power(jMax,2)*v0 - 12*Power(jMax,3)*tf*v0 + 12*af*Power(jMax,2)*vf - 12*aMax*Power(jMax,2)*vf + h1)/h2;
    profile.t[3] = -(-2*Power(a0,3)*jMax - 6*Power(a0,2)*af*jMax - 4*Power(af,3)*jMax + 12*Power(a0,2)*aMax*jMax + 12*a0*af*aMax*jMax + 12*Power(af,2)*aMax*jMax - 18*a0*Power(aMax,2)*jMax - 18*af*Power(aMax,2)*jMax + 12*Power(aMax,3)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf - 6*Power(a0,2)*Power(jMax,2)*tf + 12*a0*aMax*Power(jMax,2)*tf + 12*af*aMax*Power(jMax,2)*tf - 18*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*af*Power(jMax,2)*v0 - 12*aMax*Power(jMax,2)*v0 + 12*Power(jMax,3)*tf*v0 - 12*af*Power(jMax,2)*vf + 12*aMax*Power(jMax,2)*vf + h1)/h2;
    profile.t[4] = -(-2*Power(a0,3)*jMax + 2*Power(af,3)*jMax + 6*Power(a0,2)*aMax*jMax - 6*Power(af,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 6*af*Power(aMax,2)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf - 6*Power(a0,2)*Power(jMax,2)*tf + 12*a0*aMax*Power(jMax,2)*tf - 6*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*Power(jMax,3)*tf*v0 - h1)/h2;
    profile.t[5] = 0;
    profile.t[6] = 0;

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    if (profile.check(tf, pf, vf, af, vMax, aMax)) {
        return true;
    }

    return false;
}

bool RuckigStep2::time_up_none(Profile& profile, double vMax, double aMax, double jMax) {
    if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        profile.t[0] = tf/4;
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        double jMaxNew = (-32*(p0 - pf))/Power(tf,3);

        profile.set(p0, v0, a0, {jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
        return profile.check(tf, pf, vf, af, vMax, aMax);
    }
    
    if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON) {
        double h1 = Sqrt(Power(tf,2)*(Power(tf,2)*Power(vf,2) + 4*Power(2*p0 - 2*pf + tf*vf,2)));
        profile.t[0] = (4*p0 - 4*pf + 3*tf*vf + h1/tf)/(4.*vf);
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        double jMaxNew = (4*(-4*p0*tf + 4*pf*tf - 2*Power(tf,2)*vf + h1))/Power(tf,4);

        profile.set(p0, v0, a0, {jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
        return profile.check(tf, pf, vf, af, vMax, aMax);
    }

    if (std::abs(a0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON) {
        // Solution 1
        {
            profile.t[0] = (-4*p0 + 4*pf - tf*v0 + Sqrt(Power(tf,2)*(Power(tf,2)*Power(v0,2) + 4*Power(2*p0 - 2*pf + tf*v0,2)))/tf)/(4.*v0);
            profile.t[1] = 0;
            profile.t[2] = profile.t[0];
            profile.t[3] = 0;
            profile.t[4] = (4*p0 - 4*pf + 3*tf*v0 - Sqrt(Power(tf,2)*(Power(tf,2)*Power(v0,2) + 4*Power(2*p0 - 2*pf + tf*v0,2)))/tf)/(4.*v0);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            double jMaxNew = (-4*(4*p0*tf - 4*pf*tf + 2*Power(tf,2)*v0 + Sqrt(Power(tf,2)*(Power(tf,2)*Power(v0,2) + 4*Power(2*p0 - 2*pf + tf*v0,2)))))/Power(tf,4);

            profile.set(p0, v0, a0, {jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profiles with a3 != 0, Solution UDDU
    {
        // First acc, then constant
        {
            std::array<double, 5> polynom;
            polynom[0] = 12*Power(jMax,4);
            polynom[1] = 24*Power(jMax,3)*(a0 - jMax*tf);
            polynom[2] = 12*Power(jMax,2)*(2*Power(a0,2) - 4*a0*jMax*tf + jMax*(jMax*Power(tf,2) - 2*v0 + 2*vf));
            polynom[3] = 8*jMax*(Power(a0,3) - 3*Power(a0,2)*jMax*tf + 3*a0*Power(jMax,2)*Power(tf,2) - 6*Power(jMax,2)*(p0 - pf + tf*vf));
            polynom[4] = Power(a0,4) - 4*Power(a0,3)*jMax*tf + 6*Power(a0,2)*Power(jMax,2)*Power(tf,2) - 24*a0*Power(jMax,2)*(p0 - pf + tf*vf) + 12*Power(jMax,2)*(Power(v0 - vf,2) + jMax*tf*(2*p0 - 2*pf + tf*(v0 + vf)));
            auto roots = Roots::solveQuart(polynom);

            for (double t: roots) {
                if (t < 0.0 || t > tf) {
                    continue;
                }

                profile.t[0] = t;
                profile.t[1] = 0;
                profile.t[2] = (-Power(a0,3) + 3*Power(a0,2)*jMax*(4*t + tf) + 6*a0*jMax*(jMax*(3*Power(t,2) - 4*t*tf + Power(tf,2)) + 3*(v0 - vf)) + 6*Power(jMax,2)*(-8*p0 + 8*pf + 2*jMax*Power(t,3) - 3*jMax*Power(t,2)*tf + jMax*t*Power(tf,2) - 2*t*v0 - 3*tf*v0 + 2*t*vf - 5*tf*vf))/(6.*jMax*(-Power(a0,2) + 2*a0*jMax*tf + jMax*(jMax*Power(tf,2) + 4*v0 - 4*vf)));
                profile.t[3] = -(a0/jMax) - 2*t + tf;
                profile.t[4] = -(5*Power(a0,3) - 9*Power(a0,2)*jMax*(-2*t + tf) + 6*a0*jMax*(3*jMax*t*(t - 2*tf) - v0 + vf) + 6*Power(jMax,2)*(-8*p0 + 8*pf + 2*jMax*Power(t,3) - 3*jMax*Power(t,2)*tf - 6*t*v0 - 3*tf*v0 + 6*t*vf - 5*tf*vf))/(6.*jMax*(-Power(a0,2) + 2*a0*jMax*tf + jMax*(jMax*Power(tf,2) + 4*v0 - 4*vf)));
                profile.t[5] = 0;
                profile.t[6] = 0;

                // std::cout << t << std::endl;
                // std::cout << profile.t[0] << std::endl;
                // std::cout << profile.t[1] << std::endl;
                // std::cout << profile.t[2] << std::endl;
                // std::cout << profile.t[3] << std::endl;
                // std::cout << profile.t[4] << " " << -(5*Power(a0,3) - 9*Power(a0,2)*jMax*(-2*t + tf) + 6*a0*jMax*(3*jMax*t*(t - 2*tf) - v0 + vf) + 6*Power(jMax,2)*(-8*p0 + 8*pf + 2*jMax*Power(t,3) - 3*jMax*Power(t,2)*tf - 6*t*v0 - 3*tf*v0 + 6*t*vf - 5*tf*vf)) << std::endl;
                // std::cout << "---" << std::endl;

                profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
                if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                    return true;
                }
            }
        }

        // First constant, then acc
        {
            std::array<double, 5> polynom;
            polynom[0] = 18*Power(jMax,4)*(Power(a0,2) + 2*a0*jMax*tf - jMax*(jMax*Power(tf,2) - 4*v0 + 4*vf));
            polynom[1] = -12*Power(jMax,3)*(2*Power(a0,3) + 6*Power(a0,2)*jMax*tf + 3*a0*jMax*(jMax*Power(tf,2) + 2*v0 - 2*vf) - 3*Power(jMax,2)*(4*p0 - 4*pf + jMax*Power(tf,3) - 2*tf*v0 + 6*tf*vf));
            polynom[2] = 18*Power(jMax,2)*(Power(a0,4) + 4*Power(a0,3)*jMax*tf + 4*Power(a0,2)*jMax*(jMax*Power(tf,2) + v0 - vf) - 2*a0*Power(jMax,2)*(6*p0 - 6*pf + jMax*Power(tf,3) - 4*tf*v0 + 10*tf*vf) - Power(jMax,2)*(Power(jMax,2)*Power(tf,4) - 8*Power(v0 - vf,2) + 4*jMax*tf*(3*p0 - 3*pf + tf*v0 + 2*tf*vf)));
            polynom[3] = -6*jMax*(Power(a0,5) + 5*Power(a0,4)*jMax*tf + 4*Power(a0,3)*jMax*(2*jMax*Power(tf,2) + v0 - vf) - 12*Power(a0,2)*Power(jMax,2)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf) - 6*a0*Power(jMax,2)*(Power(jMax,2)*Power(tf,4) - 2*Power(v0 - vf,2) + 8*jMax*tf*(p0 - pf + tf*vf)) - 12*Power(jMax,3)*(jMax*Power(tf,2)*(p0 - pf + tf*v0) + (v0 - vf)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf)));
            polynom[4] = Power(a0,6) + 6*Power(a0,5)*jMax*tf - 72*Power(jMax,3)*(jMax*Power(p0 - pf + tf*v0,2) - Power(v0 - vf,3)) + 6*Power(a0,4)*jMax*(2*jMax*Power(tf,2) + v0 - vf) - 24*Power(a0,3)*Power(jMax,2)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf) - 18*Power(a0,2)*Power(jMax,2)*(Power(jMax,2)*Power(tf,4) - 2*Power(v0 - vf,2) + 8*jMax*tf*(p0 - pf + tf*vf)) - 72*a0*Power(jMax,3)*(jMax*Power(tf,2)*(p0 - pf + tf*v0) + (v0 - vf)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf));
            auto roots = Roots::solveQuart(polynom);

            for (double t: roots) {
                if (t < 0.0 || t > tf) {
                    continue;
                }

                profile.t[0] = 0;
                profile.t[1] = 0;
                profile.t[2] = t;
                profile.t[3] = (-Power(a0,4) + 2*Power(a0,3)*jMax*(t - 2*tf) - 3*Power(a0,2)*jMax*(jMax*(2*Power(t,2) - 2*t*tf + Power(tf,2)) + 2*(v0 - vf)) + 6*a0*Power(jMax,2)*(4*p0 - 4*pf + tf*(jMax*(-2*Power(t,2) + t*tf + Power(tf,2)) - 2*v0 + 6*vf)) + 6*Power(jMax,2)*(Power(jMax,2)*t*(t - tf)*Power(tf,2) - 4*Power(v0 - vf,2) + jMax*(-8*p0*t + 8*pf*t + 4*p0*tf - 4*pf*tf - 4*Power(t,2)*v0 + 3*Power(tf,2)*v0 + 4*Power(t,2)*vf - 8*t*tf*vf + Power(tf,2)*vf)))/(jMax*(Power(a0,3) + 3*Power(a0,2)*jMax*tf + 6*a0*jMax*(v0 - vf) + 6*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf))));
                profile.t[4] = (Power(a0,4) - 2*Power(a0,3)*jMax*(t - 2*tf) + 3*Power(a0,2)*jMax*(jMax*Power(t - tf,2) + 2*(v0 - vf)) + 3*a0*Power(jMax,2)*(-2*p0 + 2*pf + 2*jMax*Power(t,2)*tf - jMax*t*Power(tf,2) - jMax*Power(tf,3) - 2*t*v0 + 4*tf*v0 + 2*t*vf - 6*tf*vf) - 3*Power(jMax,2)*(Power(jMax,2)*t*(t - tf)*Power(tf,2) - 4*Power(v0 - vf,2) + 2*jMax*(-2*p0*t + 2*pf*t + p0*tf - pf*tf - 2*Power(t,2)*v0 + t*tf*v0 + Power(tf,2)*v0 + 2*Power(t,2)*vf - 3*t*tf*vf)))/(jMax*(Power(a0,3) + 3*Power(a0,2)*jMax*tf + 6*a0*jMax*(v0 - vf) + 6*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf))));
                profile.t[5] = 0;
                profile.t[6] = (Power(a0,3)*(-t + tf) + 3*Power(a0,2)*jMax*(Power(t,2) - t*tf + Power(tf,2)) + 3*a0*jMax*(-6*p0 + 6*pf - tf*(jMax*(-2*Power(t,2) + t*tf + Power(tf,2)) - 2*v0 + 8*vf)) - 3*jMax*(Power(jMax,2)*t*(t - tf)*Power(tf,2) - 4*Power(v0 - vf,2) + jMax*(-8*p0*t + 8*pf*t + 2*p0*tf - 2*pf*tf - 4*Power(t,2)*v0 + 2*Power(tf,2)*v0 + 4*Power(t,2)*vf - 8*t*tf*vf)))/(Power(a0,3) + 3*Power(a0,2)*jMax*tf + 6*a0*jMax*(v0 - vf) + 6*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf)));

                profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
                if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                    return true;
                }
            }       
        }
    }

    // Profiles with a3 != 0, Solution UDUD
    {
        // First constant, then acc
        {
            std::array<double, 5> polynom;
            polynom[0] = 18*Power(jMax,4)*Power(a0 - jMax*tf,2);
            polynom[1] = 36*Power(jMax,4)*(a0 - jMax*tf)*(jMax*Power(tf,2) + 2*v0 - 2*vf);
            polynom[2] = 6*Power(jMax,2)*(-2*Power(a0,4) + 8*Power(a0,3)*jMax*tf - 12*Power(a0,2)*Power(jMax,2)*Power(tf,2) + 6*a0*Power(jMax,2)*(2*p0 - 2*pf + jMax*Power(tf,3) + 2*tf*vf) + 3*Power(jMax,2)*(Power(jMax,2)*Power(tf,4) + 4*Power(v0 - vf,2) + 4*jMax*tf*(-p0 + pf + tf*v0 - 2*tf*vf)));
            polynom[3] = 6*jMax*(Power(a0,5) - 5*Power(a0,4)*jMax*tf + 4*Power(a0,3)*jMax*(2*jMax*Power(tf,2) - v0 + vf) - 12*Power(a0,2)*Power(jMax,2)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf) - 12*Power(jMax,3)*(jMax*Power(tf,2)*(p0 - pf + tf*v0) + (v0 - vf)*(-2*p0 + 2*pf + tf*v0 - 3*tf*vf)) - 6*a0*Power(jMax,2)*(Power(jMax,2)*Power(tf,4) - 2*Power(v0 - vf,2) - 8*jMax*tf*(p0 - pf + tf*vf)));
            polynom[4] = -Power(a0,6) + 6*Power(a0,5)*jMax*tf + 72*Power(jMax,3)*(jMax*Power(p0 - pf + tf*v0,2) + Power(v0 - vf,3)) - 6*Power(a0,4)*jMax*(2*jMax*Power(tf,2) - v0 + vf) + 24*Power(a0,3)*Power(jMax,2)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf) + 72*a0*Power(jMax,3)*(jMax*Power(tf,2)*(p0 - pf + tf*v0) + (v0 - vf)*(-2*p0 + 2*pf + tf*v0 - 3*tf*vf)) + 18*Power(a0,2)*Power(jMax,2)*(Power(jMax,2)*Power(tf,4) - 2*Power(v0 - vf,2) - 8*jMax*tf*(p0 - pf + tf*vf));
            auto roots = Roots::solveQuart(polynom);

            for (double t: roots) {
                if (t < 0.0 || t > tf) {
                    continue;
                }

                profile.t[0] = 0;
                profile.t[1] = 0;
                profile.t[2] = t;
                profile.t[3] = -((Power(a0,5) - Power(a0,4)*jMax*(4*t + 5*tf) + 2*Power(a0,3)*jMax*(jMax*(3*Power(t,2) + 8*t*tf + 4*Power(tf,2)) - 2*v0 + 2*vf) + 6*Power(a0,2)*Power(jMax,2)*(-10*p0 + 10*pf + 2*jMax*Power(t,3) - 3*jMax*Power(t,2)*tf - 3*jMax*t*Power(tf,2) + jMax*Power(tf,3) + 2*t*v0 + 2*tf*v0 - 2*t*vf - 12*tf*vf) + 12*Power(jMax,3)*(Power(jMax,2)*t*Power(t - tf,2)*Power(tf,2) + (v0 - vf)*(2*p0 - 2*pf + 2*t*v0 - 3*tf*v0 - 2*t*vf + 5*tf*vf) - jMax*tf*(2*p0*(t + 2*tf) - 2*pf*(t + 2*tf) + 3*Power(t,2)*v0 - 3*t*tf*v0 + 3*Power(tf,2)*v0 - 3*Power(t,2)*vf + 5*t*tf*vf + Power(tf,2)*vf)) - 12*a0*Power(jMax,2)*(Power(jMax,2)*Power(t - tf,2)*tf*(2*t + tf) - 3*Power(v0 - vf,2) - jMax*(2*p0*(t + 5*tf) - 2*pf*(t + 5*tf) + 3*Power(t,2)*v0 - 2*t*tf*v0 + 2*Power(tf,2)*v0 - 3*Power(t,2)*vf + 4*t*tf*vf + 8*Power(tf,2)*vf)))/(jMax*(Power(a0,4) - 4*Power(a0,3)*jMax*tf + 6*Power(a0,2)*Power(jMax,2)*Power(tf,2) - 24*a0*Power(jMax,2)*(p0 - pf + tf*vf) + 12*Power(jMax,2)*(Power(v0 - vf,2) + jMax*tf*(2*p0 - 2*pf + tf*(v0 + vf))))));
                profile.t[4] = (-2*Power(a0,4)*t + Power(a0,3)*(jMax*(3*Power(t,2) + 8*t*tf - Power(tf,2)) - 2*v0 + 2*vf) + 3*Power(a0,2)*jMax*(-6*p0 + 6*pf + 2*jMax*Power(t,3) - 3*jMax*Power(t,2)*tf - 3*jMax*t*Power(tf,2) + 2*jMax*Power(tf,3) + 2*t*v0 + 2*tf*v0 - 2*t*vf - 8*tf*vf) + 6*Power(jMax,2)*(Power(jMax,2)*t*Power(t - tf,2)*Power(tf,2) + 2*(v0 - vf)*(p0 - pf + t*v0 - tf*v0 - t*vf + 2*tf*vf) + jMax*tf*(-2*p0*(t + tf) + 2*pf*(t + tf) - 3*Power(t,2)*v0 + 3*t*tf*v0 - 2*Power(tf,2)*v0 + 3*Power(t,2)*vf - 5*t*tf*vf)) - 6*a0*jMax*(Power(jMax,2)*Power(t - tf,2)*tf*(2*t + tf) - 2*Power(v0 - vf,2) - jMax*(2*p0*(t + 3*tf) - 2*pf*(t + 3*tf) + 3*Power(t,2)*v0 - 2*t*tf*v0 + Power(tf,2)*v0 - 3*Power(t,2)*vf + 4*t*tf*vf + 5*Power(tf,2)*vf)))/(Power(a0,4) - 4*Power(a0,3)*jMax*tf + 6*Power(a0,2)*Power(jMax,2)*Power(tf,2) - 24*a0*Power(jMax,2)*(p0 - pf + tf*vf) + 12*Power(jMax,2)*(Power(v0 - vf,2) + jMax*tf*(2*p0 - 2*pf + tf*(v0 + vf))));
                profile.t[5] = 0;
                profile.t[6] = (Power(a0,5) - Power(a0,4)*jMax*(3*t + 4*tf) + Power(a0,3)*jMax*(jMax*(3*Power(t,2) + 12*t*tf + 5*Power(tf,2)) - 2*v0 + 2*vf) + 3*Power(a0,2)*Power(jMax,2)*(-14*p0 + 14*pf + 2*jMax*Power(t,3) - 3*jMax*Power(t,2)*tf - 5*jMax*t*Power(tf,2) + 2*jMax*Power(tf,3) + 2*t*v0 + 2*tf*v0 - 2*t*vf - 16*tf*vf) + 6*Power(jMax,3)*(Power(jMax,2)*t*Power(t - tf,2)*Power(tf,2) - 2*(v0 - vf)*(-p0 + pf + tf*v0 - 2*tf*vf) + jMax*tf*(6*pf*t + 2*pf*tf - 2*p0*(3*t + tf) - 3*Power(t,2)*v0 + t*tf*v0 - 2*Power(tf,2)*v0 + 3*Power(t,2)*vf - 7*t*tf*vf)) - 6*a0*Power(jMax,2)*(Power(jMax,2)*Power(t - tf,2)*tf*(2*t + tf) - 4*Power(v0 - vf,2) - jMax*(6*p0*t - 6*pf*t + 10*p0*tf - 10*pf*tf + 3*Power(t,2)*v0 - 2*t*tf*v0 + 3*Power(tf,2)*v0 - 3*Power(t,2)*vf + 8*t*tf*vf + 7*Power(tf,2)*vf)))/(jMax*(Power(a0,4) - 4*Power(a0,3)*jMax*tf + 6*Power(a0,2)*Power(jMax,2)*Power(tf,2) - 24*a0*Power(jMax,2)*(p0 - pf + tf*vf) + 12*Power(jMax,2)*(Power(v0 - vf,2) + jMax*tf*(2*p0 - 2*pf + tf*(v0 + vf)))));

                // std::cout << t << std::endl;
                // std::cout << profile.t[0] << std::endl;
                // std::cout << profile.t[1] << std::endl;
                // std::cout << profile.t[2] << std::endl;
                // std::cout << profile.t[3] << std::endl;
                // std::cout << profile.t[4] << " " << std::endl;
                // std::cout << "---" << std::endl;

                profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
                if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

bool RuckigStep2::time_down_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc1_vel(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_vel(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_vel(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc1(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_none(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_none(profile, -vMax, -aMax, -jMax);
}

bool RuckigStep2::get_profile(Profile& profile, double vMax, double aMax, double jMax) {
    // Test all cases to get ones that match
    if (pf > p0) {
        if (time_up_acc0_acc1_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1_VEL;

        } else if (time_down_acc0_acc1_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1_VEL;

        } else if (time_up_acc0_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_VEL;

        } else if (time_down_acc0_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_VEL;

        } else if (time_up_acc1_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1_VEL;

        } else if (time_down_acc1_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1_VEL;

        } else if (time_up_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_VEL;

        } else if (time_down_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_VEL;

        } else if (time_up_none(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_NONE;

        } else if (time_up_acc0(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0;

        } else if (time_up_acc1(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1;

        } else if (time_up_acc0_acc1(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1;

        } else if (time_down_none(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_NONE;

        } else if (time_down_acc0(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0;

        } else if (time_down_acc1(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1;

        } else if (time_down_acc0_acc1(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1;

        } else {
            return false;
        }

    } else {
        if (time_down_acc0_acc1_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1_VEL;

        } else if (time_up_acc0_acc1_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1_VEL;

        } else if (time_down_acc0_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_VEL;

        } else if (time_up_acc0_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_VEL;

        } else if (time_down_acc1_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1_VEL;

        } else if (time_up_acc1_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1_VEL;

        } else if (time_down_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_VEL;

        } else if (time_up_vel(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_VEL;

        } else if (time_down_none(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_NONE;

        } else if (time_down_acc0(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0;

        } else if (time_down_acc1(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1;

        } else if (time_down_acc0_acc1(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1;

        } else if (time_up_none(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_NONE;

        } else if (time_up_acc0(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0;

        } else if (time_up_acc1(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1;

        } else if (time_up_acc0_acc1(profile, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1;

        } else {
            return false;
        }
    }
    return true;
}

} // namespace ruckig
