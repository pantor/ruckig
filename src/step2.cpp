#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>
#include <ruckig/wolfram.hpp>


namespace ruckig {

Step2::Step2(double tf, double p0, double v0, double a0, double pf, double vf, double af, double vMax, double aMax, double jMax): tf(tf), p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af) {
    
}

bool Step2::time_up_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Profile UDDU
    {
        double h1 = Sqrt(3)*Sqrt(Power(aMax,2)*Power(jMax,2)*(-3*Power(a0,4) - 3*Power(af,4) + 4*Power(a0,3)*aMax - 4*Power(af,3)*aMax + 6*Power(a0,2)*(Power(af,2) + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))) - 12*a0*aMax*(Power(af,2) + 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))) + 12*Power(af,2)*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)) + 24*af*aMax*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)) + 12*(Power(aMax,4) - 2*Power(aMax,3)*jMax*tf + Power(aMax,2)*Power(jMax,2)*Power(tf,2) - Power(jMax,2)*Power(v0 - vf,2) + 2*aMax*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf)))));

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = -(-3*Power(a0,2)*aMax*jMax + 3*Power(af,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 6*af*Power(aMax,2)*jMax + 18*Power(aMax,3)*jMax - 6*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,2)*v0 - 6*aMax*Power(jMax,2)*vf + h1)/(12.*Power(aMax,2)*Power(jMax,2));
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = (-6*Power(aMax,3)*jMax + h1)/(6.*Power(aMax,2)*Power(jMax,2));
        profile.t[4] = profile.t[2];
        profile.t[5] = -(3*Power(a0,2)*aMax*jMax - 3*Power(af,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 6*af*Power(aMax,2)*jMax + 18*Power(aMax,3)*jMax - 6*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,2)*v0 + 6*aMax*Power(jMax,2)*vf + h1)/(12.*Power(aMax,2)*Power(jMax,2));
        profile.t[6] = profile.t[4] + af/jMax;

        profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(tf, pf, vf, af, vMax, aMax)) {
            return true;
        }
    }

    // Profile UDUD
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (3*Power(a0,4) + 3*Power(af,4) - 4*Power(a0,3)*aMax - 8*Power(af,3)*aMax + 24*a0*Power(aMax,3) + 24*af*aMax*(Power(aMax,2) + jMax*(v0 - vf)) - 6*Power(af,2)*(Power(aMax,2) + 2*jMax*(v0 - vf)) + 6*Power(a0,2)*(Power(af,2) - 2*af*aMax - Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf) - 12*(2*Power(aMax,4) - 2*Power(aMax,3)*jMax*tf - 2*aMax*Power(jMax,2)*(p0 - pf + tf*v0) - Power(jMax,2)*Power(v0 - vf,2) + Power(aMax,2)*jMax*(-v0 + vf)))/(12.*aMax*jMax*(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf - jMax*v0 + jMax*vf)));
        profile.t[2] = profile.t[0] + a0/jMax;
        profile.t[3] = -(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + 4*Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf)/(2.*aMax*jMax);
        profile.t[4] = profile.t[2];
        profile.t[5] = (3*Power(a0,4) + 3*Power(af,4) - 8*Power(a0,3)*aMax - 4*Power(af,3)*aMax + 24*af*Power(aMax,3) - 6*Power(af,2)*(Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*(v0 - vf)) + 6*Power(a0,2)*(Power(af,2) - Power(aMax,2) - 2*jMax*v0 + 2*jMax*vf) - 12*a0*aMax*(Power(af,2) - 2*(Power(aMax,2) + jMax*v0 - jMax*vf)) - 12*(2*Power(aMax,4) - 2*Power(aMax,3)*jMax*tf - Power(jMax,2)*Power(v0 - vf,2) + Power(aMax,2)*jMax*(-v0 + vf) + 2*aMax*Power(jMax,2)*(p0 - pf + tf*vf)))/(12.*aMax*jMax*(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*af*aMax + 2*(Power(aMax,2) - aMax*jMax*tf - jMax*v0 + jMax*vf)));
        profile.t[6] = profile.t[4] - af/jMax;

        profile.set({jMax, 0, -jMax, 0, jMax, 0, -jMax});
        if (profile.check(tf, pf, vf, af, vMax, aMax)) {
            return true;
        }
    }
    
    return false;
}

bool Step2::time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Profile UDDU
    {
        std::array<double, 5> polynom;

        polynom[0] = 1.0;
        polynom[1] = (2*(2*a0 + aMax))/jMax;
        polynom[2] = (5*Power(a0,2) + Power(af,2) + 4*a0*aMax + 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf)/Power(jMax,2);
        polynom[3] = (2*a0*(Power(a0,2) + Power(af,2) + a0*aMax + 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf))/Power(jMax,3);
        polynom[4] = (3*Power(a0,4) + 3*Power(af,4) + 4*Power(a0,3)*aMax + 8*Power(af,3)*aMax + 6*Power(af,2)*(Power(aMax,2) + 2*jMax*(v0 - vf)) + 12*jMax*(-2*aMax*jMax*(p0 - pf + tf*v0) + Power(aMax,2)*(v0 - vf) + jMax*Power(v0 - vf,2)) + 24*af*aMax*jMax*(v0 - vf) + 6*Power(a0,2)*(Power(af,2) + 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf))/(12.*Power(jMax,4));

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0 || t > tf) {
                continue;
            }
            
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = -(Power(a0,2) + Power(af,2) + 2*a0*(aMax + 2*jMax*t) + 2*(af*aMax + Power(aMax,2) + aMax*jMax*(2*t - tf) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[4] = aMax/jMax;
            profile.t[5] = (Power(a0,2) + Power(af,2) - 2*Power(aMax,2) + 4*a0*jMax*t + 2*jMax*(jMax*Power(t,2) + v0 - vf))/(2.*aMax*jMax);
            profile.t[6] = profile.t[4] + af/jMax;
            
            profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
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

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0 || t > tf) {
                continue;
            }

            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + t;
            profile.t[3] = (Power(a0,2) - Power(af,2) - 2*a0*aMax + 2*af*aMax - 2*Power(aMax,2) + 4*a0*jMax*t - 4*aMax*jMax*t + 2*Power(jMax,2)*Power(t,2) + 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf)/(2.*aMax*jMax);
            profile.t[4] = aMax/jMax;
            profile.t[5] = -(Power(a0,2) - Power(af,2) + 4*a0*jMax*t + 2*(Power(aMax,2) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[6] = profile.t[4] - af/jMax;
            
            profile.set({jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool Step2::time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {   
    // Profile UDDU
    {
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = (2*aMax)/jMax;
        polynom[2] = (Power(a0,2) - Power(af,2) - 2*a0*aMax + 2*af*aMax + Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf)/Power(jMax,2);
        polynom[3] = 0;
        polynom[4] = -(-3*Power(a0,4) - 3*Power(af,4) + 8*Power(a0,3)*aMax + 4*Power(af,3)*aMax - 12*a0*aMax*(Power(af,2) + 2*jMax*(v0 - vf)) + 6*Power(a0,2)*(Power(af,2) - Power(aMax,2) + 2*jMax*v0 - 2*jMax*vf) + 6*Power(af,2)*(Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(-v0 + vf)) + 12*jMax*(Power(aMax,2)*(v0 - vf) - jMax*Power(v0 - vf,2) + 2*aMax*jMax*(p0 - pf + tf*vf)))/(12.*Power(jMax,4));

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0 || t > tf) {
                continue;
            }

            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (Power(a0,2) - Power(af,2) - 2*Power(aMax,2) + 2*jMax*(jMax*Power(t,2) - v0 + vf))/(2.*aMax*jMax);
            profile.t[2] = aMax/jMax;
            profile.t[3] = (-Power(a0,2) + Power(af,2) + 2*a0*aMax - 2*(af*aMax + Power(aMax,2) + aMax*jMax*(2*t - tf) + jMax*(jMax*Power(t,2) - v0 + vf)))/(2.*aMax*jMax);
            profile.t[4] = t;
            profile.t[5] = 0;
            profile.t[6] = af/jMax + t;
            
            profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
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

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0 || t > tf) {
                continue;
            }

            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (Power(a0,2) + Power(af,2) - 2*(Power(aMax,2) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[2] = aMax/jMax;
            profile.t[3] = -(Power(a0,2) + Power(af,2) - 2*a0*aMax - 2*(af*aMax - Power(aMax,2) + aMax*jMax*(-2*t + tf) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[4] = t;
            profile.t[5] = 0;
            profile.t[6] = -(af/jMax) + t;
            
            profile.set({jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool Step2::time_up_vel(Profile& profile, double vMax, double aMax, double jMax) {
    // Profile UDDU
    {
        // Find root of 5th order polynom
        std::array<double, 6> polynom;
        polynom[0] = 1.0;
        polynom[1] = (15*Power(a0,2) + Power(af,2) + 4*af*jMax*tf - 16*a0*(af - jMax*tf) - 2*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf))/(4.*jMax*(a0 - af + jMax*tf));
        polynom[2] = (29*Power(a0,3) - 2*Power(af,3) - 33*Power(a0,2)*(af - jMax*tf) + 6*Power(jMax,2)*(p0 - pf + tf*v0) + 6*af*jMax*(-v0 + vf) + 6*a0*(Power(af,2) + 4*af*jMax*tf - 2*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf)))/(6.*Power(jMax,2)*(a0 - af + jMax*tf));
        polynom[3] = (61*Power(a0,4) + Power(af,4) + 8*Power(af,3)*jMax*tf - 76*Power(a0,3)*(af - jMax*tf) - 24*Power(jMax,3)*tf*(p0 - pf + tf*v0) - 16*a0*(Power(af,3) - 3*Power(jMax,2)*(p0 - pf + tf*v0) + 3*af*jMax*(v0 - vf)) + 12*Power(af,2)*jMax*(v0 - vf) + 36*Power(jMax,2)*Power(v0 - vf,2) + 24*af*Power(jMax,2)*(p0 - pf + 2*tf*v0 - tf*vf) + 30*Power(a0,2)*(Power(af,2) + 4*af*jMax*tf - 2*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf)))/(24.*Power(jMax,3)*(a0 - af + jMax*tf));
        polynom[4] = (a0*(7*Power(a0,4) + Power(af,4) + 8*Power(af,3)*jMax*tf - 10*Power(a0,3)*(af - jMax*tf) - 24*Power(jMax,3)*tf*(p0 - pf + tf*v0) - 4*a0*(Power(af,3) - 3*Power(jMax,2)*(p0 - pf + tf*v0) + 3*af*jMax*(v0 - vf)) + 12*Power(af,2)*jMax*(v0 - vf) + 36*Power(jMax,2)*Power(v0 - vf,2) + 24*af*Power(jMax,2)*(p0 - pf + 2*tf*v0 - tf*vf) + 6*Power(a0,2)*(Power(af,2) + 4*af*jMax*tf - 2*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf))))/(12.*Power(jMax,4)*(a0 - af + jMax*tf));
        polynom[5] = (7*Power(a0,6) + Power(af,6) - 12*Power(a0,5)*(af - jMax*tf) + 48*Power(af,3)*Power(jMax,2)*(p0 - pf + tf*v0) - 8*Power(a0,3)*(Power(af,3) - 3*Power(jMax,2)*(p0 - pf + tf*v0) + 3*af*jMax*(v0 - vf)) - 72*Power(jMax,3)*(jMax*Power(p0 - pf + tf*v0,2) - Power(v0 - vf,3)) + 6*Power(af,4)*jMax*(v0 - vf) + 144*af*Power(jMax,3)*(p0 - pf + tf*v0)*(v0 - vf) + 36*Power(af,2)*Power(jMax,2)*Power(v0 - vf,2) + 9*Power(a0,4)*(Power(af,2) + 4*af*jMax*tf - 2*jMax*(jMax*Power(tf,2) - 3*v0 + 3*vf)) + 3*Power(a0,2)*(Power(af,4) + 8*Power(af,3)*jMax*tf - 24*Power(jMax,3)*tf*(p0 - pf + tf*v0) + 12*Power(af,2)*jMax*(v0 - vf) + 36*Power(jMax,2)*Power(v0 - vf,2) + 24*af*Power(jMax,2)*(p0 - pf + 2*tf*v0 - tf*vf)))/(144.*Power(jMax,5)*(a0 - af + jMax*tf));

        std::array<double, 5> deriv;
        deriv[0] = 1;
        deriv[1] = 4./5 * polynom[1];
        deriv[2] = 3./5 * polynom[2];
        deriv[3] = 2./5 * polynom[3];
        deriv[4] = 1./5 * polynom[4];

        // Solve 4th order derivative analytically
        auto extremas = Roots::solveQuartMonic(deriv[1], deriv[2], deriv[3], deriv[4]);
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
            // std::cout << "tz: " << tz << " " << val_new << " " << val_current << std::endl;
            if (std::abs(val_new) < 1e-15) {
                tz_intervals.insert({tz - 2e-16, tz + 2e-16});
                tz += 1e-15;
                
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
            double t = Roots::shrinkInterval(polynom, lower, upper, 1e-16);

            // Single Newton step
            {
                double h2 = Sqrt((Power(a0,2) + Power(af,2) + 4*a0*jMax*t + 2*jMax*(jMax*Power(t,2) + v0 - vf))/Power(jMax,2));
                double orig = -pf -(2*Power(a0,3) + 4*Power(af,3) + 12*a0*jMax*t*(2*af + 2*jMax*t - 2*jMax*tf + Sqrt(2)*jMax*h2) + 3*Power(a0,2)*(2*af + 4*jMax*t - 2*jMax*tf + Sqrt(2)*jMax*h2) + 3*Sqrt(2)*Power(af,2)*jMax*h2 + 12*af*jMax*(jMax*Power(t,2) + v0 - vf) + 6*Power(jMax,2)*(-2*p0 + 2*jMax*Power(t,3) - 2*jMax*Power(t,2)*tf - 2*tf*v0 + Sqrt(2)*jMax*Power(t,2)*h2 + Sqrt(2)*v0*h2 - Sqrt(2)*h2*vf))/(12.*Power(jMax,2));
                double deriv = -((a0 + jMax*t)*(3*Sqrt(2)*(Power(a0,2) + Power(af,2)) + 12*Sqrt(2)*a0*jMax*t + 6*Sqrt(2)*Power(jMax,2)*Power(t,2) + 6*Sqrt(2)*jMax*(v0 - vf) + 2*a0*jMax*h2 + 4*af*jMax*h2 + 6*Power(jMax,2)*t*h2 - 4*Power(jMax,2)*tf*h2))/(2*Power(jMax,2)*h2);
            
                // std::cout << "t: " << t << " " << orig << " new t: " << t - orig / deriv << std::endl;
                t = t - orig / deriv;
            }

            if (t < 0.0 || t > tf) {
                continue;
            }
            
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[4] = Sqrt(Power(a0,2)/2 + Power(af,2)/2 + jMax*(2*a0*t + jMax*Power(t,2) + v0 - vf))/Abs(jMax);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4] + af/jMax;
            profile.t[3] = tf - (profile.t[0] + profile.t[2] + profile.t[4] + profile.t[6]);

            profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        // Find root of 6th order polynom
        std::array<double, 7> polynom;
        polynom[0] = 1.0;
        polynom[1] = -((-5*a0 + af + jMax*tf)/jMax);
        polynom[2] = (39*Power(a0,2) - Power(af,2) + 4*af*jMax*tf - 16*a0*(af + jMax*tf) + 2*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf))/(4.*Power(jMax,2));
        polynom[3] = (55*Power(a0,3) - 33*Power(a0,2)*(af + jMax*tf) - 6*a0*(Power(af,2) - 4*af*jMax*tf - 2*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf)) + 2*(Power(af,3) - 3*Power(jMax,2)*(p0 - pf + tf*v0) + 3*af*jMax*(-v0 + vf)))/(6.*Power(jMax,3));
        polynom[4] = (101*Power(a0,4) + Power(af,4) - 8*Power(af,3)*jMax*tf - 76*Power(a0,3)*(af + jMax*tf) - 30*Power(a0,2)*(Power(af,2) - 4*af*jMax*tf - 2*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf)) + 12*Power(jMax,2)*(2*jMax*tf*(p0 - pf + tf*v0) + 3*Power(v0 - vf,2)) - 12*Power(af,2)*jMax*(v0 - vf) + 24*af*Power(jMax,2)*(p0 - pf + 2*tf*v0 - tf*vf) + 16*a0*(Power(af,3) - 3*Power(jMax,2)*(p0 - pf + tf*v0) + 3*af*jMax*(-v0 + vf)))/(24.*Power(jMax,4));
        polynom[5] = (a0*(11*Power(a0,4) + Power(af,4) - 8*Power(af,3)*jMax*tf - 10*Power(a0,3)*(af + jMax*tf) - 6*Power(a0,2)*(Power(af,2) - 4*af*jMax*tf - 2*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf)) + 12*Power(jMax,2)*(2*jMax*tf*(p0 - pf + tf*v0) + 3*Power(v0 - vf,2)) - 12*Power(af,2)*jMax*(v0 - vf) + 24*af*Power(jMax,2)*(p0 - pf + 2*tf*v0 - tf*vf) + 4*a0*(Power(af,3) - 3*Power(jMax,2)*(p0 - pf + tf*v0) + 3*af*jMax*(-v0 + vf))))/(12.*Power(jMax,5));
        polynom[6] = (11*Power(a0,6) - Power(af,6) - 12*Power(a0,5)*(af + jMax*tf) - 48*Power(af,3)*Power(jMax,2)*(p0 - pf + tf*v0) - 9*Power(a0,4)*(Power(af,2) - 4*af*jMax*tf - 2*jMax*(jMax*Power(tf,2) + 3*v0 - 3*vf)) + 72*Power(jMax,3)*(jMax*Power(p0 - pf + tf*v0,2) + Power(v0 - vf,3)) + 6*Power(af,4)*jMax*(v0 - vf) + 144*af*Power(jMax,3)*(p0 - pf + tf*v0)*(v0 - vf) - 36*Power(af,2)*Power(jMax,2)*Power(v0 - vf,2) + 8*Power(a0,3)*(Power(af,3) - 3*Power(jMax,2)*(p0 - pf + tf*v0) + 3*af*jMax*(-v0 + vf)) + 3*Power(a0,2)*(Power(af,4) - 8*Power(af,3)*jMax*tf + 12*Power(jMax,2)*(2*jMax*tf*(p0 - pf + tf*v0) + 3*Power(v0 - vf,2)) - 12*Power(af,2)*jMax*(v0 - vf) + 24*af*Power(jMax,2)*(p0 - pf + 2*tf*v0 - tf*vf)))/(144.*Power(jMax,6));

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
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[4] = Sqrt((vf - vPlat)/jMax); 
            profile.t[5] = 0;
            profile.t[6] = profile.t[4] - af/jMax;
            profile.t[3] = tf - (profile.t[0] + profile.t[2] + profile.t[4] + profile.t[6]);

            profile.set({jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool Step2::time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    if (std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        profile.t[0] = (Power(aMax,2)*Power(tf,2) - Power(v0 - vf,2) + 2*aMax*(2*p0 - 2*pf + tf*(v0 + vf)))/(2.*Power(aMax,2)*tf);
        profile.t[1] = -(Power(aMax,2)*Power(tf,2) - 2*Power(v0 - vf,2) + aMax*(8*p0 - 8*pf + 5*tf*v0 + 3*tf*vf))/(2.*Power(aMax,2)*tf);
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = -(Power(aMax,2)*Power(tf,2) - 2*Power(v0 - vf,2) + aMax*(8*p0 - 8*pf + 3*tf*v0 + 5*tf*vf))/(2.*Power(aMax,2)*tf);
        profile.t[6] = profile.t[0];
        double jMaxNew = aMax/profile.t[0];

        profile.set({jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
        return profile.check(tf, pf, vf, af, vMax, aMax) && std::abs(jMaxNew) <= std::abs(jMax) + 1e-12;
    }

    double h3 = Sqrt(aMax*(6*Power(aMax,5)*Power(tf,2) + 3*Power(af,4)*(2*p0 - 2*pf + tf*(aMax*tf + 2*v0)) + 3*Power(a0,4)*(2*p0 - 2*pf + tf*(aMax*tf + 2*vf)) + 6*Power(af,2)*aMax*(Power(aMax,2)*Power(tf,2) + 2*Power(v0 - vf,2) + aMax*(-4*p0 + 4*pf + tf*v0 - 5*tf*vf)) + 4*Power(af,3)*(2*Power(aMax,2)*Power(tf,2) + Power(v0 - vf,2) + 2*aMax*(p0 - pf + 2*tf*v0 - tf*vf)) - 4*Power(a0,3)*(2*Power(aMax,2)*Power(tf,2) + Power(v0 - vf,2) + 2*aMax*(p0 - pf - tf*v0 + 2*tf*vf)) - 12*af*Power(aMax,2)*(-Power(v0 - vf,2) + aMax*(4*p0 - 4*pf + tf*v0 + 3*tf*vf)) - 6*Power(a0,2)*(Power(af,2)*(2*p0 - 2*pf + tf*(v0 + vf)) + 2*af*aMax*(2*p0 - 2*pf + tf*(v0 + vf)) - aMax*(Power(aMax,2)*Power(tf,2) + 2*Power(v0 - vf,2) + aMax*(-4*p0 + 4*pf - 5*tf*v0 + tf*vf))) + 12*a0*aMax*(Power(af,2)*(2*p0 - 2*pf + tf*(v0 + vf)) + 2*af*aMax*(2*p0 - 2*pf + tf*(v0 + vf)) + aMax*(-Power(v0 - vf,2) + aMax*(4*p0 - 4*pf + 3*tf*v0 + tf*vf)))));
    double h1 = (3*Power(a0,2)*aMax*tf + 3*Power(af,2)*aMax*tf - 6*a0*Power(aMax,2)*tf + 6*af*Power(aMax,2)*tf + 6*Power(aMax,3)*tf - 3*Power(a0,2)*v0 + 3*Power(af,2)*v0 + 6*a0*aMax*v0 + 6*af*aMax*v0 + 3*Power(a0,2)*vf - 3*Power(af,2)*vf - 6*a0*aMax*vf - 6*af*aMax*vf - Sqrt(6)*h3);
    double h5 = (3*Power(a0,3) - 3*Power(af,3) + Power(a0,2)*(3*af - 4*aMax) - 4*Power(af,2)*aMax + 12*af*Power(aMax,2) + 24*Power(aMax,3) - a0*(3*Power(af,2) + 16*af*aMax + 12*Power(aMax,2)));
    double h7 = 2*(a0 - af)*aMax*h5;
    double jMaxNew = -(-12*Power(a0,2)*aMax*tf - 12*Power(af,2)*aMax*tf + 24*a0*Power(aMax,2)*tf - 24*af*Power(aMax,2)*tf - 24*Power(aMax,3)*tf + 12*Power(a0,2)*v0 - 12*Power(af,2)*v0 - 24*a0*aMax*v0 - 24*af*aMax*v0 - 12*Power(a0,2)*vf + 12*Power(af,2)*vf + 24*a0*aMax*vf + 24*af*aMax*vf - 4*Sqrt(9*Power(2*Power(aMax,3)*tf + Power(af,2)*(aMax*tf + v0 - vf) + 2*af*aMax*(aMax*tf + v0 - vf) + Power(a0,2)*(aMax*tf - v0 + vf) - 2*a0*aMax*(aMax*tf - v0 + vf),2) + 3*(a0 - af)*h5*(Power(aMax,2)*Power(tf,2) - Power(v0 - vf,2) + 2*aMax*(2*p0 - 2*pf + tf*(v0 + vf)))))/(24.*(Power(aMax,2)*Power(tf,2) - Power(v0 - vf,2) + 2*aMax*(2*p0 - 2*pf + tf*(v0 + vf))));

    profile.t[0] = (2*(a0 - aMax)*h1)/((a0 - af)*h5);
    profile.t[1] = (-6*Power(a0,2)*Power(af,2)*aMax*tf + 6*Power(af,4)*aMax*tf - 4*Power(a0,3)*Power(aMax,2)*tf - 12*Power(a0,2)*af*Power(aMax,2)*tf + 16*Power(af,3)*Power(aMax,2)*tf + 18*Power(a0,2)*Power(aMax,3)*tf + 30*Power(af,2)*Power(aMax,3)*tf - 36*a0*Power(aMax,4)*tf + 36*af*Power(aMax,4)*tf + 48*Power(aMax,5)*tf + 4*Power(a0,3)*aMax*v0 - 12*a0*Power(af,2)*aMax*v0 + 8*Power(af,3)*aMax*v0 - 24*Power(a0,2)*Power(aMax,2)*v0 - 24*a0*af*Power(aMax,2)*v0 + 48*Power(af,2)*Power(aMax,2)*v0 + 24*a0*Power(aMax,3)*v0 + 72*af*Power(aMax,3)*v0 - 4*Power(a0,3)*aMax*vf + 12*a0*Power(af,2)*aMax*vf - 8*Power(af,3)*aMax*vf + 24*Power(a0,2)*Power(aMax,2)*vf + 24*a0*af*Power(aMax,2)*vf - 48*Power(af,2)*Power(aMax,2)*vf - 24*a0*Power(aMax,3)*vf - 72*af*Power(aMax,3)*vf + Sqrt(6)*Power(a0,2)*h3 - Sqrt(6)*Power(af,2)*h3 + 2*Sqrt(6)*a0*aMax*h3 - 2*Sqrt(6)*af*aMax*h3 - 8*Sqrt(6)*Power(aMax,2)*h3)/h7;
    profile.t[2] = profile.t[0] + a0/jMaxNew;
    profile.t[3] = 0;
    profile.t[4] = profile.t[2];
    profile.t[5] = (6*Power(a0,4)*aMax*tf - 6*Power(a0,2)*Power(af,2)*aMax*tf - 16*Power(a0,3)*Power(aMax,2)*tf + 12*a0*Power(af,2)*Power(aMax,2)*tf + 4*Power(af,3)*Power(aMax,2)*tf + 30*Power(a0,2)*Power(aMax,3)*tf + 18*Power(af,2)*Power(aMax,3)*tf - 36*a0*Power(aMax,4)*tf + 36*af*Power(aMax,4)*tf + 48*Power(aMax,5)*tf + 8*Power(a0,3)*aMax*v0 - 12*Power(a0,2)*af*aMax*v0 + 4*Power(af,3)*aMax*v0 - 48*Power(a0,2)*Power(aMax,2)*v0 + 24*a0*af*Power(aMax,2)*v0 + 24*Power(af,2)*Power(aMax,2)*v0 + 72*a0*Power(aMax,3)*v0 + 24*af*Power(aMax,3)*v0 - 8*Power(a0,3)*aMax*vf + 12*Power(a0,2)*af*aMax*vf - 4*Power(af,3)*aMax*vf + 48*Power(a0,2)*Power(aMax,2)*vf - 24*a0*af*Power(aMax,2)*vf - 24*Power(af,2)*Power(aMax,2)*vf - 72*a0*Power(aMax,3)*vf - 24*af*Power(aMax,3)*vf - Sqrt(6)*Power(a0,2)*h3 + Sqrt(6)*Power(af,2)*h3 + 2*Sqrt(6)*a0*aMax*h3 - 2*Sqrt(6)*af*aMax*h3 - 8*Sqrt(6)*Power(aMax,2)*h3)/h7;
    profile.t[6] = profile.t[4] + af/jMaxNew;
    
    profile.set({jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
    return profile.check(tf, pf, vf, af, vMax, aMax) && std::abs(jMaxNew) <= std::abs(jMax) + 1e-12;
}

bool Step2::time_up_acc1(Profile& profile, double vMax, double aMax, double jMax) {
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

        profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
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

        profile.set({jMax, 0, -jMax, 0, jMax, 0, -jMax});
        if (profile.check(tf, pf, vf, af, vMax, aMax)) {
            return true;
        }
    }
    return false;
}

bool Step2::time_up_acc0(Profile& profile, double vMax, double aMax, double jMax) {
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

    profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
    if (profile.check(tf, pf, vf, af, vMax, aMax)) {
        return true;
    }

    return false;
}

bool Step2::time_up_none(Profile& profile, double vMax, double aMax, double jMax) {
    if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        double jMaxNew = -32*(p0 - pf)/Power(tf,3);

        profile.t[0] = tf/4;
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        profile.set({jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
        return profile.check(tf, pf, vf, af, vMax, aMax) && std::abs(jMaxNew) <= std::abs(jMax) + 1e-12;
    }
    
    if (std::abs(v0) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON) {
        double h1 = Sqrt(Power(tf,2)*Power(vf,2) + Power(4*(p0 - pf) + tf*vf,2));
        double jMaxNew = 4*(-4*(p0 - pf) - 2*tf*vf + h1)/Power(tf,3);
        
        profile.t[0] = (4*(p0 - pf) + 3*tf*vf + h1)/(4*vf);
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        profile.set({jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
        if (profile.check(tf, pf, vf, af, vMax, aMax) && std::abs(jMaxNew) <= std::abs(jMax) + 1e-12) {
            return true;
        }
    }

    /* if (std::abs(a0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON) {
        // Solution 1
        {
            double jMaxNew = (-4*(4*p0*tf - 4*pf*tf + 2*Power(tf,2)*v0 + Sqrt(Power(tf,2)*(Power(tf,2)*Power(v0,2) + 4*Power(2*p0 - 2*pf + tf*v0,2)))))/Power(tf,4);

            profile.t[0] = (-4*p0 + 4*pf - tf*v0 + Sqrt(Power(tf,2)*(Power(tf,2)*Power(v0,2) + 4*Power(2*p0 - 2*pf + tf*v0,2)))/tf)/(4.*v0);
            profile.t[1] = 0;
            profile.t[2] = profile.t[0];
            profile.t[3] = 0;
            profile.t[4] = (4*p0 - 4*pf + 3*tf*v0 - Sqrt(Power(tf,2)*(Power(tf,2)*Power(v0,2) + 4*Power(2*p0 - 2*pf + tf*v0,2)))/tf)/(4.*v0);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            profile.set({jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
            if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                return true;
            }
        }
    } */

    if (std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        // Solution 1
        {
            double h1 = Sqrt(Power(tf,2)*Power(v0 - vf,2) + 4*Power(2*p0 - 2*pf + tf*(v0 + vf),2));
            double jMaxNew = (-4*(4*p0*tf - 4*pf*tf + 2*Power(tf,2)*v0 + 2*Power(tf,2)*vf + Sqrt(Power(tf,2)*(16*Power(p0,2) + 16*Power(pf,2) - 16*pf*tf*(v0 + vf) + Power(tf,2)*(5*Power(v0,2) + 6*v0*vf + 5*Power(vf,2)) + 16*p0*(-2*pf + tf*(v0 + vf))))))/Power(tf,4);
        
            profile.t[0] = -(4*(p0 - pf) + tf*v0 + 3*tf*vf - h1)/(4*(v0 - vf));
            profile.t[1] = 0;
            profile.t[2] = profile.t[0];
            profile.t[3] = 0;
            profile.t[4] = (4*(p0 - pf) + 3*tf*v0 + tf*vf - h1)/(4*(v0 - vf));
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];
            
            profile.set({jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
            if (profile.check(tf, pf, vf, af, vMax, aMax) && std::abs(jMaxNew) <= std::abs(jMax) + 1e-12) {
                return true;
            }
        }
    }

    /* std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (4*(8*p0 - 8*pf + tf*(a0*tf - af*tf + 4*(v0 + vf))))/Power(tf,3);
    polynom[2] = (-2*(Power(a0,2)*Power(tf,2) + Power(af,2)*Power(tf,2) + 16*af*(p0 - pf + tf*v0) + 8*Power(v0 - vf,2) - 2*a0*(8*p0 - 8*pf - 3*af*Power(tf,2) + 8*tf*vf)))/Power(tf,4);
    polynom[3] = (-4*Power(a0 - af,3))/(3.*Power(tf,3));
    polynom[4] = -Power(a0 - af,4)/(3.*Power(tf,4));
    auto roots = Roots::solveQuartMonic(polynom);

    for (double t: roots) {
        if (t < 0.0 || t > tf) {
            continue;
        }

        profile.t[0] = (19*Power(a0,4)*Power(tf,2) - 41*Power(af,4)*Power(tf,2) + Power(a0,3)*tf*(-64*af*tf + 51*t*Power(tf,2) + 12*v0 - 12*vf) + 3*Power(af,2)*(192*p0*t*tf - 192*pf*t*tf + 21*Power(t,2)*Power(tf,4) + 200*t*Power(tf,2)*v0 - 16*Power(v0,2) - 8*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) - 36*t*(v0 - vf)*(32*p0*t*tf - 32*pf*t*tf + Power(t,2)*Power(tf,4) + 16*t*Power(tf,2)*v0 - 16*Power(v0,2) + 16*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) - 3*Power(a0,2)*(192*p0*t*tf - 192*pf*t*tf - 10*Power(af,2)*Power(tf,2) - 85*af*t*Power(tf,3) + 27*Power(t,2)*Power(tf,4) + 44*af*tf*v0 - 40*t*Power(tf,2)*v0 + 16*Power(v0,2) - 44*af*tf*vf + 232*t*Power(tf,2)*vf - 32*v0*vf + 16*Power(vf,2)) + 3*Power(af,3)*tf*(7*t*Power(tf,2) + 36*(-v0 + vf)) - 3*af*t*(5*Power(t,2)*Power(tf,5) + 32*t*Power(tf,3)*(v0 + 4*vf) - 16*tf*(29*Power(v0,2) - 34*v0*vf + 5*Power(vf,2)) + 32*p0*(5*t*Power(tf,2) + 12*(-v0 + vf)) - 32*pf*(5*t*Power(tf,2) + 12*(-v0 + vf))) + a0*(56*Power(af,3)*Power(tf,2) + 3*Power(af,2)*tf*(83*t*Power(tf,2) + 76*(v0 - vf)) + 6*af*(3*Power(t,2)*Power(tf,4) + 168*t*Power(tf,2)*(v0 - vf) + 16*Power(v0 - vf,2)) - 3*t*(7*Power(t,2)*Power(tf,5) + 32*p0*(7*t*Power(tf,2) + 12*(v0 - vf)) - 32*pf*(7*t*Power(tf,2) + 12*(v0 - vf)) + 32*t*Power(tf,3)*(5*v0 + 2*vf) - 16*tf*(7*Power(v0,2) - 38*v0*vf + 31*Power(vf,2)))))/(24.*Power(a0 - af,3)*(a0*tf + af*tf + 2*v0 - 2*vf));
        profile.t[1] = 0;
        profile.t[2] = (-77*Power(a0,5)*Power(tf,2) + Power(a0,4)*tf*(109*af*tf - 93*t*Power(tf,2) - 180*v0 + 180*vf) + Power(a0,3)*(1728*p0*t*tf - 1728*pf*t*tf + 94*Power(af,2)*Power(tf,2) - 804*af*t*Power(tf,3) + 207*Power(t,2)*Power(tf,4) + 432*af*tf*v0 - 168*t*Power(tf,2)*v0 - 48*Power(v0,2) - 432*af*tf*vf + 1896*t*Power(tf,2)*vf + 96*v0*vf - 48*Power(vf,2)) - a0*(Power(af,4)*Power(tf,2) + 12*Power(af,3)*tf*(31*t*Power(tf,2) + 12*(v0 - vf)) - 108*t*(v0 - vf)*(32*p0*t*tf - 32*pf*t*tf + Power(t,2)*Power(tf,4) + 16*t*Power(tf,2)*v0 - 16*Power(v0,2) + 16*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) + 3*Power(af,2)*(576*p0*t*tf - 576*pf*t*tf + 81*Power(t,2)*Power(tf,4) + 1000*t*Power(tf,2)*v0 + 48*Power(v0,2) - 424*t*Power(tf,2)*vf - 96*v0*vf + 48*Power(vf,2)) - 6*af*t*(416*p0*t*Power(tf,2) - 416*pf*t*Power(tf,2) + 13*Power(t,2)*Power(tf,5) + 384*pf*(v0 - vf) + 384*p0*(-v0 + vf) + 32*t*Power(tf,3)*(5*v0 + 8*vf) - 16*tf*(49*Power(v0,2) - 74*v0*vf + 25*Power(vf,2)))) + Power(a0,2)*(-166*Power(af,3)*Power(tf,2) - 6*Power(af,2)*tf*(169*t*Power(tf,2) + 36*(v0 - vf)) + 3*af*(192*p0*t*tf - 192*pf*t*tf + 33*Power(t,2)*Power(tf,4) - 1048*t*Power(tf,2)*v0 + 48*Power(v0,2) + 1240*t*Power(tf,2)*vf - 96*v0*vf + 48*Power(vf,2)) + 3*t*(17*Power(t,2)*Power(tf,5) + 32*p0*(17*t*Power(tf,2) + 36*(v0 - vf)) - 32*pf*(17*t*Power(tf,2) + 36*(v0 - vf)) + 32*t*Power(tf,3)*(13*v0 + 4*vf) - 16*tf*(17*Power(v0,2) - 106*v0*vf + 89*Power(vf,2)))) + af*(41*Power(af,4)*Power(tf,2) + 36*t*(v0 - vf)*(32*p0*t*tf - 32*pf*t*tf + Power(t,2)*Power(tf,4) + 16*t*Power(tf,2)*v0 - 16*Power(v0,2) + 16*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) + Power(af,2)*(-576*p0*t*tf + 576*pf*t*tf - 63*Power(t,2)*Power(tf,4) - 600*t*Power(tf,2)*v0 + 48*Power(v0,2) + 24*t*Power(tf,2)*vf - 96*v0*vf + 48*Power(vf,2)) - 3*Power(af,3)*tf*(7*t*Power(tf,2) + 36*(-v0 + vf)) + 3*af*t*(5*Power(t,2)*Power(tf,5) + 32*t*Power(tf,3)*(v0 + 4*vf) - 16*tf*(29*Power(v0,2) - 34*v0*vf + 5*Power(vf,2)) + 32*p0*(5*t*Power(tf,2) + 12*(-v0 + vf)) - 32*pf*(5*t*Power(tf,2) + 12*(-v0 + vf)))))/(24.*Power(a0 - af,4)*(a0*tf + af*tf + 2*v0 - 2*vf));
        profile.t[3] = 0;
        profile.t[4] = (41*Power(a0,5)*Power(tf,2) + Power(a0,4)*tf*(-(af*tf) + 21*t*Power(tf,2) + 108*(v0 - vf)) + Power(a0,3)*(-576*p0*t*tf + 576*pf*t*tf - 166*Power(af,2)*Power(tf,2) + 372*af*t*Power(tf,3) - 63*Power(t,2)*Power(tf,4) - 144*af*tf*v0 + 24*t*Power(tf,2)*v0 + 48*Power(v0,2) + 144*af*tf*vf - 600*t*Power(tf,2)*vf - 96*v0*vf + 48*Power(vf,2)) + Power(a0,2)*(94*Power(af,3)*Power(tf,2) - 3*af*(576*p0*t*tf - 576*pf*t*tf + 81*Power(t,2)*Power(tf,4) - 424*t*Power(tf,2)*v0 + 48*Power(v0,2) + 1000*t*Power(tf,2)*vf - 96*v0*vf + 48*Power(vf,2)) + 6*Power(af,2)*tf*(169*t*Power(tf,2) + 36*(-v0 + vf)) - 3*t*(5*Power(t,2)*Power(tf,5) + 32*p0*(5*t*Power(tf,2) + 12*(v0 - vf)) - 32*pf*(5*t*Power(tf,2) + 12*(v0 - vf)) + 32*t*Power(tf,3)*(4*v0 + vf) - 16*tf*(5*Power(v0,2) - 34*v0*vf + 29*Power(vf,2)))) + a0*(109*Power(af,4)*Power(tf,2) + 12*Power(af,3)*tf*(67*t*Power(tf,2) + 36*(v0 - vf)) - 36*t*(v0 - vf)*(32*p0*t*tf - 32*pf*t*tf + Power(t,2)*Power(tf,4) + 16*t*Power(tf,2)*v0 - 16*Power(v0,2) + 16*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) + 3*Power(af,2)*(192*p0*t*tf - 192*pf*t*tf + 33*Power(t,2)*Power(tf,4) + 1240*t*Power(tf,2)*v0 + 48*Power(v0,2) - 1048*t*Power(tf,2)*vf - 96*v0*vf + 48*Power(vf,2)) - 6*af*t*(416*p0*t*Power(tf,2) - 416*pf*t*Power(tf,2) + 13*Power(t,2)*Power(tf,5) + 384*p0*(v0 - vf) - 384*pf*(v0 - vf) + 32*t*Power(tf,3)*(8*v0 + 5*vf) - 16*tf*(25*Power(v0,2) - 74*v0*vf + 49*Power(vf,2)))) - af*(77*Power(af,4)*Power(tf,2) - 3*Power(af,2)*(576*p0*t*tf - 576*pf*t*tf + 69*Power(t,2)*Power(tf,4) + 632*t*Power(tf,2)*v0 - 16*Power(v0,2) - 56*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) + 108*t*(v0 - vf)*(32*p0*t*tf - 32*pf*t*tf + Power(t,2)*Power(tf,4) + 16*t*Power(tf,2)*v0 - 16*Power(v0,2) + 16*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) - 3*Power(af,3)*tf*(31*t*Power(tf,2) + 60*(-v0 + vf)) + 3*af*t*(17*Power(t,2)*Power(tf,5) + 32*t*Power(tf,3)*(4*v0 + 13*vf) - 16*tf*(89*Power(v0,2) - 106*v0*vf + 17*Power(vf,2)) + 32*p0*(17*t*Power(tf,2) + 36*(-v0 + vf)) - 32*pf*(17*t*Power(tf,2) + 36*(-v0 + vf)))))/(24.*Power(a0 - af,4)*(a0*tf + af*tf + 2*v0 - 2*vf));
        profile.t[5] = 0;
        profile.t[6] = (41*Power(a0,4)*Power(tf,2) - 19*Power(af,4)*Power(tf,2) + Power(a0,3)*tf*(-56*af*tf + 21*t*Power(tf,2) + 108*(v0 - vf)) + 3*Power(af,3)*tf*(17*t*Power(tf,2) - 4*v0 + 4*vf) - 36*t*(v0 - vf)*(32*p0*t*tf - 32*pf*t*tf + Power(t,2)*Power(tf,4) + 16*t*Power(tf,2)*v0 - 16*Power(v0,2) + 16*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) - 3*Power(a0,2)*(192*p0*t*tf - 192*pf*t*tf + 10*Power(af,2)*Power(tf,2) - 83*af*t*Power(tf,3) + 21*Power(t,2)*Power(tf,4) + 76*af*tf*v0 - 8*t*Power(tf,2)*v0 - 16*Power(v0,2) - 76*af*tf*vf + 200*t*Power(tf,2)*vf + 32*v0*vf - 16*Power(vf,2)) + 3*Power(af,2)*(192*p0*t*tf - 192*pf*t*tf + 27*Power(t,2)*Power(tf,4) + 232*t*Power(tf,2)*v0 + 16*Power(v0,2) - 40*t*Power(tf,2)*vf - 32*v0*vf + 16*Power(vf,2)) - 3*af*t*(7*Power(t,2)*Power(tf,5) + 32*t*Power(tf,3)*(2*v0 + 5*vf) - 16*tf*(31*Power(v0,2) - 38*v0*vf + 7*Power(vf,2)) + 32*p0*(7*t*Power(tf,2) + 12*(-v0 + vf)) - 32*pf*(7*t*Power(tf,2) + 12*(-v0 + vf))) + a0*(64*Power(af,3)*Power(tf,2) + 3*Power(af,2)*tf*(85*t*Power(tf,2) + 44*(v0 - vf)) - 6*af*(3*Power(t,2)*Power(tf,4) - 168*t*Power(tf,2)*(v0 - vf) + 16*Power(v0 - vf,2)) - 3*t*(5*Power(t,2)*Power(tf,5) + 32*p0*(5*t*Power(tf,2) + 12*(v0 - vf)) - 32*pf*(5*t*Power(tf,2) + 12*(v0 - vf)) + 32*t*Power(tf,3)*(4*v0 + vf) - 16*tf*(5*Power(v0,2) - 34*v0*vf + 29*Power(vf,2)))))/(24.*Power(a0 - af,3)*(a0*tf + af*tf + 2*v0 - 2*vf));

        std::cout << t << std::endl;
        std::cout << profile.t[0] << std::endl;
        std::cout << profile.t[1] << std::endl;
        std::cout << profile.t[2] << std::endl;
        std::cout << profile.t[3] << std::endl;
        std::cout << profile.t[4] << " " << -(5*Power(a0,3) - 9*Power(a0,2)*jMax*(-2*t + tf) + 6*a0*jMax*(3*jMax*t*(t - 2*tf) - v0 + vf) + 6*Power(jMax,2)*(-8*p0 + 8*pf + 2*jMax*Power(t,3) - 3*jMax*Power(t,2)*tf - 6*t*v0 - 3*tf*v0 + 6*t*vf - 5*tf*vf)) << std::endl;
        std::cout << "---" << std::endl;

        profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(tf, pf, vf, af, vMax, aMax)) {
            return true;
        }
    } */


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

                profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
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

                // Solution 2 with aPlat
                profile.t[0] = 0;
                profile.t[1] = 0;
                profile.t[2] = t;
                profile.t[3] = (-Power(a0,4) + 2*Power(a0,3)*jMax*(t - 2*tf) - 3*Power(a0,2)*jMax*(jMax*(2*Power(t,2) - 2*t*tf + Power(tf,2)) + 2*(v0 - vf)) + 6*a0*Power(jMax,2)*(4*p0 - 4*pf + tf*(jMax*(-2*Power(t,2) + t*tf + Power(tf,2)) - 2*v0 + 6*vf)) + 6*Power(jMax,2)*(Power(jMax,2)*t*(t - tf)*Power(tf,2) - 4*Power(v0 - vf,2) + jMax*(-8*p0*t + 8*pf*t + 4*p0*tf - 4*pf*tf - 4*Power(t,2)*v0 + 3*Power(tf,2)*v0 + 4*Power(t,2)*vf - 8*t*tf*vf + Power(tf,2)*vf)))/(jMax*(Power(a0,3) + 3*Power(a0,2)*jMax*tf + 6*a0*jMax*(v0 - vf) + 6*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf))));
                profile.t[4] = (2*a0*jMax - 2*Power(jMax,2)*t + Sqrt(2)*Sqrt(Power(jMax,2)*(Power(a0,2) + Power(af,2) + 2*a0*jMax*profile.t[3] - 2*jMax*(jMax*t*profile.t[3] - v0 + vf))))/(2.*Power(jMax,2));
                profile.t[5] = 0;
                profile.t[6] = (2*af*jMax + Sqrt(2)*Sqrt(Power(jMax,2)*(Power(a0,2) + Power(af,2) + 2*a0*jMax*profile.t[3] - 2*jMax*(jMax*t*profile.t[3] - v0 + vf))))/(2.*Power(jMax,2));
                
                profile.set({jMax, 0, -jMax, 0, -jMax, 0, jMax});
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
            polynom[0] = 1.0;
            polynom[1] = (4*af*tf - 2*jMax*Power(tf,2) + 4*v0 - 4*vf)/(a0 - af + jMax*tf);
            polynom[2] = (-2*Power(a0,4) - 2*Power(af,4) + 8*Power(af,3)*jMax*tf + 6*Power(af,2)*Power(jMax,2)*Power(tf,2) + 8*Power(a0,3)*(af - jMax*tf) - 12*Power(a0,2)*Power(af - jMax*tf,2) - 12*af*Power(jMax,2)*(p0 - pf + jMax*Power(tf,3) - 2*tf*v0 + 3*tf*vf) + 2*a0*(4*Power(af,3) - 12*Power(af,2)*jMax*tf + 9*af*Power(jMax,2)*Power(tf,2) - 3*Power(jMax,2)*(-2*p0 + 2*pf + jMax*Power(tf,3) - 2*tf*vf)) + 3*Power(jMax,2)*(Power(jMax,2)*Power(tf,4) + 4*Power(v0 - vf,2) - 4*jMax*tf*(-p0 + pf + tf*v0 - 2*tf*vf)))/(3.*Power(jMax,2)*Power(a0 - af + jMax*tf,2));
            polynom[3] = (-Power(a0,5) + Power(af,5) - Power(af,4)*jMax*tf + 5*Power(a0,4)*(af - jMax*tf) - 2*Power(a0,3)*(5*Power(af,2) - 8*af*jMax*tf + 2*jMax*(2*jMax*Power(tf,2) + v0 - vf)) - 4*Power(af,3)*jMax*(jMax*Power(tf,2) - v0 + vf) + 12*Power(af,2)*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf)) - 12*af*Power(jMax,2)*(-Power(v0 - vf,2) + jMax*tf*(2*p0 - 2*pf + 3*tf*v0 - tf*vf)) + 2*Power(a0,2)*(5*Power(af,3) - 9*Power(af,2)*jMax*tf + 6*af*jMax*(v0 - vf) + 6*Power(jMax,2)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf)) + 12*Power(jMax,3)*(jMax*Power(tf,2)*(p0 - pf + tf*v0) + (v0 - vf)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf)) + a0*(-5*Power(af,4) + 8*Power(af,3)*jMax*tf + 12*Power(af,2)*jMax*(jMax*Power(tf,2) - v0 + vf) - 24*af*Power(jMax,2)*(2*p0 - 2*pf + jMax*Power(tf,3) + 2*tf*vf) + 6*Power(jMax,2)*(Power(jMax,2)*Power(tf,4) - 2*Power(v0 - vf,2) + 8*jMax*tf*(p0 - pf + tf*vf))))/(3.*Power(jMax,3)*Power(a0 - af + jMax*tf,2));
            polynom[4] = -(Power(a0,6) + Power(af,6) - 6*Power(a0,5)*(af - jMax*tf) + 48*Power(af,3)*Power(jMax,2)*(p0 - pf + tf*v0) - 72*Power(jMax,3)*(jMax*Power(p0 - pf + tf*v0,2) - Power(v0 - vf,3)) + 3*Power(a0,4)*(5*Power(af,2) - 8*af*jMax*tf + 2*jMax*(2*jMax*Power(tf,2) + v0 - vf)) + 6*Power(af,4)*jMax*(v0 - vf) + 144*af*Power(jMax,3)*(p0 - pf + tf*v0)*(v0 - vf) + 36*Power(af,2)*Power(jMax,2)*Power(v0 - vf,2) - 4*Power(a0,3)*(5*Power(af,3) - 9*Power(af,2)*jMax*tf + 6*af*jMax*(v0 - vf) + 6*Power(jMax,2)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf)) + 3*Power(a0,2)*(5*Power(af,4) - 8*Power(af,3)*jMax*tf - 12*Power(af,2)*jMax*(jMax*Power(tf,2) - v0 + vf) + 24*af*Power(jMax,2)*(2*p0 - 2*pf + jMax*Power(tf,3) + 2*tf*vf) - 6*Power(jMax,2)*(Power(jMax,2)*Power(tf,4) - 2*Power(v0 - vf,2) + 8*jMax*tf*(p0 - pf + tf*vf))) - 6*a0*(Power(af,5) - Power(af,4)*jMax*tf - 4*Power(af,3)*jMax*(jMax*Power(tf,2) - v0 + vf) + 12*Power(af,2)*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf)) - 12*af*Power(jMax,2)*(-Power(v0 - vf,2) + jMax*tf*(2*p0 - 2*pf + 3*tf*v0 - tf*vf)) + 12*Power(jMax,3)*(jMax*Power(tf,2)*(p0 - pf + tf*v0) + (v0 - vf)*(2*p0 - 2*pf - tf*v0 + 3*tf*vf))))/(18.*Power(jMax,4)*Power(a0 - af + jMax*tf,2));
       
            auto roots = Roots::solveQuartMonic(polynom);

            for (double t: roots) {
                if (t < 0.0 || t > tf) {
                    continue;
                }

                profile.t[0] = t;
                profile.t[1] = -((-Power(a0,5) + Power(af,5) - Power(af,4)*jMax*(4*t + tf) + Power(a0,4)*(5*af - 4*jMax*t - 5*jMax*tf) - 2*Power(a0,3)*(5*Power(af,2) - 8*af*jMax*(t + tf) + jMax*(jMax*(t + 2*tf)*(3*t + 2*tf) + 2*(v0 - vf))) + 2*Power(af,3)*jMax*(jMax*(3*Power(t,2) + 2*t*tf + Power(tf,2)) + 2*(v0 - vf)) + 6*Power(af,2)*Power(jMax,2)*(10*p0 - 10*pf + jMax*(2*Power(t,3) - 9*Power(t,2)*tf + 7*t*Power(tf,2) - Power(tf,3)) - 2*t*v0 + 10*tf*v0 + 2*t*vf) + 12*Power(jMax,3)*(Power(jMax,2)*t*Power(t - tf,2)*Power(tf,2) + (v0 - vf)*(2*p0 - 2*pf + 2*t*v0 - 3*tf*v0 - 2*t*vf + 5*tf*vf) + jMax*tf*(2*p0*(t + 2*tf) - 2*pf*(t + 2*tf) + 3*Power(t,2)*v0 - 3*t*tf*v0 + 3*Power(tf,2)*v0 - 3*Power(t,2)*vf + 5*t*tf*vf + Power(tf,2)*vf)) + 2*Power(a0,2)*(5*Power(af,3) - 3*Power(af,2)*jMax*(4*t + 3*tf) + 3*af*jMax*(3*jMax*(Power(t,2) + 2*t*tf - Power(tf,2)) + 2*(v0 - vf)) + 3*Power(jMax,2)*(10*p0 - 10*pf + jMax*(2*Power(t,3) - 3*Power(t,2)*tf - 3*t*Power(tf,2) + Power(tf,3)) - 2*(t + tf)*v0 + 2*(t + 6*tf)*vf)) - 12*af*Power(jMax,2)*(2*Power(jMax,2)*t*(t - 2*tf)*(t - tf)*tf - 3*Power(v0 - vf,2) + jMax*(2*p0*(t + 4*tf) - 2*pf*(t + 4*tf) + 3*(Power(t,2) - 2*t*tf + 3*Power(tf,2))*v0 - (3*Power(t,2) - 8*t*tf + Power(tf,2))*vf)) + a0*(-5*Power(af,4) + 8*Power(af,3)*jMax*(2*t + tf) - 6*Power(af,2)*jMax*(jMax*(3*t - 2*tf)*(t + 2*tf) + 2*(v0 - vf)) - 24*af*Power(jMax,2)*(5*p0 - 5*pf + jMax*(Power(t,3) - 3*Power(t,2)*tf + 2*Power(tf,3)) - t*v0 + 2*tf*v0 + t*vf + 3*tf*vf) + 12*Power(jMax,2)*(Power(jMax,2)*Power(t - tf,2)*tf*(2*t + tf) - 3*Power(v0 - vf,2) + jMax*(2*p0*(t + 5*tf) - 2*pf*(t + 5*tf) + 3*Power(t,2)*v0 - 2*t*tf*v0 + 2*Power(tf,2)*v0 - 3*Power(t,2)*vf + 4*t*tf*vf + 8*Power(tf,2)*vf))))/(jMax*(Power(a0,4) + Power(af,4) - 4*Power(af,3)*jMax*tf + 6*Power(af,2)*Power(jMax,2)*Power(tf,2) - 4*Power(a0,3)*(af - jMax*tf) + 6*Power(a0,2)*Power(af - jMax*tf,2) + 24*af*Power(jMax,2)*(p0 - pf + tf*v0) - 4*a0*(Power(af,3) - 3*Power(af,2)*jMax*tf + 6*Power(jMax,2)*(p0 - pf + tf*vf)) + 12*Power(jMax,2)*(Power(v0 - vf,2) - jMax*tf*(2*p0 - 2*pf + tf*(v0 + vf))))));
                profile.t[2] = (-2*Power(a0,4)*t + 2*Power(af,4)*(-t + tf) + Power(af,3)*(3*jMax*Power(t,2) + 2*jMax*t*tf - 4*jMax*Power(tf,2) + 2*v0 - 2*vf) + Power(a0,3)*(t*(8*af - 3*jMax*t) - 2*(af + 4*jMax*t)*tf + jMax*Power(tf,2) - 2*v0 + 2*vf) + 3*Power(af,2)*jMax*(6*p0 - 6*pf + jMax*t*(2*t - 7*tf)*(t - tf) - 2*t*v0 + 6*tf*v0 + 2*t*vf) + 6*Power(jMax,2)*(Power(jMax,2)*t*Power(t - tf,2)*Power(tf,2) + 2*(v0 - vf)*(p0 - pf + t*v0 - tf*v0 - t*vf + 2*tf*vf) + jMax*tf*(2*p0*(t + tf) - 2*pf*(t + tf) + 3*Power(t,2)*v0 - 3*t*tf*v0 + 2*Power(tf,2)*v0 - 3*Power(t,2)*vf + 5*t*tf*vf)) - 6*af*jMax*(2*Power(jMax,2)*t*(t - 2*tf)*(t - tf)*tf - 2*Power(v0 - vf,2) + jMax*(2*p0*(t + 2*tf) - 2*pf*(t + 2*tf) + 3*Power(t,2)*v0 - 6*t*tf*v0 + 6*Power(tf,2)*v0 - 3*Power(t,2)*vf + 8*t*tf*vf - 2*Power(tf,2)*vf)) + 3*Power(a0,2)*(Power(jMax,2)*(t - 2*tf)*(2*t - tf)*(t + tf) + 2*af*(af*(-2*t + tf) + v0 - vf) + jMax*(6*p0 - 6*pf + 3*af*(Power(t,2) + 2*t*tf - 2*Power(tf,2)) - 2*(t + tf)*v0 + 2*(t + 4*tf)*vf)) + a0*(Power(af,3)*(8*t - 6*tf) - 3*Power(af,2)*(jMax*(t - tf)*(3*t + 7*tf) + 2*(v0 - vf)) - 12*af*jMax*(3*p0 - 3*pf + jMax*(Power(t,3) - 3*Power(t,2)*tf + 2*Power(tf,3)) - t*v0 + tf*v0 + t*vf + 2*tf*vf) + 6*jMax*(Power(jMax,2)*Power(t - tf,2)*tf*(2*t + tf) - 2*Power(v0 - vf,2) + jMax*(2*p0*(t + 3*tf) - 2*pf*(t + 3*tf) + 3*Power(t,2)*v0 - 2*t*tf*v0 + Power(tf,2)*v0 - 3*Power(t,2)*vf + 4*t*tf*vf + 5*Power(tf,2)*vf))))/(Power(a0,4) + Power(af,4) - 4*Power(af,3)*jMax*tf + 6*Power(af,2)*Power(jMax,2)*Power(tf,2) - 4*Power(a0,3)*(af - jMax*tf) + 6*Power(a0,2)*Power(af - jMax*tf,2) + 24*af*Power(jMax,2)*(p0 - pf + tf*v0) - 4*a0*(Power(af,3) - 3*Power(af,2)*jMax*tf + 6*Power(jMax,2)*(p0 - pf + tf*vf)) + 12*Power(jMax,2)*(Power(v0 - vf,2) - jMax*tf*(2*p0 - 2*pf + tf*(v0 + vf))));
                profile.t[3] = 0;
                profile.t[4] = -((Power(a0,5) - Power(af,5) + Power(af,4)*jMax*(3*t + 2*tf) + Power(a0,4)*(-5*af + 3*jMax*t + 4*jMax*tf) + Power(a0,3)*(10*Power(af,2) - 2*af*jMax*(6*t + 7*tf) + jMax*(3*jMax*Power(t,2) + 12*jMax*t*tf + 5*jMax*Power(tf,2) + 2*v0 - 2*vf)) - Power(af,3)*jMax*(3*jMax*Power(t,2) + 6*jMax*t*tf + 2*jMax*Power(tf,2) + 2*v0 - 2*vf) - 3*Power(af,2)*Power(jMax,2)*(14*p0 - 14*pf + jMax*t*(2*Power(t,2) - 9*t*tf + 5*Power(tf,2)) - 2*t*v0 + 14*tf*v0 + 2*t*vf) + 6*Power(jMax,3)*(-(Power(jMax,2)*t*Power(t - tf,2)*Power(tf,2)) + jMax*tf*(6*pf*t + 2*pf*tf - 2*p0*(3*t + tf) - 3*Power(t,2)*v0 + t*tf*v0 - 2*Power(tf,2)*v0 + t*(3*t - 7*tf)*vf) + 2*(v0 - vf)*(-p0 + pf + tf*v0 - 2*tf*vf)) + 6*af*Power(jMax,2)*(2*Power(jMax,2)*t*(t - 2*tf)*(t - tf)*tf - 4*Power(v0 - vf,2) + jMax*(6*p0*t - 6*pf*t + 8*p0*tf - 8*pf*tf + 3*Power(t,2)*v0 - 2*t*tf*v0 + 8*Power(tf,2)*v0 - 3*Power(t,2)*vf + 8*t*tf*vf)) + Power(a0,2)*(-10*Power(af,3) + 18*Power(af,2)*jMax*(t + tf) - 3*af*jMax*(3*jMax*Power(t,2) + 10*jMax*t*tf - 4*jMax*Power(tf,2) + 2*v0 - 2*vf) + 3*Power(jMax,2)*(-14*p0 + 14*pf + jMax*(-2*Power(t,3) + 3*Power(t,2)*tf + 5*t*Power(tf,2) - 2*Power(tf,3)) + 2*(t + tf)*v0 - 2*(t + 8*tf)*vf)) + a0*(5*Power(af,4) - 2*Power(af,3)*jMax*(6*t + 5*tf) + 3*Power(af,2)*jMax*(3*jMax*Power(t,2) + 8*jMax*t*tf - 5*jMax*Power(tf,2) + 2*v0 - 2*vf) + 12*af*Power(jMax,2)*(7*p0 - 7*pf + jMax*(Power(t,3) - 3*Power(t,2)*tf + 2*Power(tf,3)) - t*v0 + 3*tf*v0 + t*vf + 4*tf*vf) - 6*Power(jMax,2)*(Power(jMax,2)*Power(t - tf,2)*tf*(2*t + tf) - 4*Power(v0 - vf,2) + jMax*(6*p0*t - 6*pf*t + 10*p0*tf - 10*pf*tf + 3*Power(t,2)*v0 - 2*t*tf*v0 + 3*Power(tf,2)*v0 - 3*Power(t,2)*vf + 8*t*tf*vf + 7*Power(tf,2)*vf))))/(jMax*(Power(a0,4) + Power(af,4) - 4*Power(af,3)*jMax*tf + 6*Power(af,2)*Power(jMax,2)*Power(tf,2) - 4*Power(a0,3)*(af - jMax*tf) + 6*Power(a0,2)*Power(af - jMax*tf,2) + 24*af*Power(jMax,2)*(p0 - pf + tf*v0) - 4*a0*(Power(af,3) - 3*Power(af,2)*jMax*tf + 6*Power(jMax,2)*(p0 - pf + tf*vf)) + 12*Power(jMax,2)*(Power(v0 - vf,2) - jMax*tf*(2*p0 - 2*pf + tf*(v0 + vf))))));
                profile.t[5] = 0;
                profile.t[6] = 0;
 
                profile.set({jMax, 0, -jMax, 0, jMax, 0, -jMax});
                if (profile.check(tf, pf, vf, af, vMax, aMax)) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

bool Step2::time_down_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, -vMax, -aMax, -jMax);
}

bool Step2::time_down_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc1_vel(profile, -vMax, -aMax, -jMax);
}

bool Step2::time_down_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_vel(profile, -vMax, -aMax, -jMax);
}

bool Step2::time_down_vel(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_vel(profile, -vMax, -aMax, -jMax);
}

bool Step2::time_down_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, -vMax, -aMax, -jMax);
}

bool Step2::time_down_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc1(profile, -vMax, -aMax, -jMax);
}

bool Step2::time_down_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_acc0(profile, -vMax, -aMax, -jMax);
}

bool Step2::time_down_none(Profile& profile, double vMax, double aMax, double jMax) {
    return time_up_none(profile, -vMax, -aMax, -jMax);
}

bool Step2::get_profile(Profile& profile, double vMax, double aMax, double jMax) {
    profile.a[0] = a0;
    profile.v[0] = v0;
    profile.p[0] = p0;

    // Test all cases to get ones that match
    if (pf > p0) {
        if (time_up_acc0_acc1_vel(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc0_acc1_vel(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc0_vel(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc0_vel(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc1_vel(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc1_vel(profile, vMax, aMax, jMax)) {
        } else if (time_up_vel(profile, vMax, aMax, jMax)) {
        } else if (time_down_vel(profile, vMax, aMax, jMax)) {
        } else if (time_up_none(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc0(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc1(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc0_acc1(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc0(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc1(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc0_acc1(profile, vMax, aMax, jMax)) {
        } else if (time_down_none(profile, vMax, aMax, jMax)) {
        } else {
            return false;
        }

    } else {
        if (time_down_acc0_acc1_vel(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc0_acc1_vel(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc0_vel(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc0_vel(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc1_vel(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc1_vel(profile, vMax, aMax, jMax)) {
        } else if (time_down_vel(profile, vMax, aMax, jMax)) {
        } else if (time_up_vel(profile, vMax, aMax, jMax)) {
        } else if (time_down_none(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc0(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc1(profile, vMax, aMax, jMax)) {
        } else if (time_down_acc0_acc1(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc0(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc1(profile, vMax, aMax, jMax)) {
        } else if (time_up_acc0_acc1(profile, vMax, aMax, jMax)) {
        } else if (time_up_none(profile, vMax, aMax, jMax)) {
        } else {
            return false;
        }
    }
    return true;
}

} // namespace ruckig
