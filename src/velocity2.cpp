#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>


namespace ruckig {

Velocity2::Velocity2(double tf, double p0, double v0, double a0, double vf, double af, double aMax, double aMin, double jMax): p0(p0), v0(v0), a0(a0), tf(tf), vf(vf), af(af), aMax(aMax), aMin(aMin), jMax(jMax)  { }

bool Velocity2::time_acc0(Profile& profile, double aMax, double aMin, double jMax) {
    // UD Solution 1/2
    {
        const double h1 = Sqrt(-a0*a0 - af*af + 2*(a0 + af)*jMax*tf + 2*a0*af + jMax*(jMax*tf*tf - 4*(vf - v0)))/Abs(jMax);
        
        profile.t[0] = (af - a0 + jMax*tf - jMax*h1)/(2*jMax);
        profile.t[1] = h1;
        profile.t[2] = tf - (profile.t[0] + h1);
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::ACC0>(jMax, aMax, aMin)) {
            profile.pf = profile.p[7];
            return true;
        }
    }

    // UU Solution
    {
        const double h1 = (a0 - af + jMax*tf);

        profile.t[0] = -((a0*a0 + af*af - 2*a0*af - 2*jMax*(vf - v0 - a0*tf))/(2*jMax*h1));
        profile.t[1] = h1/jMax;
        profile.t[2] = 0;
        profile.t[3] = 0;
        profile.t[4] = tf - (profile.t[0] + profile.t[1]);
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDUD, Limits::ACC0>(jMax, aMax, aMin)) {
            profile.pf = profile.p[7];
            return true;
        }
    }

    return false;
}

bool Velocity2::time_none([[maybe_unused]] Profile& profile, [[maybe_unused]] double aMax, [[maybe_unused]] double aMin, [[maybe_unused]] double jMax) {
    if (std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON && std::abs(vf - v0) < DBL_EPSILON) {
        profile.t[0] = 0;
        profile.t[1] = 0;
        profile.t[2] = 0;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, aMax, aMin)) {
            profile.pf = profile.p[7];
            return true;
        }
    }
    
    // UD Solution 1/2
    {        
        profile.t[0] = -(2*(af*tf + v0 - vf))/(a0 - af);
        profile.t[1] = tf - profile.t[0];
        profile.t[2] = 0;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        double jf = (a0 - af)*(a0 - af)/(2*(af*tf + v0 - vf));

        if (std::abs(jf) < std::abs(jMax) + 1e-12 && profile.check<JerkSigns::UDDU, Limits::NONE>(jf, aMax, aMin)) {
            profile.pf = profile.p[7];
            return true;
        }
    }
    
    return false;
}

bool Velocity2::get_profile(Profile& profile) {
    profile.a[0] = a0;
    profile.v[0] = v0;
    profile.p[0] = p0;
    profile.af = af;
    profile.vf = vf;

    // Test all cases to get ones that match
    // However we should guess which one is correct and try them first...
    if (vf > v0) {
        return check_all(profile, aMax, aMin, jMax) || check_all(profile, aMin, aMax, -jMax);

    } else {
        return check_all(profile, aMin, aMax, -jMax) || check_all(profile, aMax, aMin, jMax);
    }
}

} // namespace ruckig
