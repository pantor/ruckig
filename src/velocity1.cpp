#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>


namespace ruckig {

Velocity1::Velocity1(double p0, double v0, double a0, double vf, double af, double aMax, double aMin, double jMax): p0(p0), v0(v0), a0(a0), vf(vf), af(af), aMax(aMax), aMin(aMin), jMax(jMax) { }

void Velocity1::time_acc0(Profile& profile, double aMax, double aMin, double jMax) {
    // UD
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (a0*a0 + af*af - 2*aMax*aMax + 2*jMax*(vf - v0))/(2*aMax*jMax);
        profile.t[2] = (-af + aMax)/jMax;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::ACC0>(jMax, aMax, aMin)) {
            add_profile(profile, jMax);
        }
    }

    // UU
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (a0*a0 - af*af + 2*jMax*(vf - v0))/(2*aMax*jMax);
        profile.t[2] = 0;
        profile.t[3] = 0;
        profile.t[4] = (af - aMax)/jMax;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDUD, Limits::ACC0>(jMax, aMax, aMin)) {
            add_profile(profile, jMax);
        }
    }
}

void Velocity1::time_none(Profile& profile, double aMax, double aMin, double jMax) {
    double h1 = Sqrt((a0*a0 + af*af)/2 + jMax*(vf - v0));

    // Solution 1
    {
        profile.t[0] = -(a0 + h1)/jMax;
        profile.t[1] = 0;
        profile.t[2] = -(af + h1)/jMax;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, aMax, aMin)) {
            add_profile(profile, jMax);
        }
    }

    // Solution 2
    {
        profile.t[0] = (-a0 + h1)/jMax;
        profile.t[1] = 0;
        profile.t[2] = (-af + h1)/jMax;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = 0;

        if (profile.check<JerkSigns::UDDU, Limits::NONE>(jMax, aMax, aMin)) {
            add_profile(profile, jMax);
        }
    }
}

bool Velocity1::get_profile(const Profile& input, Block& block) {
    Profile profile = input;
    profile.a[0] = a0;
    profile.v[0] = v0;
    profile.p[0] = p0;
    profile.af = af;
    profile.vf = vf;
    valid_profile_counter = 0;

    if (std::abs(v0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        time_none(profile, aMax, aMin, jMax);

    } else {
        time_none(profile, aMax, aMin, jMax);
        time_acc0(profile, aMax, aMin, jMax);
        time_none(profile, aMin, aMax, -jMax);
        time_acc0(profile, aMin, aMax, -jMax);
    }

    return Block::calculate_block(block, valid_profiles, valid_profile_counter);
}

} // namespace ruckig
