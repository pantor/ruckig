#pragma once

#include <array>
#include <optional>


namespace ruckig {

using Limits = Profile::Limits;
using JerkSigns = Profile::JerkSigns;

//! Mathematical equations for Step 1 in velocity interface: Extremal profiles
class VelocityStep1 {
    double p0, v0, a0;
    double vf, af;
    double _aMax, _aMin, _jMax;

    // Max 3 valid profiles
    std::array<Profile, 3> valid_profiles;
    size_t valid_profile_counter;

    void time_acc0(Profile& profile, double aMax, double aMin, double jMax);
    void time_none(Profile& profile, double aMax, double aMin, double jMax);

    inline void add_profile(Profile profile, double jMax) {
        profile.pf = profile.p[7];
        profile.direction = (jMax > 0) ? Profile::Direction::UP : Profile::Direction::DOWN;
        valid_profiles[valid_profile_counter] = profile;
        ++valid_profile_counter;
    }

public:
    explicit VelocityStep1(double p0, double v0, double a0, double vf, double af, double aMax, double aMin, double jMax);

    bool get_profile(const Profile& input, Block& block);
};


//! Mathematical equations for Step 2 in velocity interface: Time synchronization
class VelocityStep2 {
    double p0, v0, a0;
    double tf, vf, af; 
    double _aMax, _aMin, _jMax;

    bool time_acc0(Profile& profile, double aMax, double aMin, double jMax);
    bool time_none(Profile& profile, double aMax, double aMin, double jMax);

    inline bool check_all(Profile& profile, double aMax, double aMin, double jMax) {
        return time_acc0(profile, aMax, aMin, jMax) || time_none(profile, aMax, aMin, jMax);
    }

public:
    explicit VelocityStep2(double tf, double p0, double v0, double a0, double vf, double af, double aMax, double aMin, double jMax);

    bool get_profile(Profile& profile);
};

} // namespace ruckig
