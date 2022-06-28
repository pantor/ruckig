#pragma once

#include <array>
#include <optional>


namespace ruckig {

using Limits = Profile::Limits;
using JerkSigns = Profile::JerkSigns;

//! Mathematical equations for Step 1 in velocity interface: Extremal profiles
class VelocityStep1 {
    double a0, af;
    double _aMax, _aMin, _jMax;

    // Pre-calculated expressions
    double vd;

    // Max 3 valid profiles
    std::array<Profile, 3> valid_profiles;
    size_t valid_profile_counter;

    void time_acc0(Profile& profile, double aMax, double aMin, double jMax, bool return_after_found);
    void time_none(Profile& profile, double aMax, double aMin, double jMax, bool return_after_found);

    inline void add_profile(const Profile& profile) {
        valid_profiles[valid_profile_counter] = profile;
        valid_profiles[valid_profile_counter].pf = profile.p.back();
        ++valid_profile_counter;
    }

public:
    explicit VelocityStep1(double v0, double a0, double vf, double af, double aMax, double aMin, double jMax);

    bool get_profile(const Profile& input, Block& block);
};


//! Mathematical equations for Step 2 in velocity interface: Time synchronization
class VelocityStep2 {
    double a0, tf, af; 
    double _aMax, _aMin, _jMax;

    // Pre-calculated expressions
    double vd;

    bool time_acc0(Profile& profile, double aMax, double aMin, double jMax);
    bool time_none(Profile& profile, double aMax, double aMin, double jMax);

    inline bool check_all(Profile& profile, double aMax, double aMin, double jMax) {
        return time_acc0(profile, aMax, aMin, jMax) || time_none(profile, aMax, aMin, jMax);
    }

public:
    explicit VelocityStep2(double tf, double v0, double a0, double vf, double af, double aMax, double aMin, double jMax);

    bool get_profile(Profile& profile);
};

} // namespace ruckig
