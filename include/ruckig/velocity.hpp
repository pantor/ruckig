#pragma once

#include <array>
#include <optional>


namespace ruckig {

//! Mathematical equations for Step 1 in velocity interface: Extremal profiles
class VelocityStep1 {
    using ReachedLimits = Profile::ReachedLimits;
    using JerkSigns = Profile::JerkSigns;

    const double a0, af;
    const double _aMax, _aMin, _jMax;

    // Pre-calculated expressions
    double vd;

    // Max 3 valid profiles
    using ProfileIter = std::array<Profile, 3>::iterator;
    std::array<Profile, 3> valid_profiles;

    void time_acc0(ProfileIter& profile, double aMax, double aMin, double jMax, bool return_after_found) const;
    void time_none(ProfileIter& profile, double aMax, double aMin, double jMax, bool return_after_found) const;

    // Only for zero-limits case
    bool time_all_single_step(ProfileIter& profile, double aMax, double aMin, double jMax) const;

    inline void add_profile(ProfileIter& profile) const {
        const auto prev_profile = profile;
        ++profile;
        profile->set_boundary(*prev_profile);
    }

public:
    explicit VelocityStep1(double v0, double a0, double vf, double af, double aMax, double aMin, double jMax);

    bool get_profile(const Profile& input, Block& block);
};


//! Mathematical equations for Step 2 in velocity interface: Time synchronization
class VelocityStep2 {
    using ReachedLimits = Profile::ReachedLimits;
    using JerkSigns = Profile::JerkSigns;

    const double a0, tf, af;
    const double _aMax, _aMin, _jMax;

    // Pre-calculated expressions
    double vd, ad;

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
