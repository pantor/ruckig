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

    if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1_VEL>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_up_acc1_vel(Profile& profile, double vMax, double aMax, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }
    
    const double h1 = Sqrt(a0_a0/2 + jMax*(vMax - v0))/Abs(jMax);

    profile.t[0] = (-a0 + h1*jMax)/jMax;
    profile.t[1] = 0;
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*af_p4 + 8*aMax*(af_p3 - a0_p3) + 24*aMax*jMax*(a0*v0 - af*vf) + 6*af_af*(aMax_aMax - 2*jMax*vf) - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(vf + vMax) + jMax*(-vf_vf + vMax*vMax) - aMax*h1*(a0_a0 - 2*jMax*(v0 + vMax))))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = aMax/jMax;
    profile.t[5] = (af_af/2 - aMax_aMax + jMax*(vMax - vf))/(aMax*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::ACC1_VEL>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_up_acc0_vel(Profile& profile, double vMax, double aMax, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }

    const double h1 = Sqrt(af_af/2 + jMax*(vMax - vf))/Abs(jMax);

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax_aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = profile.t[0] + a0/jMax;
    profile.t[3] = (3*a0_p4 + 8*(af_p3 - a0_p3)*aMax + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 12*af_af*aMax*h1*jMax - 12*jMax*(-2*aMax*jMax*pd + aMax_aMax*(v0 + vMax) + jMax*(-v0_v0 + vMax*vMax) + 2*aMax*(vf + vMax)*h1*jMax))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::ACC0_VEL>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_up_vel(Profile& profile, double vMax, double aMax, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }

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

    if (profile.check<JerkSigns::UDDU, Limits::VEL>(pf, vf, af, jMax, vMax, aMax)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_up_acc0_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    const double h1 = Sqrt((a0_p4 + af_p4)/2 + 4./3*(af_p3 - a0_p3)*aMax + 4*aMax*jMax*(a0*v0 - af*vf) + a0_a0*(aMax_aMax - 2*jMax*v0) + af_af*(aMax_aMax - 2*jMax*vf) + aMax_aMax*(aMax_aMax - 2*jMax*(v0 + vf)) + 2*jMax_jMax*(v0_v0 + vf_vf + 2*aMax*pd));
    
    if (!std::isnan(h1)) {
        const double h2 = a0_a0 - 3*aMax_aMax - 2*jMax*v0;
        const double h3 = ((af_af - a0_a0)/2 - jMax*(vf - v0))/(aMax*jMax);

        // UDDU: Solution 2
        {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (h2 - h1)/(2*aMax*jMax);
            profile.t[2] = profile.t[0] + a0/jMax;
            profile.t[3] = 0;
            profile.t[4] = profile.t[2];
            profile.t[5] = profile.t[1] + h3;
            profile.t[6] = profile.t[4] + af/jMax;

            if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1>(pf, vf, af, jMax, vMax, aMax)) {
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
            profile.t[5] = profile.t[1] + h3;
            profile.t[6] = profile.t[4] + af/jMax;

            if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, jMax);
            }
        }
    }
}

void Step1::time_up_acc1(Profile& profile, double vMax, double aMax, double jMax) {
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (2*(2*a0 + aMax))/jMax;
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
        profile.t[2] = (a0 + aMax)/jMax + t;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = ((a0_a0 + af_af)/2 - aMax_aMax + 2*a0*jMax*t + jMax_jMax*t*t - jMax*(vf - v0))/(aMax*jMax);
        profile.t[6] = (af + aMax)/jMax;
            
        if (profile.check<JerkSigns::UDDU, Limits::ACC1>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, jMax);
        }
    }
}

void Step1::time_up_acc0(Profile& profile, double vMax, double aMax, double jMax) {
    // Combined UDDU and UDUD
    // UDUD Strategy t7 == 0 => is equal to UDDU
    std::array<double, 5> polynom;
    polynom[0] = 1.0;
    polynom[1] = (-2*aMax)/jMax;
    polynom[2] = (-af_af + aMax_aMax + 2*jMax*vf)/jMax_jMax;
    polynom[3] = 0;
    polynom[4] = (-3*a0_p4 + 3*af_p4 + 8*a0_p3*aMax - 8*af_p3*aMax - 24*a0*aMax*jMax*v0 - 6*a0_a0*(aMax_aMax - 2*jMax*v0) + 24*af*aMax*jMax*vf + 6*af_af*(aMax_aMax - 2*jMax*vf) + 12*jMax*(-2*aMax*jMax*pd - aMax_aMax*(vf - v0) + jMax*(vf_vf - v0_v0)))/(12*jMax_jMax*jMax_jMax);

    auto roots = Roots::solveQuartMonic(polynom);
    for (double t: roots) {
        if (t < 0.0) {
            continue;
        }

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (a0_a0 - af_af + 2*jMax*(-2*aMax*t + jMax*t*t - v0 + vf))/(2*aMax*jMax);
        profile.t[2] = t;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = (af - aMax + jMax*t)/jMax;

        if (profile.check<JerkSigns::UDDU, Limits::ACC0>(pf, vf, af, jMax, vMax, aMax)) {
            add_profile(profile, jMax);
        }
    }
}

void Step1::time_up_none(Profile& profile, double vMax, double aMax, double jMax) {
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

            if (profile.check<JerkSigns::UDDU, Limits::NONE>(pf, vf, af, jMax, vMax, aMax)) {
                add_profile(profile, jMax);
            }
            return;
        }
    }

    // UDDU / UDUD modern
    {
        // UDUD Strategy: t7 == 0 (equals UDDU)
        std::array<double, 5> polynom;
        polynom[0] = 1.0;
        polynom[1] = 0;
        polynom[2] = (-2*(a0_a0 + af_af - 2*jMax*(v0 + vf)))/jMax_jMax;
        polynom[3] = (4*(a0_p3 - af_p3 - 3*jMax_jMax*pd - 3*a0*jMax*v0 + 3*af*jMax*vf))/(3*jMax*jMax_jMax);
        polynom[4] = -Power(a0_a0 - af_af + 2*jMax*(vf - v0),2)/(4*jMax_jMax*jMax_jMax);

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0) {
                continue;
            }

            profile.t[0] = (a0_a0 - af_af - 4*a0*jMax*t + 2*jMax*(jMax*t*t - v0 + vf))/(4*jMax_jMax*t);
            profile.t[1] = 0;
            profile.t[2] = t;
            profile.t[3] = 0;
            profile.t[4] = 0;
            profile.t[5] = 0;
            profile.t[6] = (-a0 + af + jMax*(t - profile.t[0]))/jMax;

            if (profile.check<JerkSigns::UDDU, Limits::NONE>(pf, vf, af, jMax, vMax, aMax)) {
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

bool Step1::calculate_block(Block& block) const {
    // if (valid_profile_counter > 0) {
    //     std::cout << "---\n " << valid_profile_counter << std::endl;
    //     for (size_t i = 0; i < valid_profile_counter; ++i) {
    //         std::cout << valid_profiles[i].t_sum[6] << " " << valid_profiles[i].to_string() << std::endl;
    //     }
    // }

    if (valid_profile_counter == 1) {
        block = Block(valid_profiles[0]);
        return true;
    
    } else if (valid_profile_counter % 2 == 0) {
        return false;
    }

    // Find index of fastest profile
    auto idx_min_it = std::min_element(valid_profiles.cbegin(), valid_profiles.cbegin() + valid_profile_counter, [](auto& a, auto& b) { return a.t_sum[6] + a.t_brake.value_or(0.0) < b.t_sum[6] + b.t_brake.value_or(0.0); });
    size_t idx_min = std::distance(valid_profiles.cbegin(), idx_min_it);
    
    block = Block(valid_profiles[idx_min]);

    if (valid_profile_counter == 3) {
        size_t idx_else_1 = (idx_min + 1) % 3;
        size_t idx_else_2 = (idx_min + 2) % 3;
        
        add_interval(block.a, idx_else_1, idx_else_2);
        return true;

    } else if (valid_profile_counter == 5) {
        size_t idx_else_1 = (idx_min + 1) % 5;
        size_t idx_else_2 = (idx_min + 2) % 5;
        size_t idx_else_3 = (idx_min + 3) % 5;
        size_t idx_else_4 = (idx_min + 4) % 5;

        if (valid_profiles[idx_else_1].direction == valid_profiles[idx_else_2].direction) {
            add_interval(block.a, idx_else_1, idx_else_2);
            add_interval(block.b, idx_else_3, idx_else_4);
        } else {
            add_interval(block.a, idx_else_1, idx_else_4);
            add_interval(block.b, idx_else_2, idx_else_3);
        }
        return true;
    }
    
    return false;
}

bool Step1::get_profile(const Profile& input, Block& block) {
    Profile profile = input;
    profile.a[0] = a0;
    profile.v[0] = v0;
    profile.p[0] = p0;
    valid_profile_counter = 0;

    if (std::abs(pf - p0) < DBL_EPSILON && std::abs(v0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        time_up_none(profile, vMax, aMax, jMax);
    
    } else if (pf > p0) {
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
        time_down_acc0(profile, vMin, aMax, jMax);
        time_down_acc1(profile, vMin, aMax, jMax);
        time_down_acc0_acc1(profile, vMin, aMax, jMax);
        time_up_none(profile, vMax, aMax, jMax);
        time_up_acc0(profile, vMax, aMax, jMax);
        time_up_acc1(profile, vMax, aMax, jMax);
        time_up_acc0_acc1(profile, vMax, aMax, jMax);
    }

    return calculate_block(block);
}

} // namespace ruckig
