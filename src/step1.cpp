#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>


namespace ruckig {

Step1::Step1(double p0, double v0, double a0, double pf, double vf, double af, double vMax, double vMin, double aMax, double aMin, double jMax): p0(p0), v0(v0), a0(a0), pf(pf), vf(vf), af(af), vMax(vMax), vMin(vMin), aMax(aMax), aMin(aMin), jMax(jMax) {
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

void Step1::time_acc0_acc1_vel(Profile& profile, double vMax, double aMax, double aMin, double jMax) {
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax*aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = aMax/jMax;
    profile.t[3] = (3*a0_p4*aMin - 3*af_p4*aMax + 8*(af_p3 - a0_p3)*aMax*aMin + 24*aMax*aMin*jMax*(a0*v0 - af*vf) + 6*a0_a0*aMin*(aMax*aMax - 2*jMax*v0) - 6*af_af*aMax*(aMin*aMin - 2*jMax*vf) - 12*jMax*(aMax*aMax*aMin*(v0 + vMax) + aMin*jMax*(vMax*vMax - v0_v0) - aMax*(2*aMin*jMax*pd + aMin*aMin*(vf + vMax) + jMax*(vMax*vMax - vf_vf))))/(24*aMax*aMin*jMax_jMax*vMax);
    profile.t[4] = -aMin/jMax;
    profile.t[5] = -(af_af/2 - aMin*aMin - jMax*(vf - vMax))/(aMin*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1_VEL>(pf, vf, af, jMax, vMax, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_acc1_vel(Profile& profile, double vMax, double aMax, double aMin, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }
    
    const double h1 = Sqrt(a0_a0/2 + jMax*(vMax - v0))/Abs(jMax);

    profile.t[0] = -a0/jMax + h1;
    profile.t[1] = 0;
    profile.t[2] = h1;
    profile.t[3] = -(3*af_p4 - 8*aMin*(af_p3 - a0_p3) - 24*aMin*jMax*(a0*v0 - af*vf) + 6*af_af*(aMin*aMin - 2*jMax*vf) - 12*jMax*(2*aMin*jMax*pd + aMin*aMin*(vf + vMax) + jMax*(vMax*vMax - vf_vf) + aMin*h1*(a0_a0 - 2*jMax*(v0 + vMax))))/(24*aMin*jMax_jMax*vMax);
    profile.t[4] = -aMin/jMax;
    profile.t[5] = -(af_af/2 - aMin*aMin + jMax*(vMax - vf))/(aMin*jMax);
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::ACC1_VEL>(pf, vf, af, jMax, vMax, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_acc0_vel(Profile& profile, double vMax, double aMax, double aMin, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }

    const double h1 = Sqrt(af_af/2 + jMax*(vMax - vf))/Abs(jMax);

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (a0_a0/2 - aMax*aMax - jMax*(v0 - vMax))/(aMax*jMax);
    profile.t[2] = aMax/jMax;
    profile.t[3] = (3*a0_p4 + 8*(af_p3 - a0_p3)*aMax + 24*aMax*jMax*(a0*v0 - af*vf) + 6*a0_a0*(aMax*aMax - 2*jMax*v0) - 12*jMax*(-2*aMax*jMax*pd + aMax*aMax*(v0 + vMax) + jMax*(vMax*vMax - v0_v0) + (2*(vf + vMax)*jMax - af_af)*aMax*h1))/(24*aMax*jMax_jMax*vMax);
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::ACC0_VEL>(pf, vf, af, jMax, vMax, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_vel(Profile& profile, double vMax, double aMax, double aMin, double jMax) {
    if ((jMax > 0 && has_up_vel) || (jMax < 0 && has_down_vel)) {
        return;
    }

    const double h1 = Sqrt(af_af/2 + jMax*(vMax - vf))/Abs(jMax);
    const double h2 = Sqrt(a0_a0/2 + jMax*(vMax - v0))/Abs(jMax);

    // Solution 3/4
    profile.t[0] = h2 - a0/jMax;
    profile.t[1] = 0;
    profile.t[2] = h2;
    profile.t[3] = (af_p3 - a0_p3)/(3*jMax_jMax*vMax) + (a0*v0 - af*vf + (af_af*h1 + a0_a0*h2)/2)/(jMax*vMax) - (v0/vMax + 1.0)*h2 - (vf/vMax + 1.0)*h1 + pd/vMax;
    profile.t[4] = h1;
    profile.t[5] = 0;
    profile.t[6] = profile.t[4] + af/jMax;

    if (profile.check<JerkSigns::UDDU, Limits::VEL>(pf, vf, af, jMax, vMax, aMax, aMin)) {
        add_profile(profile, jMax);
    }
}

void Step1::time_acc0_acc1(Profile& profile, double vMax, double aMax, double aMin, double jMax) {
    const double h1 = Sqrt(3)*Sqrt((aMax - aMin)*(3*Power(af,4)*aMax - 3*Power(a0,4)*aMin + 8*Power(a0,3)*aMax*aMin - 8*Power(af,3)*aMax*aMin - 24*a0*aMax*aMin*jMax*v0 - 6*Power(a0,2)*aMin*(Power(aMax,2) - 2*jMax*v0) + 24*af*aMax*aMin*jMax*vf + 6*Power(af,2)*aMax*(Power(aMin,2) - 2*jMax*vf) + 3*(Power(aMax,3)*Power(aMin,2) - 4*aMin*Power(jMax,2)*Power(v0,2) - Power(aMax,2)*(Power(aMin,3) - 4*aMin*jMax*v0) + 4*aMax*jMax*(2*aMin*jMax*(p0 - pf) - Power(aMin,2)*vf + jMax*Power(vf,2)))))*Abs(jMax)/(3*(aMax - aMin)*jMax);

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

            if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1>(pf, vf, af, jMax, vMax, aMax, aMin)) {
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

            if (profile.check<JerkSigns::UDDU, Limits::ACC0_ACC1>(pf, vf, af, jMax, vMax, aMax, aMin)) {
                add_profile(profile, jMax);
            }
        }
    }
}

void Step1::time_acc1(Profile& profile, double vMax, double aMax, double aMin, double jMax) {
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

        // Single Newton step (regarding pd)
        {
            double h1 = jMax*t*t + v0;
            double orig = -pd + (3*(a0_p4 - af_p4) + 8*af_p3*aMin + 8*a0_p3*(-aMin + 3*jMax*t) + 24*a0*jMax*(-aMin + 2*jMax*t)*(-aMin*t + h1) + 6*a0_a0*(aMin*aMin + 2*jMax*(-4*aMin*t + 5*h1 - 4*v0)) - 24*af*aMin*jMax*vf - 6*af_af*(aMin*aMin - 2*jMax*vf) + 12*jMax*(-2*aMin*jMax*t*(h1 + v0) + aMin*aMin*(h1 + vf) + jMax*(h1*h1 - vf_vf)))/(-24*aMin*jMax_jMax);
            double deriv = -((a0 - aMin + jMax*t)*(a0_a0 - aMin*jMax*t + a0*(-aMin + 4*jMax*t) + 2*jMax*h1))/(aMin*jMax);

            t -= orig / deriv;
        }

        // Corresponds to inverse ACC0
        // if (t < DBL_EPSILON && (-af_af + 2*jMax*(vf - v0) + aMin*aMin) < DBL_EPSILON) {
        //     continue;
        // }

        profile.t[0] = t;
        profile.t[1] = 0;
        profile.t[2] = (a0 - aMin)/jMax + t;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = -((a0_a0 + af_af)/2 - aMin*aMin + jMax*t*(2*a0 + jMax*t) - jMax*(vf - v0))/(aMin*jMax);
        profile.t[6] = (af - aMin)/jMax;

        if (profile.check<JerkSigns::UDDU, Limits::ACC1>(pf, vf, af, jMax, vMax, aMax, aMin)) {
            add_profile(profile, jMax);
        }
    }
}

void Step1::time_acc0(Profile& profile, double vMax, double aMax, double aMin, double jMax) {
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

        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (a0_a0 - af_af + 2*jMax*(-2*aMax*t + jMax*t*t + vf - v0))/(2*aMax*jMax);
        profile.t[2] = t;
        profile.t[3] = 0;
        profile.t[4] = 0;
        profile.t[5] = 0;
        profile.t[6] = (af - aMax)/jMax + t;

        if (profile.check<JerkSigns::UDDU, Limits::ACC0>(pf, vf, af, jMax, vMax, aMax, aMin)) {
            add_profile(profile, jMax);
        }
    }
}

void Step1::time_none(Profile& profile, double vMax, double aMax, double aMin, double jMax) {
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

            if (profile.check<JerkSigns::UDDU, Limits::NONE>(pf, vf, af, jMax, vMax, aMax, aMin)) {
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

                if (profile.check<JerkSigns::UDDU, Limits::NONE>(pf, vf, af, jMax, vMax, aMax, aMin)) {
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

                if (profile.check<JerkSigns::UDDU, Limits::NONE>(pf, vf, af, jMax, vMax, aMax, aMin)) {
                    add_profile(profile, jMax);
                }
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
        polynom[3] = (4*(a0_p3 - af_p3 + 3*jMax*(af*vf - a0*v0 - jMax*pd)))/(3*jMax*jMax_jMax);
        polynom[4] = -Power(a0_a0 - af_af + 2*jMax*(vf - v0),2)/(4*jMax_jMax*jMax_jMax);

        auto roots = Roots::solveQuartMonic(polynom);
        for (double t: roots) {
            if (t < 0.0) {
                continue;
            }

            profile.t[0] = (a0_a0 - af_af + 2*jMax*(jMax*t*t - 2*a0*t - v0 + vf))/(4*jMax_jMax*t);
            profile.t[1] = 0;
            profile.t[2] = t;
            profile.t[3] = 0;
            profile.t[4] = 0;
            profile.t[5] = 0;
            profile.t[6] = (af - a0)/jMax + t - profile.t[0];

            if (profile.check<JerkSigns::UDDU, Limits::NONE>(pf, vf, af, jMax, vMax, aMax, aMin)) {
                add_profile(profile, jMax);
            }
        }
    }
}

bool Step1::calculate_block(Block& block) const {
    // if (valid_profile_counter > 0)
    // {
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
    profile.af = af;
    profile.vf = vf;
    profile.pf = pf;
    valid_profile_counter = 0;

    if (std::abs(pf - p0) < DBL_EPSILON && std::abs(v0) < DBL_EPSILON && std::abs(vf) < DBL_EPSILON && std::abs(a0) < DBL_EPSILON && std::abs(af) < DBL_EPSILON) {
        time_none(profile, vMax, aMax, aMin, jMax);
    
    } else {
        time_acc0_acc1_vel(profile, vMax, aMax, aMin, jMax);
        time_acc0_acc1_vel(profile, vMin, aMin, aMax, -jMax);
        time_acc1_vel(profile, vMax, aMax, aMin, jMax);
        time_acc1_vel(profile, vMin, aMin, aMax, -jMax);
        time_acc0_vel(profile, vMax, aMax, aMin, jMax);
        time_acc0_vel(profile, vMin, aMin, aMax, -jMax);
        time_vel(profile, vMax, aMax, aMin, jMax);
        time_vel(profile, vMin, aMin, aMax, -jMax);
        time_none(profile, vMax, aMax, aMin, jMax);
        time_acc0(profile, vMax, aMax, aMin, jMax);
        time_acc1(profile, vMax, aMax, aMin, jMax);
        time_acc0_acc1(profile, vMax, aMax, aMin, jMax);
        time_none(profile, vMin, aMin, aMax, -jMax);
        time_acc0(profile, vMin, aMin, aMax, -jMax);
        time_acc1(profile, vMin, aMin, aMax, -jMax);
        time_acc0_acc1(profile, vMin, aMin, aMax, -jMax);
    }

    return calculate_block(block);
}

} // namespace ruckig
