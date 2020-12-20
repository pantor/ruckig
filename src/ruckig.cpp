#include <complex>
#include <iomanip>

#include <ruckig/ruckig.hpp>
#include <ruckig/roots.hpp>
#include <ruckig/wolfram.hpp>


namespace ruckig {

void Profile::set(double p0, double v0, double a0, const std::array<double, 7>& j) {
    this->j = j;
    t_sum[0] = t[0];
    a[0] = a0;
    v[0] = v0;
    p[0] = p0;

    for (size_t i = 0; i < 6; i += 1) {
        t_sum[i+1] = t_sum[i] + t[i+1];
    }
    for (size_t i = 0; i < 7; i += 1) {
        a[i+1] = a[i] + t[i] * j[i];
        v[i+1] = v[i] + t[i] * a[i] + std::pow(t[i], 2) * j[i] / 2;
        p[i+1] = p[i] + t[i] * v[i] + std::pow(t[i], 2) * a[i] / 2 + std::pow(t[i], 3) * j[i] / 6;
    }
}

bool Profile::check(double pf, double vf, double vMax, double aMax) const {
    // Velocity and acceleration limits can be broken in the beginnging if the initial velocity and acceleration are too high
    // std::cout << std::setprecision(15) << "target: " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << std::endl;
    return std::all_of(t.begin(), t.end(), [](double tm){ return tm >= 0; })
        && std::all_of(v.begin() + 3, v.end(), [vMax](double vm){ return std::abs(vm) < std::abs(vMax) + 1e-9; })
        && std::all_of(a.begin() + 2, a.end(), [aMax](double am){ return std::abs(am) < std::abs(aMax) + 1e-9; })
        && std::abs(p[7] - pf) < 2e-7 && std::abs(v[7] - vf) < 2e-7;
}

std::tuple<double, double, double> Profile::integrate(double t, double p0, double v0, double a0, double j) {
    const double p_new = p0 + t * (v0 + t * (a0 / 2 + t * j / 6));
    const double v_new = v0 + t * (a0 + t * j / 2);
    const double a_new = a0 + t * j;
    return {p_new, v_new, a_new};
}

bool RuckigStep1::time_up_acc0_acc1_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (Power(a0,2) - 2*Power(aMax,2) - 2*jMax*v0 + 2*jMax*vMax)/(2*aMax*jMax);
    profile.t[2] = aMax/jMax;
    profile.t[3] = (3*Power(a0,4) - 8*Power(a0,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) - 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(v0 + vf + 2*vMax) - jMax*(Power(v0,2) + Power(vf,2) - 2*Power(vMax,2))))/(24.*aMax*Power(jMax,2)*vMax);
    profile.t[4] = profile.t[2];
    profile.t[5] = -aMax/jMax + (vMax-vf)/aMax;
    profile.t[6] = profile.t[2];

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(pf, vf, vMax, aMax);
}

bool RuckigStep2::time_up_acc0_acc1_vel(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    // Profile UDDU
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (3*aMax*jMax*(Power(a0,2) + 2*a0*aMax - 6*Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*(-v0 + vf)) - Sqrt(3)*Sqrt(Power(aMax,2)*(-3*Power(a0,4) + 4*Power(a0,3)*aMax + 12*Power(a0,2)*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)) - 24*a0*aMax*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)) + 12*(Power(aMax,4) - 2*Power(aMax,3)*jMax*tf + Power(aMax,2)*Power(jMax,2)*Power(tf,2) - Power(jMax,2)*Power(v0 - vf,2) + 2*aMax*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf)))))*Abs(jMax))/(12.*Power(aMax,2)*Power(jMax,2));
        profile.t[2] = aMax/jMax;
        profile.t[3] = (-6*Power(aMax,3)*jMax + Sqrt(3)*Sqrt(Power(aMax,2)*(-3*Power(a0,4) + 4*Power(a0,3)*aMax + 12*Power(a0,2)*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)) - 24*a0*aMax*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)) + 12*(Power(aMax,4) - 2*Power(aMax,3)*jMax*tf + Power(aMax,2)*Power(jMax,2)*Power(tf,2) - Power(jMax,2)*Power(v0 - vf,2) + 2*aMax*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf)))))*Abs(jMax))/(6.*Power(aMax,2)*Power(jMax,2));
        profile.t[4] = profile.t[2];
        profile.t[5] = -(3*aMax*jMax*(Power(a0,2) - 2*a0*aMax + 6*Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(-v0 + vf)) + Sqrt(3)*Sqrt(Power(aMax,2)*(-3*Power(a0,4) + 4*Power(a0,3)*aMax + 12*Power(a0,2)*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)) - 24*a0*aMax*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)) + 12*(Power(aMax,4) - 2*Power(aMax,3)*jMax*tf + Power(aMax,2)*Power(jMax,2)*Power(tf,2) - Power(jMax,2)*Power(v0 - vf,2) + 2*aMax*Power(jMax,2)*(2*p0 - 2*pf + tf*(v0 + vf)))))*Abs(jMax))/(12.*Power(aMax,2)*Power(jMax,2));
        profile.t[6] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    // Profile UDUD
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = (3*Power(a0,4) - 4*Power(a0,3)*aMax + 24*a0*Power(aMax,3) - 6*Power(a0,2)*(Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*(v0 - vf)) - 12*(2*Power(aMax,4) - 2*Power(aMax,3)*jMax*tf - 2*aMax*Power(jMax,2)*(p0 - pf + tf*v0) - Power(jMax,2)*Power(v0 - vf,2) + Power(aMax,2)*jMax*(-v0 + vf)))/(12.*aMax*jMax*(Power(a0,2) - 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf))));
        profile.t[2] = aMax/jMax;
        profile.t[3] = -(Power(a0,2) - 2*a0*aMax + 4*Power(aMax,2) - 2*aMax*jMax*tf - 2*jMax*v0 + 2*jMax*vf)/(2.*aMax*jMax);
        profile.t[4] = profile.t[2];
        profile.t[5] = (3*Power(a0,4) - 8*Power(a0,3)*aMax - 6*Power(a0,2)*(Power(aMax,2) + 2*jMax*(v0 - vf)) + 24*a0*(Power(aMax,3) + aMax*jMax*(v0 - vf)) - 12*(2*Power(aMax,4) - 2*Power(aMax,3)*jMax*tf - Power(jMax,2)*Power(v0 - vf,2) + Power(aMax,2)*jMax*(-v0 + vf) + 2*aMax*Power(jMax,2)*(p0 - pf + tf*vf)))/(12.*aMax*jMax*(Power(a0,2) - 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf))));
        profile.t[6] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }
    
    return false;
}

bool RuckigStep1::time_up_acc1_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    profile.t[0] = (-2*a0*jMax + Sqrt(2)*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*Abs(jMax))/(2*Power(jMax,2));
    profile.t[1] = 0;
    profile.t[2] = Sqrt(Power(a0,2)/2 + jMax*(-v0 + vMax))/Abs(jMax);
    profile.t[3] = (-2*jMax*(2*Power(a0,3)*aMax - 6*a0*aMax*jMax*v0 + 3*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(vf + vMax) + jMax*(-Power(vf,2) + Power(vMax,2)))) + 3*Sqrt(2)*aMax*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*(Power(a0,2) - 2*jMax*(v0 + vMax))*Abs(jMax))/(12.*aMax*Power(jMax,3)*vMax);
    profile.t[4] = aMax/jMax;
    profile.t[5] = (-Power(aMax,2)/jMax - vf + vMax)/aMax;
    profile.t[6] = profile.t[4];

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(pf, vf, vMax, aMax);
}

bool RuckigStep2::time_up_acc1_vel(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    // Profile UDDU
    {
        std::array<double, 5> polynom;
        polynom[0] = 12*Power(jMax,4);
        polynom[1] = 24*(2*a0 + aMax)*Power(jMax,3);
        polynom[2] = 12*Power(jMax,2)*(5*Power(a0,2) + 4*a0*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(v0 - vf));
        polynom[3] = 24*a0*jMax*(Power(a0,2) + a0*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(v0 - vf));
        polynom[4] = 3*Power(a0,4) + 4*Power(a0,3)*aMax + 6*Power(a0,2)*(Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(v0 - vf)) + 12*jMax*(-2*aMax*jMax*(p0 - pf + tf*v0) + Power(aMax,2)*(v0 - vf) + jMax*Power(v0 - vf,2));

        auto roots = Roots::solveQuart(polynom);
        for (double t: roots) {
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + t;
            profile.t[3] = -(Power(a0,2) + 2*a0*(aMax + 2*jMax*t) + 2*(Power(aMax,2) + aMax*jMax*(2*t - tf) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[4] = aMax/jMax;
            profile.t[5] = (Power(a0,2) - 2*Power(aMax,2) + 4*a0*jMax*t + 2*jMax*(jMax*Power(t,2) + v0 - vf))/(2.*aMax*jMax);
            profile.t[6] = profile.t[4];
            
            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        std::array<double, 5> polynom;
        polynom[0] = 12*Power(jMax,4);
        polynom[1] = 24*(2*a0 - aMax)*Power(jMax,3);
        polynom[2] = 12*Power(jMax,2)*(5*Power(a0,2) - 4*a0*aMax - Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*(v0 - vf));
        polynom[3] = 24*a0*jMax*(Power(a0,2) - a0*aMax - Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*(v0 - vf));
        polynom[4] = 3*Power(a0,4) - 4*Power(a0,3)*aMax + 12*jMax*(2*aMax*jMax*(p0 - pf + tf*v0) + jMax*Power(v0 - vf,2) + Power(aMax,2)*(-v0 + vf)) - 6*Power(a0,2)*(Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(-v0 + vf));

        auto roots = Roots::solveQuart(polynom);
        for (double t: roots) {
            profile.t[0] = t;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + t;
            profile.t[3] = (Power(a0,2) - 2*Power(aMax,2) - 2*a0*(aMax - 2*jMax*t) + 2*aMax*jMax*(-2*t + tf) + 2*jMax*(jMax*Power(t,2) + v0 - vf))/(2.*aMax*jMax);
            profile.t[4] = aMax/jMax;
            profile.t[5] = -(Power(a0,2) + 4*a0*jMax*t + 2*(Power(aMax,2) + jMax*(jMax*Power(t,2) + v0 - vf)))/(2.*aMax*jMax);
            profile.t[6] = profile.t[4];
            
            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool RuckigStep1::time_up_acc0_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (Power(a0,2) - 2*Power(aMax,2) + 2*jMax*(-v0 + vMax))/(2*aMax*jMax);
    profile.t[2] = aMax/jMax;
    profile.t[3] = ((3*Power(a0,4) - 8*Power(a0,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) - 12*jMax*(Power(aMax,2)*(v0 + vMax) + jMax*(-Power(v0,2) + Power(vMax,2)) + 2*aMax*(jMax*(p0 - pf) + SqrtComplex(jMax)*SqrtComplex(-vf + vMax)*(vf + vMax))))/(24.*aMax*Power(jMax,2)*vMax)).real();
    profile.t[4] = Sqrt((-vf + vMax)/jMax);
    profile.t[5] = 0;
    profile.t[6] = profile.t[4];

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(pf, vf, vMax, aMax);
}

bool RuckigStep2::time_up_acc0_vel(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    const double h1 = Power(a0,2) - 2*a0*aMax + Power(aMax,2) - 2*aMax*jMax*tf + 2*jMax*(-v0 + vf);
    const double h2 = 3*Power(a0,4) - 8*Power(a0,3)*aMax + 24*a0*aMax*jMax*(v0 - vf) + 6*Power(a0,2)*(Power(aMax,2) + 2*jMax*(-v0 + vf)) - 12*jMax*(Power(aMax,2)*(v0 - vf) - jMax*Power(v0 - vf,2) + 2*aMax*jMax*(p0 - pf + tf*vf));
    const double h3 = 144*(Power(h1,2) + h2)*Power(jMax,4);
        
    // Profile UDDU
    {
        const double h4 = 1728*(2*Power(h1,3) + 9*Power(aMax,2)*h2 - 6*h1*h2)*Power(jMax,6);
        const auto h5 = PowerComplex(h4 + SqrtComplex(-4*Power(h3,3) + Power(h4,2)),1./3);
        const auto h6 = SqrtComplex((2*Power(2,1./3)*h3 + h5*(Power(2,2./3)*h5 + 24*(3*Power(aMax,2) - 2*h1)*Power(jMax,2)))/(h5*Power(jMax,4)))/(6.*Sqrt(2));
        const auto h7 = (-2*(Power(aMax,3) - aMax*h1))/(h6*Power(jMax,3));
        const auto h8 = -(2*Power(2,1./3)*h3 + h5*(Power(2,2./3)*h5 + 48*(-3*Power(aMax,2) + 2*h1)*Power(jMax,2)))/(72.*h5*Power(jMax,4));

        // Solution 3
        {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = ((Power(a0,2)*aMax - 2*a0*aMax*(aMax - h6*jMax) - jMax*(PowerComplex(h6,2)*SqrtComplex(h7 + h8)*Power(jMax,2) + Power(aMax,2)*(h6 + 2*tf) + aMax*(PowerComplex(h6,2)*jMax - h6*jMax*(SqrtComplex(h7 + h8) + 2*tf) + 2*(v0 - vf))))/(2.*aMax*h6*Power(jMax,2))).real();
            profile.t[2] = aMax/jMax;
            profile.t[3] = ((-(Power(a0,2)*aMax) + 2*a0*Power(aMax,2) + PowerComplex(h6,2)*SqrtComplex(h7 + h8)*Power(jMax,3) - Power(aMax,2)*jMax*(h6 - 2*tf) + aMax*jMax*(-(PowerComplex(h6,2)*jMax) + h6*SqrtComplex(h7 + h8)*jMax + 2*v0 - 2*vf))/(2.*aMax*h6*Power(jMax,2))).real();
            profile.t[4] = -(aMax - h6*jMax + SqrtComplex(h7 + h8)*jMax).real()/(2.*jMax);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            // std::cout << profile.t[0] << std::endl;
            // std::cout << profile.t[1] << std::endl;
            // std::cout << profile.t[2] << std::endl;
            // std::cout << profile.t[3] << std::endl;
            // std::cout << profile.t[4] << std::endl;
            // std::cout << profile.t[5] << std::endl;
            // std::cout << profile.t[6] << std::endl;
            
            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }
    }

    // Profile UDUD
    {
        const double h4 = 1728*(-2*Power(h1,3) + 9*Power(aMax,2)*h2 + 6*h1*h2)*Power(jMax,6);
        const auto h5 = PowerComplex(h4 + SqrtComplex(-4*Power(h3,3) + Power(h4,2)),1./3);
        const auto h6 = SqrtComplex((2*Power(2,0.3333333333333333)*h3 + h5*(Power(2,0.6666666666666666)*h5 + 24*(3*Power(aMax,2) + 2*h1)*Power(jMax,2)))/(h5*Power(jMax,4)))/(6.*Sqrt(2));
        const auto h7 = (2*aMax*(Power(aMax,2) + h1))/(h6*Power(jMax,3));
        const auto h8 = (-2*Power(2,0.3333333333333333)*h3 + h5*(-(Power(2,0.6666666666666666)*h5) + 48*(3*Power(aMax,2) + 2*h1)*Power(jMax,2)))/(72.*h5*Power(jMax,4));

        // Solution 2
        {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = ((Power(a0,2)*aMax + 2*Power(aMax,3) + PowerComplex(h6,2)*SqrtComplex(-h7 + h8)*Power(jMax,3) - 2*a0*aMax*(aMax - h6*jMax) - Power(aMax,2)*jMax*(5.*h6 + 2*tf) + aMax*jMax*(PowerComplex(h6,2)*jMax - h6*SqrtComplex(-h7 + h8)*jMax + 2.*h6*jMax*tf - 2*v0 + 2*vf))/(2.*aMax*h6*Power(jMax,2))).real();
            profile.t[2] = aMax/jMax;
            profile.t[3] = (-(Power(a0,2)*aMax - 2*a0*Power(aMax,2) + 2*Power(aMax,3) + PowerComplex(h6,2)*SqrtComplex(-h7 + h8)*Power(jMax,3) + Power(aMax,2)*jMax*(h6 - 2*tf) + aMax*jMax*(-(PowerComplex(h6,2)*jMax) + h6*SqrtComplex(-h7 + h8)*jMax - 2*v0 + 2*vf))/(2.*aMax*h6*Power(jMax,2))).real();
            profile.t[4] = (aMax - h6*jMax + SqrtComplex(-h7 + h8)*jMax).real()/(2.*jMax);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            // std::cout << profile.t[0] << std::endl;
            // std::cout << profile.t[1] << std::endl;
            // std::cout << profile.t[2] << std::endl;
            // std::cout << profile.t[3] << std::endl;
            // std::cout << profile.t[4] << std::endl;
            // std::cout << profile.t[5] << std::endl;
            // std::cout << profile.t[6] << std::endl;
            // std::cout << "---" << std::endl;
            
            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool RuckigStep1::time_up_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    profile.t[0] = ((-2*a0*jMax + Sqrt(2)*Sqrt(Power(a0,2) + 2*jMax*(-v0 + vMax))*Abs(jMax))/(2*Power(jMax,2)));
    profile.t[1] = 0;
    profile.t[2] = Sqrt(Power(a0,2)/2 + jMax*(-v0 + vMax))/Abs(jMax);
    profile.t[3] = ((-4*jMax*(Power(a0,3) + 3*Power(jMax,2)*(p0 - pf) - 3*a0*jMax*v0 + 3*jMax*SqrtComplex(jMax)*SqrtComplex(-vf + vMax)*(vf + vMax)) + 3*Sqrt(2)*SqrtComplex(Power(a0,2) + 2*jMax*(-v0 + vMax))*(Power(a0,2) - 2*jMax*(v0 + vMax))*Abs(jMax))/(12.*Power(jMax,3)*vMax)).real();
    profile.t[4] = Sqrt((-vf + vMax)/jMax);
    profile.t[5] = 0;
    profile.t[6] = profile.t[4];

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(pf, vf, vMax, aMax);
}

bool RuckigStep2::time_up_vel(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
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

        // Solve 4th order derivative analytically
        auto extremas = Roots::solveQuart(5 * polynom[0], 4 * polynom[1], 3 * polynom[2], 2 * polynom[3], polynom[4]);
        std::set<std::tuple<double, double>> tz_intervals;

        double tz_min {0.0};
        double tz_max = std::min<double>(tf, (tf - a0/jMax) / 2);
        double tz_current {tz_min};

        for (double tz: extremas) {
            if (tz <= 0.0 || tz >= tz_max) {
                continue;
            }

            // Check that polynom(lower) and polynom(upper) have different signs (should only happen at first and last boundary)
            // std::cout << "tz_current: " << tz_current << " " << eval(polynom, tz_current) << std::endl;
            // std::cout << "tz: " << tz << " " << eval(polynom, tz) << std::endl;
            double val_current = eval(polynom, tz_current);
            double val_new = eval(polynom, tz);
            if (std::abs(val_new) < 1e-15) {
                tz_intervals.insert({tz - 1e-12, tz + 1e-12});
                tz += 1e-14;
            } else if (val_current * val_new < 0) {
                tz_intervals.insert({tz_current, tz});
            }
            tz_current = tz;
        }
        if (eval(polynom, tz_current) * eval(polynom, tz_max) < 0) {
            tz_intervals.insert({tz_current, tz_max});
        }

        for (auto interval: tz_intervals) {
            // Use safe Newton method
            double lower = std::get<0>(interval);
            double upper = std::get<1>(interval);
            double tz = Roots::shrinkInterval(polynom.data(), 6, lower, upper, 1e-16);

            double vPlat = Power(a0,2)/(2.*jMax) + 2*a0*tz + jMax*Power(tz,2) + v0;
            double h1 = 2*Sqrt((vPlat - vf)/jMax);

            profile.t[0] = tz;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + tz;
            profile.t[3] = -(4*Power(a0,3) + 12*a0*jMax*(v0 + jMax*tz*(3*tz + h1)) + 3*Power(a0,2)*jMax*(8*tz + h1) + 6*Power(jMax,2)*(2*p0 - 2*pf + 2*jMax*Power(tz,3) + 4*tz*v0 + jMax*Power(tz,2)*h1 + v0*h1 + h1*vf))/(6.*jMax*(Power(a0,2) + 4*a0*jMax*tz + 2*jMax*(jMax*Power(tz,2) + v0)));;
            profile.t[4] = Sqrt((-vf + vPlat)/jMax);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            // std::cout << profile.t[0] << std::endl;
            // std::cout << profile.t[1] << std::endl;
            // std::cout << profile.t[2] << std::endl;
            // std::cout << profile.t[3] << std::endl;
            // std::cout << profile.t[4] << std::endl;
            // std::cout << profile.t[5] << std::endl;
            // std::cout << profile.t[6] << std::endl;
            // std::cout << "---" << std::endl;

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
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
            if (eval(deriv, dd_tz_current) * eval(deriv, tz) < 0) {
                dd_tz_intervals.insert({dd_tz_current, tz});
            }
            dd_tz_current = tz;
        }
        if (eval(deriv, dd_tz_current) * eval(deriv, tz_max) < 0) {
            dd_tz_intervals.insert({dd_tz_current, tz_max});
        }

        std::set<std::tuple<double, double>> tz_intervals;
        double tz_current {tz_min};

        for (auto interval: dd_tz_intervals) {
            double lower = std::get<0>(interval);
            double upper = std::get<1>(interval);
            double tz = Roots::shrinkInterval(deriv.data(), 6, lower, upper, 1e-14);

            if (tz <= 0.0 || tz >= tz_max) {
                continue;
            }
            // Check that polynom(lower) and polynom(upper) have different signs (should only happen at first and last boundary)
            if (eval(polynom, tz_current) * eval(polynom, tz) < 0) {
                tz_intervals.insert({tz_current, tz});
            }
            tz_current = tz;
        }
        if (eval(polynom, tz_current) * eval(polynom, tz_max) < 0) {
            tz_intervals.insert({tz_current, tz_max});
        }

        for (auto interval: tz_intervals) {
            // Use safe Newton method
            double lower = std::get<0>(interval);
            double upper = std::get<1>(interval);
            double tz = Roots::shrinkInterval(polynom.data(), 7, lower, upper, 1e-14);

            double vPlat = Power(a0,2)/(2.*jMax) + 2*a0*tz + jMax*Power(tz,2) + v0;
            // std::cout << "BE CAREFUL " << vPlat << std::endl;

            profile.t[0] = tz;
            profile.t[1] = 0;
            profile.t[2] = a0/jMax + tz;
            profile.t[3] = tf - (2*tz + a0/jMax + 2 * Sqrt((vf - vPlat)/jMax));
            profile.t[4] = Sqrt((vf - vPlat)/jMax); 
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            // std::cout << profile.t[4] << " " << Sqrt((vf - vPlat)/jMax) << std::endl;

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }
    }

    return false;
}

bool RuckigStep1::time_up_acc0_acc1(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    const double h1 = Abs(aMax)*Abs(jMax)*Sqrt(6*(3*Power(a0,4) - 8*Power(a0,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) + 6*(Power(aMax,4) + 4*aMax*Power(jMax,2)*(-p0 + pf) - 2*Power(aMax,2)*jMax*(v0 + vf) + 2*Power(jMax,2)*(Power(v0,2) + Power(vf,2)))));

    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = (6*Power(a0,2)*aMax*jMax - 18*Power(aMax,3)*jMax - 12*aMax*Power(jMax,2)*v0 + h1)/(12.*Power(aMax,2)*Power(jMax,2));
    profile.t[2] = aMax/jMax;
    profile.t[3] = 0;
    profile.t[4] = profile.t[2];
    profile.t[5] = (-18*Power(aMax,3)*jMax - 12*aMax*Power(jMax,2)*vf + h1)/(12.*Power(aMax,2)*Power(jMax,2));
    profile.t[6] = profile.t[2];

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(pf, vf, vMax, aMax);
}

bool RuckigStep2::time_up_acc0_acc1(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    if (std::abs(a0) < 1e-16) {
        profile.t[0] = (Power(aMax,2)*Power(tf,2) - Power(v0 - vf,2) + 2*aMax*(2*p0 - 2*pf + tf*(v0 + vf)))/(2.*Power(aMax,2)*tf);
        profile.t[1] = -(Power(aMax,2)*Power(tf,2) - 2*Power(v0 - vf,2) + aMax*(8*p0 - 8*pf + 5*tf*v0 + 3*tf*vf))/(2.*Power(aMax,2)*tf);
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = -(Power(aMax,2)*Power(tf,2) - 2*Power(v0 - vf,2) + aMax*(8*p0 - 8*pf + 3*tf*v0 + 5*tf*vf))/(2.*Power(aMax,2)*tf);
        profile.t[6] = profile.t[0];
        jMax = aMax/profile.t[0];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        return profile.check(pf, vf, vMax, aMax);
    }
    
    double h1 = -12*(2*Power(aMax,3)*tf + Power(a0,2)*(aMax*tf - v0 + vf) - 2*a0*aMax*(aMax*tf - v0 + vf));
    double h2 = -12*(Power(aMax,2)*Power(tf,2) - Power(v0 - vf,2) + 2*aMax*(2*p0 - 2*pf + tf*(v0 + vf)));
    double h3 = 3*Power(a0,3) - 4*Power(a0,2)*aMax - 12*a0*Power(aMax,2) + 24*Power(aMax,3);
    double h4 = h1 - Sqrt(Power(h1,2) - 4*a0*h2*h3);
    double h5 = h1 + Sqrt(Power(h1,2) - 4*a0*h2*h3);
    double h6 = 2*aMax*h4*(2*p0 - 2*pf + tf*(v0 + vf));
    double h7 = 2*a0*aMax*h2*(a0*tf + 2*v0 - 2*vf);

    profile.t[0] = (6*(a0 - aMax)*(-h6 + h7 + 4*Power(aMax,3)*h2*tf - Power(aMax,2)*tf*(4*a0*h2 + h4*tf) + (-2*Power(a0,2)*h2 + h4*(v0 - vf))*(v0 - vf)))/(a0*h2*h3);
    profile.t[1] = (24*Power(aMax,2)*(-h6 + 4*Power(aMax,3)*h2*tf - Power(aMax,2)*h4*Power(tf,2) + h4*Power(v0 - vf,2)) - 6*a0*aMax*(-h6 + 16*Power(aMax,3)*h2*tf - Power(aMax,2)*(h4*Power(tf,2) + 12*h2*(v0 - vf)) + h4*Power(v0 - vf,2)) - 3*Power(a0,4)*h2*(aMax*tf - v0 + vf) - 4*Power(a0,3)*aMax*h2*(aMax*tf - v0 + vf) + 3*Power(a0,2)*(h6 + 16*Power(aMax,3)*h2*tf - h4*Power(v0 - vf,2) + Power(aMax,2)*(h4*Power(tf,2) + 20*h2*(-v0 + vf))))/(2.*a0*aMax*h2*h3);
    profile.t[2] = (-6*aMax*(-h6 + h7 + 4*Power(aMax,3)*h2*tf - Power(aMax,2)*tf*(4*a0*h2 + h4*tf) + (-2*Power(a0,2)*h2 + h4*(v0 - vf))*(v0 - vf)))/(a0*h2*h3);
    profile.t[3] = 0;
    profile.t[4] = profile.t[2];
    profile.t[5] = (24*Power(aMax,2)*(-h6 + 4*Power(aMax,3)*h2*tf - Power(aMax,2)*h4*Power(tf,2) + h4*Power(v0 - vf,2)) - 6*a0*aMax*(-h6 + 16*Power(aMax,3)*h2*tf - Power(aMax,2)*(h4*Power(tf,2) + 20*h2*(v0 - vf)) + h4*Power(v0 - vf,2)) + 3*Power(a0,2)*(-h6 + 24*Power(aMax,3)*h2*tf - Power(aMax,2)*(h4*Power(tf,2) + 28*h2*(v0 - vf)) + h4*Power(v0 - vf,2)) + 3*Power(a0,4)*h2*(3*aMax*tf - v0 + vf) - 4*Power(a0,3)*aMax*h2*(7*aMax*tf - 5*v0 + 5*vf))/(2.*a0*aMax*h2*h3);
    profile.t[6] = profile.t[2];

    jMax = h4/(2.*h2);

    // std::cout << profile.t[0] << std::endl;
    // std::cout << profile.t[1] << std::endl;
    // std::cout << profile.t[2] << std::endl;
    // std::cout << profile.t[4] << std::endl;
    // std::cout << profile.t[5] << std::endl;
    // std::cout << jMax << std::endl;

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    return profile.check(pf, vf, vMax, aMax);
}

bool RuckigStep1::time_up_acc1(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    // if (a0 < 2e-15 && v0 < 2e-15) {

    // }
    
    const double h1 = 5*Power(a0,2) + 6*a0*aMax + Power(aMax,2) + 2*jMax*v0;
    const double h2 = 2*a0 + aMax;
    const double h3 = 3*Power(a0,4) + 8*Power(a0,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) + 2*jMax*v0) + 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(v0 + vf) + jMax*(Power(v0,2) - Power(vf,2)));
    const double h4 = (a0 + aMax)*(Power(a0,2) + a0*aMax + 2*jMax*v0);
    const double h5 = 4*Power(a0,4) + 8*Power(a0,3)*aMax + Power(aMax,4) + 24*aMax*Power(jMax,2)*(p0 - pf) - 24*a0*aMax*jMax*v0 + 4*Power(a0,2)*(Power(aMax,2) - 4*jMax*v0) + Power(aMax,2)*jMax*(-8*v0 + 12*vf) + 4*Power(jMax,2)*(4*Power(v0,2) - 3*Power(vf,2));
    const double h6 = (2*Power(h1,3) - 6*h1*(h3 + 6*h2*h4) + 9*(Power(h2,2)*h3 + 12*Power(h4,2)));
    const auto h7 = 12*Power(jMax,2)*PowerComplex(h6 + SqrtComplex(Power(h6,2) - 4*Power(h5,3)),1./3);
    const auto h8 = SqrtComplex((4*Power(2,1./3)*h5)/h7 + (Power(2,2./3)*h7 + 24*(-2*h1 + 3*Power(h2,2))*Power(jMax,2))/(72.*Power(jMax,4)));
    const auto h9 = SqrtComplex((-576*Power(2,1./3)*h5)/h7 - (2*Power(2,2./3)*h7)/Power(jMax,4) - (96.*(h1*(3*h2 + 2.*h8*jMax) - 3.*(Power(h2,3) + 2*h4 + Power(h2,2)*h8*jMax)))/(h8*Power(jMax,3)));
    const auto h10 = SqrtComplex((-576*Power(2,1./3)*h5)/h7 - (2*Power(2,2./3)*h7)/Power(jMax,4) + (96.*(3*h1*h2 - 3*Power(h2,3) - 6*h4 - 2.*h1*h8*jMax + 3*Power(h2,2)*h8*jMax))/(h8*Power(jMax,3)));

    // Solution 2
    {
        profile.t[0] = -h2/(2*jMax) + (h9 - 12.*h8).real()/24.;
        profile.t[1] = 0;
        profile.t[2] = -aMax/(2*jMax) + (h9 - 12.*h8).real()/24.;
        profile.t[3] = 0;
        profile.t[4] = aMax/jMax;
        profile.t[5] = -((12*Power(a0,2)*aMax + jMax*(12*Power(aMax,2)*h8 + aMax*(-12.*PowerComplex(h8,2)*jMax + h8*jMax*h9 - 24*v0) + h8*jMax*(h8*jMax*h9 + 24*vf)))/(24.*aMax*h8*Power(jMax,2))).real();
        profile.t[6] = profile.t[4];

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];

        // std::cout << "HERE" << std::endl;

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    // Solution 4
    {
        profile.t[0] = -h2/(2.*jMax) + (h10 + 12.*h8).real()/24.;
        profile.t[1] = 0;
        profile.t[2] = -aMax/(2.*jMax) + (h10 + 12.*h8).real()/24.;
        profile.t[3] = 0;
        profile.t[4] = aMax/jMax;
        profile.t[5] = -(Power(aMax,2) + aMax*(-(Power(a0,2)/(h8*jMax)) + h10*jMax/12. + h8*jMax + (2*v0)/h8) + jMax*(-(h10*h8*jMax)/12. + 2*vf)).real()/(2.*aMax*jMax);
        profile.t[6] = profile.t[4];

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    return false;
}

bool RuckigStep2::time_up_acc1(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    // a3 != 0
    // Case UDDU, Solution 2
    {
        profile.t[0] = 0;
        profile.t[1] = 0;
        profile.t[2] = -(-2*Power(a0,3)*jMax - 6*Power(a0,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf - 6*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*Power(jMax,3)*tf*vf - Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) + 3*Power(a0,2)*aMax + 3*a0*Power(aMax,2) + 3*Power(aMax,2)*jMax*tf - 3*Power(jMax,2)*(2*p0 - 2*pf + tf*(aMax*tf + 2*vf)),2) - 3*(Power(a0,2) + 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)))*(Power(a0,4) + 4*Power(a0,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 12*jMax*(-2*aMax*jMax*(p0 - pf + tf*v0) + Power(aMax,2)*(v0 - vf) + jMax*Power(v0 - vf,2)) - 12*a0*jMax*(-(Power(aMax,2)*tf) + jMax*(2*p0 - 2*pf + aMax*Power(tf,2) + 2*tf*vf))))))/(6.*Power(jMax,2)*(Power(a0,2) + 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))));
        profile.t[3] = -(4*Power(a0,3)*jMax + 12*Power(a0,2)*aMax*jMax + 18*a0*Power(aMax,2)*jMax + 12*Power(aMax,3)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf - 12*a0*aMax*Power(jMax,2)*tf - 18*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*a0*Power(jMax,2)*v0 + 12*aMax*Power(jMax,2)*v0 - 12*a0*Power(jMax,2)*vf - 12*aMax*Power(jMax,2)*vf + 12*Power(jMax,3)*tf*vf + Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) + 3*Power(a0,2)*aMax + 3*a0*Power(aMax,2) + 3*Power(aMax,2)*jMax*tf - 3*Power(jMax,2)*(2*p0 - 2*pf + tf*(aMax*tf + 2*vf)),2) - 3*(Power(a0,2) + 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)))*(Power(a0,4) + 4*Power(a0,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 12*jMax*(-2*aMax*jMax*(p0 - pf + tf*v0) + Power(aMax,2)*(v0 - vf) + jMax*Power(v0 - vf,2)) - 12*a0*jMax*(-(Power(aMax,2)*tf) + jMax*(2*p0 - 2*pf + aMax*Power(tf,2) + 2*tf*vf))))))/(6.*Power(jMax,2)*(Power(a0,2) + 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))));
        profile.t[4] = -(-4*Power(a0,3)*jMax - 12*Power(a0,2)*aMax*jMax - 18*a0*Power(aMax,2)*jMax - 12*Power(aMax,3)*jMax - 12*Power(jMax,3)*p0 + 12*Power(jMax,3)*pf + 12*a0*aMax*Power(jMax,2)*tf + 18*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) - 12*a0*Power(jMax,2)*v0 - 12*aMax*Power(jMax,2)*v0 + 12*a0*Power(jMax,2)*vf + 12*aMax*Power(jMax,2)*vf - 12*Power(jMax,3)*tf*vf + Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) + 3*Power(a0,2)*aMax + 3*a0*Power(aMax,2) + 3*Power(aMax,2)*jMax*tf - 3*Power(jMax,2)*(2*p0 - 2*pf + tf*(aMax*tf + 2*vf)),2) - 3*(Power(a0,2) + 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)))*(Power(a0,4) + 4*Power(a0,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 12*jMax*(-2*aMax*jMax*(p0 - pf + tf*v0) + Power(aMax,2)*(v0 - vf) + jMax*Power(v0 - vf,2)) - 12*a0*jMax*(-(Power(aMax,2)*tf) + jMax*(2*p0 - 2*pf + aMax*Power(tf,2) + 2*tf*vf))))))/(6.*Power(jMax,2)*(Power(a0,2) + 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))));
        profile.t[5] = -(2*Power(a0,3)*jMax + 12*Power(a0,2)*aMax*jMax + 18*a0*Power(aMax,2)*jMax + 12*Power(aMax,3)*jMax - 12*Power(jMax,3)*p0 + 12*Power(jMax,3)*pf - 6*Power(a0,2)*Power(jMax,2)*tf - 12*a0*aMax*Power(jMax,2)*tf - 18*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) + 12*aMax*Power(jMax,2)*v0 - 12*Power(jMax,3)*tf*v0 - 12*aMax*Power(jMax,2)*vf - Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) + 3*Power(a0,2)*aMax + 3*a0*Power(aMax,2) + 3*Power(aMax,2)*jMax*tf - 3*Power(jMax,2)*(2*p0 - 2*pf + tf*(aMax*tf + 2*vf)),2) - 3*(Power(a0,2) + 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf)))*(Power(a0,4) + 4*Power(a0,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 12*jMax*(-2*aMax*jMax*(p0 - pf + tf*v0) + Power(aMax,2)*(v0 - vf) + jMax*Power(v0 - vf,2)) - 12*a0*jMax*(-(Power(aMax,2)*tf) + jMax*(2*p0 - 2*pf + aMax*Power(tf,2) + 2*tf*vf))))))/(6.*Power(jMax,2)*(Power(a0,2) + 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(v0 - vf))));
        profile.t[6] = (aMax/jMax);

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    // Case UDUD, Solution 1
    {
        profile.t[0] = 0;
        profile.t[1] = 0;
        profile.t[2] = -(-2*Power(a0,3)*jMax + 6*Power(a0,2)*aMax*jMax - 6*a0*Power(aMax,2)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf + 6*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) + 12*Power(jMax,3)*tf*vf + Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) - 3*Power(a0,2)*aMax + 3*a0*Power(aMax,2) - 3*jMax*(Power(aMax,2)*tf + jMax*(2*p0 - 2*pf - aMax*Power(tf,2) + 2*tf*vf)),2) - 3*(Power(a0,2) - 2*a0*aMax + 2*jMax*(aMax*tf + v0 - vf))*(Power(a0,4) - 4*Power(a0,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 12*jMax*(2*aMax*jMax*(p0 - pf + tf*v0) + jMax*Power(v0 - vf,2) + Power(aMax,2)*(-v0 + vf)) - 12*a0*jMax*(Power(aMax,2)*tf + jMax*(2*p0 - 2*pf - aMax*Power(tf,2) + 2*tf*vf))))))/(6.*Power(jMax,2)*(Power(a0,2) - 2*a0*aMax + 2*jMax*(aMax*tf + v0 - vf)));
        profile.t[3] = (Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) - 3*Power(a0,2)*aMax + 3*a0*Power(aMax,2) - 3*jMax*(Power(aMax,2)*tf + jMax*(2*p0 - 2*pf - aMax*Power(tf,2) + 2*tf*vf)),2) - 3*(Power(a0,2) - 2*a0*aMax + 2*jMax*(aMax*tf + v0 - vf))*(Power(a0,4) - 4*Power(a0,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 12*jMax*(2*aMax*jMax*(p0 - pf + tf*v0) + jMax*Power(v0 - vf,2) + Power(aMax,2)*(-v0 + vf)) - 12*a0*jMax*(Power(aMax,2)*tf + jMax*(2*p0 - 2*pf - aMax*Power(tf,2) + 2*tf*vf))))))/(3.*Power(jMax,2)*(Power(a0,2) - 2*a0*aMax + 2*jMax*(aMax*tf + v0 - vf)));
        profile.t[4] = -(4*Power(a0,3)*jMax - 12*Power(a0,2)*aMax*jMax + 6*a0*Power(aMax,2)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf + 12*a0*aMax*Power(jMax,2)*tf - 6*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) + 12*a0*Power(jMax,2)*v0 - 12*aMax*Power(jMax,2)*v0 - 12*a0*Power(jMax,2)*vf + 12*aMax*Power(jMax,2)*vf + 12*Power(jMax,3)*tf*vf + Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) - 3*Power(a0,2)*aMax + 3*a0*Power(aMax,2) - 3*jMax*(Power(aMax,2)*tf + jMax*(2*p0 - 2*pf - aMax*Power(tf,2) + 2*tf*vf)),2) - 3*(Power(a0,2) - 2*a0*aMax + 2*jMax*(aMax*tf + v0 - vf))*(Power(a0,4) - 4*Power(a0,3)*aMax + 6*Power(a0,2)*Power(aMax,2) + 12*jMax*(2*aMax*jMax*(p0 - pf + tf*v0) + jMax*Power(v0 - vf,2) + Power(aMax,2)*(-v0 + vf)) - 12*a0*jMax*(Power(aMax,2)*tf + jMax*(2*p0 - 2*pf - aMax*Power(tf,2) + 2*tf*vf))))))/(6.*Power(jMax,2)*(Power(a0,2) - 2*a0*aMax + 2*jMax*(aMax*tf + v0 - vf)));
        profile.t[5] = (Power(a0,3) + Power(a0,2)*(-6*aMax + 3*jMax*tf) + 6*a0*(Power(aMax,2) + jMax*(v0 - vf)) - 6*jMax*(aMax*(aMax*tf + 2*v0 - 2*vf) - jMax*(2*p0 - 2*pf + tf*(v0 + vf))))/(3.*jMax*(Power(a0,2) - 2*a0*aMax + 2*jMax*(aMax*tf + v0 - vf)));
        profile.t[6] = (aMax/jMax);

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, jMax, 0, -jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }
    return false;
}

bool RuckigStep1::time_up_acc0(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    if (std::abs(vf) < 1e-16) {
        const auto h1 = Sqrt(3)*SqrtComplex(-3*Power(a0,4) + 8*Power(a0,3)*aMax - 24*a0*aMax*jMax*v0 - 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) + 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*v0 - jMax*Power(v0,2)));
        const auto h2 = SqrtComplex(9*Power(aMax,2) - Complex(0,6)*h1);

        // Solution 2
        {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (3*Power(a0,2) - 3*Power(aMax,2) - 6*jMax*v0 - Complex(0,1)*h1 - aMax*h2).real()/(6.*aMax*jMax);
            profile.t[2] = aMax/jMax;
            profile.t[3] = 0;
            profile.t[4] = (-3*aMax + h2.real())/(6.*jMax);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }

        // Solution 1
        {
            profile.t[0] = (-a0 + aMax)/jMax;
            profile.t[1] = (3*Power(a0,2) - 3*Power(aMax,2) - 6*jMax*v0 - Complex(0,1)*h1 + aMax*h2).real()/(6.*aMax*jMax);
            profile.t[2] = aMax/jMax;
            profile.t[3] = 0;
            profile.t[4] = -(3*aMax + h2.real())/(6.*jMax);
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }

        return false;
    }

    const double h1 = Power(aMax,2) + 2*jMax*vf;
    const double h2 = 3*Power(a0,4) - 8*Power(a0,3)*aMax + 24*a0*aMax*jMax*v0 + 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) - 12*jMax*(2*aMax*jMax*(p0 - pf) + Power(aMax,2)*(v0 + vf) + jMax*(-Power(v0,2) + Power(vf,2)));
    const double h3 = Power(jMax,4)*(-3*Power(a0,4) + 8*Power(a0,3)*aMax + Power(aMax,4) + 24*aMax*Power(jMax,2)*(p0 - pf) - 24*a0*aMax*jMax*v0 - 6*Power(a0,2)*(Power(aMax,2) - 2*jMax*v0) + 4*Power(aMax,2)*jMax*(3*v0 - 2*vf) + 4*Power(jMax,2)*(-3*Power(v0,2) + 4*Power(vf,2)));
    const double h4 = 1728*Power(jMax,6)*(-2*Power(h1,3) - 6*h1*(h2 - 12*Power(aMax,2)*jMax*vf) + 9*Power(aMax,2)*(h2 - 48*Power(jMax,2)*Power(vf,2)));
    const auto h5 = Power(h4 + Sqrt(-11943936*Power(h3,3) + Power(h4,2)),1./3);
    const auto h6 = SqrtComplex((-4*Power(2,1./3)*h3)/(h5*Power(jMax,4)) - h5/(36.*Power(2,1./3)*Power(jMax,4)) + Power(aMax,2)/Power(jMax,2) - (2*h1)/(3.*Power(jMax,2)));
    const auto h7 = SqrtComplex((288*Power(2,1./3)*h3*h6 + h5*(Power(2,2./3)*h5*h6 + 48*jMax*(3*Power(aMax,3) - 3*aMax*h1 + 3*Power(aMax,2)*h6*jMax - 2*h1*h6*jMax + 12*aMax*jMax*vf)))/(h5*h6))/(6.*Sqrt(2)*Power(Abs(jMax),2));
    const auto h8 = SqrtComplex((288*Power(2,1./3)*h3*h6 + h5*(Power(2,2./3)*h5*h6 + 48*jMax*(-3*Power(aMax,3) + 3*Power(aMax,2)*h6*jMax - 2*h1*h6*jMax + 3*aMax*(h1 - 4*jMax*vf))))/(h5*h6))/(6.*Sqrt(2)*Power(Abs(jMax),2));

    // std::cout << "HERE" << std::endl;

    // Solution 2
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = -(-Power(a0,2) + Power(aMax,2) + jMax*(h6*h7*jMax + 2*v0) + aMax*(-h6*jMax + h7*jMax - (2*vf)/h6)).real()/(2.*aMax*jMax);
        profile.t[2] = aMax/jMax;
        profile.t[3] = 0;
        profile.t[4] = (-aMax/jMax - h6 + h7).real()/2;
        profile.t[5] = 0;
        profile.t[6] = profile.t[4];

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    // Solution 4
    {
        profile.t[0] = (-a0 + aMax)/jMax;
        profile.t[1] = -(-Power(a0,2) + Power(aMax,2) + jMax*(-(h6*h8*jMax) + 2*v0) + aMax*(h6*jMax + h8*jMax + (2*vf)/h6)).real()/(2.*aMax*jMax);
        profile.t[2] = aMax/jMax;
        profile.t[3] = 0;
        profile.t[4] = (-aMax/jMax + h6 + h8).real()/2;
        profile.t[5] = 0;
        profile.t[6] = profile.t[4];

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    return false;
}

bool RuckigStep2::time_up_acc0(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    // a3 != 0

    // Solution 1
    profile.t[0] = (-a0 + aMax)/jMax;
    profile.t[1] = -(4*Power(a0,3)*jMax - 12*Power(a0,2)*aMax*jMax + 18*a0*Power(aMax,2)*jMax - 12*Power(aMax,3)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf - 12*a0*aMax*Power(jMax,2)*tf + 18*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) - 12*a0*Power(jMax,2)*v0 + 12*aMax*Power(jMax,2)*v0 + 12*a0*Power(jMax,2)*vf - 12*aMax*Power(jMax,2)*vf + 12*Power(jMax,3)*tf*vf + Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) - 6*Power(aMax,3) + 9*Power(aMax,2)*jMax*tf + 3*a0*aMax*(3*aMax - 2*jMax*tf) + Power(a0,2)*(-6*aMax + 3*jMax*tf) - 6*Power(jMax,2)*(p0 - pf + tf*v0) - 3*aMax*jMax*(jMax*Power(tf,2) - 2*v0 + 2*vf),2) - 9*Power(Power(a0,2) - 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)),3))))/(6.*Power(jMax,2)*(-Power(a0,2) + 2*a0*aMax - 2*Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf));
    profile.t[2] = (2*Power(a0,3)*jMax - 12*Power(a0,2)*aMax*jMax + 18*a0*Power(aMax,2)*jMax - 12*Power(aMax,3)*jMax - 12*Power(jMax,3)*p0 + 12*Power(jMax,3)*pf + 6*Power(a0,2)*Power(jMax,2)*tf - 12*a0*aMax*Power(jMax,2)*tf + 18*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) + 12*aMax*Power(jMax,2)*v0 - 12*Power(jMax,3)*tf*v0 - 12*aMax*Power(jMax,2)*vf + Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) - 6*Power(aMax,3) + 9*Power(aMax,2)*jMax*tf + 3*a0*aMax*(3*aMax - 2*jMax*tf) + Power(a0,2)*(-6*aMax + 3*jMax*tf) - 6*Power(jMax,2)*(p0 - pf + tf*v0) - 3*aMax*jMax*(jMax*Power(tf,2) - 2*v0 + 2*vf),2) - 9*Power(Power(a0,2) - 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)),3))))/(6.*Power(jMax,2)*(-Power(a0,2) + 2*a0*aMax - 2*Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf));
    profile.t[3] = (-2*Power(a0,3)*jMax + 12*Power(a0,2)*aMax*jMax - 18*a0*Power(aMax,2)*jMax + 12*Power(aMax,3)*jMax + 12*Power(jMax,3)*p0 - 12*Power(jMax,3)*pf - 6*Power(a0,2)*Power(jMax,2)*tf + 12*a0*aMax*Power(jMax,2)*tf - 18*Power(aMax,2)*Power(jMax,2)*tf + 6*aMax*Power(jMax,3)*Power(tf,2) - 12*aMax*Power(jMax,2)*v0 + 12*Power(jMax,3)*tf*v0 + 12*aMax*Power(jMax,2)*vf + Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) - 6*Power(aMax,3) + 9*Power(aMax,2)*jMax*tf + 3*a0*aMax*(3*aMax - 2*jMax*tf) + Power(a0,2)*(-6*aMax + 3*jMax*tf) - 6*Power(jMax,2)*(p0 - pf + tf*v0) - 3*aMax*jMax*(jMax*Power(tf,2) - 2*v0 + 2*vf),2) - 9*Power(Power(a0,2) - 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)),3))))/(6.*Power(jMax,2)*(-Power(a0,2) + 2*a0*aMax - 2*Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf));
    profile.t[4] = -(2*Power(a0,3)*jMax - 6*Power(a0,2)*aMax*jMax + 6*a0*Power(aMax,2)*jMax - 12*Power(jMax,3)*p0 + 12*Power(jMax,3)*pf + 6*Power(a0,2)*Power(jMax,2)*tf - 12*a0*aMax*Power(jMax,2)*tf + 6*Power(aMax,2)*Power(jMax,2)*tf - 6*aMax*Power(jMax,3)*Power(tf,2) - 12*Power(jMax,3)*tf*v0 + Sqrt(2)*Sqrt(Power(jMax,2)*(2*Power(Power(a0,3) - 6*Power(aMax,3) + 9*Power(aMax,2)*jMax*tf + 3*a0*aMax*(3*aMax - 2*jMax*tf) + Power(a0,2)*(-6*aMax + 3*jMax*tf) - 6*Power(jMax,2)*(p0 - pf + tf*v0) - 3*aMax*jMax*(jMax*Power(tf,2) - 2*v0 + 2*vf),2) - 9*Power(Power(a0,2) - 2*a0*aMax + 2*(Power(aMax,2) - aMax*jMax*tf + jMax*(-v0 + vf)),3))))/(6.*Power(jMax,2)*(-Power(a0,2) + 2*a0*aMax - 2*Power(aMax,2) + 2*aMax*jMax*tf + 2*jMax*v0 - 2*jMax*vf));
    profile.t[5] = 0;
    profile.t[6] = 0;

    profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
    if (profile.check(pf, vf, vMax, aMax)) {
        return true;
    }

    return false;
}

bool RuckigStep1::time_up_none(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    if (std::abs(v0) < 1e-14 && std::abs(a0) < 1e-14 && std::abs(vf) < 1e-14) {
        profile.t[0] = Power((pf - p0)/(2*jMax),1./3);
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        return profile.check(pf, vf, vMax, aMax);
    }

    if (std::abs(a0) < 1e-14 && std::abs(p0 - pf) < 1e-14) {
        // Solution 2

        profile.t[0] = Sqrt((1 + Sqrt(5))/2.)*Sqrt(-(v0/jMax));
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = ((-1 + Sqrt(5))*Sqrt((1 + Sqrt(5))/2.)*Sqrt(-(v0/jMax)))/2;
        profile.t[5] = 0;
        profile.t[6] = profile.t[4];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        return profile.check(pf, vf, vMax, aMax);
    }

    if (std::abs(a0) < 1e-14 && std::abs(vf) < 1e-14 && std::abs(p0 - pf) < 2e-5) {
        // Solution 2, plus epsilon in timing

        double h4 = Sqrt(1 + Sqrt(5));
        double h1 = (16*(27*jMax*Power(p0 - pf,2) + 26*Power(v0,3)))/Power(jMax,3);
        auto h2 = PowerComplex(h1 + SqrtComplex(Power(h1,2) - (256*Power(v0,6))/Power(jMax,6)),1./3);
        double h3 = Sqrt(-(v0/jMax));
        double h5 = (2*(-4*(4 + 3*Sqrt(5)) + Power(h4,6))*v0)/jMax;
        double h6 = (8.*(-4*p0 + 4*pf + Sqrt(2)*(5 + Sqrt(5))*h3*h4*v0))/jMax;
        auto h7 = (Power(2,2./3)*h2 - (16*v0)/jMax + (8*Power(2,1./3)*Power(v0,2))/(h2*Power(jMax,2)))/6.;
        auto h8 = (-(Sqrt(2)*h3*Power(h4,3)) + 2.*SqrtComplex(h7) + SqrtComplex((Sqrt(2)*h3*Power(h4,3)*h5 + h6)/SqrtComplex(h7) - 4.*h7 + (16*(4 + 3*Sqrt(5))*v0)/jMax - (6*(1 + Sqrt(5))*Power(h4,4)*v0)/jMax))/4.;

        auto eps1 = -(3*Sqrt(2)*(1 + Sqrt(5))*h3*h4*PowerComplex(h8,2)*jMax + 4.*(PowerComplex(h8,3)*jMax + 4*(p0 - pf) - 3.*h8*v0*(1 + Sqrt(5))))/(8.*v0);
        auto eps2 = h8;
        auto eps4 = (3*Sqrt(2)*(1 + Sqrt(5))*h3*h4*PowerComplex(h8,2)*jMax + 4.*PowerComplex(h8,3)*jMax + 16*(p0 - pf) - 4.*h8*v0*(1 + 3*Sqrt(5)))/(8.*v0);

        profile.t[0] = Sqrt((1 + Sqrt(5))/2.)*h3 + eps1.real();
        profile.t[1] = 0;
        profile.t[2] = profile.t[0] - eps1.real() + eps2.real();
        profile.t[3] = 0;
        profile.t[4] = ((-1 + Sqrt(5))*Sqrt((1 + Sqrt(5))/2.)*h3)/2;
        profile.t[5] = 0;
        profile.t[6] = profile.t[4] + eps4.real();

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        return profile.check(pf, vf, vMax, aMax);
    }

    if (std::abs(v0) < 1e-14 && std::abs(vf) < 1e-14) {
        // Solution 2
        {
            const double h1 = Power(a0,3) + 3*Power(jMax,2)*(p0 - pf);
            const double h2 = -Power(a0,8) + 192*Power(a0,5)*Power(jMax,2)*(p0 - pf) + 288*Power(a0,2)*Power(jMax,4)*Power(p0 - pf,2);
            const double h3 = Power(a0,2)*jMax*(Power(a0,3) + 3*Power(jMax,2)*(p0 - pf));
            const double h4 = 17*Power(a0,6) + 48*Power(a0,3)*Power(jMax,2)*(p0 - pf) + 72*Power(jMax,4)*Power(p0 - pf,2);
            const double h5 = 3*(-576*Power(a0,2)*Power(h3,2) + 96*Power(a0,4)*h1*h3*jMax + 3*Power(a0,12)*Power(jMax,2) + (12*Power(a0,6) + 16*Power(h1,2))*h4*Power(jMax,2));
            const double h6 = 6*Power(3*Power(jMax,4)*(h5 + Sqrt(Power(h5,2) - 3*Power(h2,3)*Power(jMax,4))),1./3);
            const double h7 = (h2/h6 + h6/(108.*Power(jMax,4)))/Power(a0,2);
            const double h8 = Sqrt(-9*h7 + (3*Power(a0,6) + 4*Power(h1,2))/(Power(a0,4)*Power(jMax,2)))/3.;
            const double h9 = (8*h1*(-27 + (8*Power(h1,2))/Power(a0,6)))/(27.*Power(jMax,3));
            const double h10 = (-6*h8 + Sqrt(36*h7 - (9*h9)/h8 + (8*(3*Power(a0,6) + 4*Power(h1,2)))/(Power(a0,4)*Power(jMax,2))) + (4*h1)/(Power(a0,2)*jMax))/12.;

            profile.t[0] = (-6*h8 + Sqrt(36*h7 - (9*h9)/h8 + (8*(3*Power(a0,6) + 4*Power(h1,2)))/(Power(a0,4)*Power(jMax,2))) - (8*a0)/jMax + (12*jMax*(p0 - pf))/Power(a0,2))/12.;
            profile.t[1] = 0;
            profile.t[2] = h10;
            profile.t[3] = 0;
            profile.t[4] = (-12*Power(a0,7) + 17*Power(a0,6)*h10*jMax + 12*Power(a0,5)*Power(h10,2)*Power(jMax,2) - 18*Power(a0,4)*Power(jMax,2)*(Power(h10,3)*jMax + 2*p0 - 2*pf) + 48*Power(a0,3)*h10*Power(jMax,3)*(p0 - pf) + 36*Power(a0,2)*Power(h10,2)*Power(jMax,4)*(p0 - pf) + 72*h10*Power(jMax,5)*Power(p0 - pf,2))/(-(Power(a0,6)*jMax) + 48*Power(a0,3)*Power(jMax,3)*(p0 - pf) + 72*Power(jMax,5)*Power(p0 - pf,2));
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }
    }

    if (std::abs(a0) < 1e-14 && std::abs(vf) < 1e-14) {
        // Solution 2
        {
            const double h1 = jMax*Power(p0 - pf,2) - Power(v0,3);
            const double h2 = p0 - pf;
            const double h3 = 2*Power(jMax,3)*(54*h1*Power(h2,2)*jMax - 36*h1*Power(v0,3) + 180*Power(h2,2)*jMax*Power(v0,3) + Power(v0,6));
            const double h4 = -(Power(jMax,2)*v0*(12*jMax*Power(p0 - pf,2) + 11*Power(v0,3)));
            const double h5 = Power(h3 + Sqrt(Power(h3,2) - 4*Power(h4,3)),1./3);
            // Be careful of numerical stability of h6
            const long double h6 = (2*PowerLong(2,1./3)*h4 + PowerLong(2,2./3)*PowerLong(h5,2))/(6.*h5*PowerLong(jMax,2)*v0);
            const long double h7 = PowerLong(h2,2)/PowerLong(v0,2) - (2*v0)/(3.*jMax);
            const long double h8 = (-24*h2)/jMax - (8*PowerLong(h2,3))/PowerLong(v0,3);
            const auto h9 = SqrtComplexLong(-h6 + 2*h7 - h8/(4.l*SqrtComplexLong(h6 + h7)));
            const auto h10_c = -((long double)h2/v0 + SqrtComplexLong(h6 + h7) - h9)/(2.l);
            const double h10 = h10_c.real();

            profile.t[0] = h10;
            profile.t[1] = 0;
            profile.t[2] = (-p0 + pf + (-SqrtComplexLong(h6 + h7) + h9).real()*v0)/(2*v0);
            profile.t[3] = 0;
            profile.t[4] = (Power(h10,2)*jMax*(-p0 + pf)*v0 - Power(h10,3)*jMax*Power(v0,2) + 2*(-p0 + pf)*Power(v0,2) + h10*(jMax*Power(p0 - pf,2) - Power(v0,3)))/(jMax*Power(p0 - pf,2) + Power(v0,3));
            profile.t[5] = 0;
            profile.t[6] = profile.t[4];

            profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
            if (profile.check(pf, vf, vMax, aMax)) {
                return true;
            }
        }
    }

    double h1 = 2*Power(a0,3) - 3*Power(jMax,2)*(p0 - pf) - 3.*a0*jMax*(v0 - 2*vf);
    long double h2 = PowerLong(a0,2) - 2.l*jMax*(v0 - vf);
    long double h3 = PowerLong(a0,5) - 24*PowerLong(a0,2)*PowerLong(jMax,2)*(p0 - pf) - 24*PowerLong(jMax,3)*(p0 - pf)*v0 + 4*PowerLong(a0,3)*jMax*(v0 + 3*vf) + 12*a0*PowerLong(jMax,2)*(PowerLong(v0,2) + 2*v0*vf - PowerLong(vf,2));
    long double h4 = 3*PowerLong(a0,4) - 24*a0*PowerLong(jMax,2)*(p0 - pf) - 4*PowerLong(jMax,2)*PowerLong(v0 - vf,2) + 4*PowerLong(a0,2)*jMax*(v0 + 5*vf);
    long double h5 = PowerLong(a0,6) - 48*PowerLong(a0,3)*PowerLong(jMax,2)*(p0 - pf) - 144*a0*PowerLong(jMax,3)*(p0 - pf)*v0 + 6*PowerLong(a0,4)*jMax*(v0 + 3*vf) + 36*PowerLong(a0,2)*PowerLong(jMax,2)*(PowerLong(v0,2) + 2*v0*vf - PowerLong(vf,2)) - 72*PowerLong(jMax,3)*(jMax*PowerLong(p0 - pf,2) - (v0 - vf)*PowerLong(v0 + vf,2));
    long double h17 = jMax*(-PowerLong(a0,6) + 48*PowerLong(a0,3)*PowerLong(jMax,2)*(p0 - pf) - 144*a0*PowerLong(jMax,3)*(p0 - pf)*v0 + 6*PowerLong(a0,4)*jMax*(v0 - 3*vf) - 36*PowerLong(a0,2)*PowerLong(jMax,2)*(PowerLong(v0,2) - 2*v0*vf - PowerLong(vf,2)) + 72*PowerLong(jMax,3)*(jMax*PowerLong(p0 - pf,2) + PowerLong(v0 - vf,2)*(v0 + vf)));
    long double h6 = -PowerLong(a0,8) + 192*PowerLong(a0,5)*PowerLong(jMax,2)*(p0 - pf) + 8*PowerLong(a0,6)*jMax*(v0 - 5*vf) + 1152*a0*PowerLong(jMax,4)*(p0 - pf)*v0*(v0 + vf) - 192*PowerLong(a0,3)*PowerLong(jMax,3)*(p0 - pf)*(5*v0 + 2*vf) - 120*PowerLong(a0,4)*PowerLong(jMax,2)*(PowerLong(v0,2) - 2*v0*vf - 3*PowerLong(vf,2)) + 96*PowerLong(a0,2)*PowerLong(jMax,3)*(3*jMax*PowerLong(p0 - pf,2) + 5*PowerLong(v0,3) - 3*PowerLong(v0,2)*vf - 15*v0*PowerLong(vf,2) + PowerLong(vf,3)) - 48*PowerLong(jMax,4)*(12*jMax*PowerLong(p0 - pf,2)*(v0 + vf) + PowerLong(v0 - vf,2)*(11*PowerLong(v0,2) + 26.l*v0*vf + 11*PowerLong(vf,2)));
    long double h8 = 4*PowerLong(h1,2)/(9*h2) - h4/3;
    long double h9 = -2*(h3 - h1*h4/h2 + PowerLong(2*h1,3)/PowerLong(3*h2,2))/3;
    long double h7 = 3*(36*h2*PowerLong(h3,2) + 4*PowerLong(2*h1,2)*h5 + 3*h4*(PowerLong(h4,2) - 8*h1*h3 - 4*h2*h5));

    // Important: Numerical stability of h10
    const auto h10_x = PowerLong(h6,3)/PowerLong(h7,2);
    std::complex<long double> h10 = PowerComplexLong(3 * h7 * (1.l - SqrtComplexLong(1.l - 3.l*h10_x)),1./3);
    if (std::abs(h10_x) < 1e-7) {
        h10 = PowerComplexLong(h10_x,1./3)*PowerComplexLong(9*h7/2,1./3) + PowerComplexLong(h10_x,4./3)*PowerComplexLong(9*h7/2,1./3)/4.l + 5.l*PowerComplexLong(h10_x,7./3)*PowerComplexLong(9*h7/2,1./3)/16.l;
    }

    const auto h11 = (h6/h10 + h10/3.l)/6.l;
    const auto h12 = SqrtComplexLong((h2*(h10 - 3.l*(2*h4 - h6/h10)) + 2*PowerLong(2*h1,2))/18.l)/h2;

    // std::cout << std::setprecision(15) << "---" << std::endl;
    // std::cout << "h1-6: " << h1 << " " << h2 << " " << h3 << " " << h4 << " " << h5 << " " << h6 << std::endl;
    // std::cout << "h7-12: " << h7 << " " << h8 << " " << h9 << " " << h10 << " " << h11 << " " << h12 << std::endl;

    const auto h12_a = SqrtComplexLong((2*h8 - h11 + h9/h12)/h2);
    const auto h12_b = SqrtComplexLong((2*h8 - h11 - h9/h12)/h2);

    const auto h13_c = (h12 - h12_a)/(2.l*jMax) - h1/(3*h2*jMax);
    const auto h14_c = (h12 + h12_a)/(2.l*jMax) - h1/(3*h2*jMax);
    const auto h15_c = (-h12 + h12_b)/(2.l*jMax) - h1/(3*h2*jMax);
    const auto h16_c = (-h12 - h12_b)/(2.l*jMax) - h1/(3*h2*jMax);
    // std::cout << "h13-16: " << h13_c << " " << h14_c << " " << h15_c << " " << h16_c << std::endl;

    // Solution 3
    if (h13_c.real() > 0.0 && std::abs(h13_c.imag()) < 2e-8) {
        const double h13 = std::abs(h13_c);

        profile.t[0] = h13;
        profile.t[1] = 0;
        profile.t[2] = (-4*PowerLong(a0,3) + 3.l*((h2*(h12 - h12_a)) + 2*PowerLong(jMax,2)*(p0 - pf)) + 6l*a0*(h2 + jMax*(v0 - 2*vf))).real()/(6.l*h2*jMax);
        profile.t[3] = 0;
        profile.t[4] = -((Power(a0,7) + 13*Power(a0,6)*h13*jMax + 72*Power(jMax,4)*(-(h13*(jMax*Power(p0 - pf,2) - Power(v0 - vf,3))) + Power(h13,2)*jMax*(p0 - pf)*(v0 - vf) + 2*(p0 - pf)*v0*(v0 - vf) + Power(h13,3)*jMax*Power(v0 - vf,2)) + 6*Power(a0,5)*jMax*(7*Power(h13,2)*jMax + v0 + 3*vf) - 12*Power(a0,3)*Power(jMax,2)*(10*h13*jMax*(p0 - pf) - Power(v0,2) + Power(h13,2)*jMax*(13*v0 - 16*vf) - 2*v0*vf + 3*Power(vf,2)) + 6*Power(a0,4)*Power(jMax,2)*(3*Power(h13,3)*jMax - 8*p0 + 8*pf + h13*(v0 + 19*vf)) - 36*Power(a0,2)*Power(jMax,3)*(Power(h13,2)*jMax*(p0 - pf) + 2*(-p0 + pf)*v0 + 2*Power(h13,3)*jMax*(v0 - vf) + h13*(3*Power(v0,2) + 2*v0*vf - 3*Power(vf,2))) - 72*a0*Power(jMax,3)*(Power(v0,3) + Power(v0,2)*vf - 3*v0*Power(vf,2) + Power(vf,3) + jMax*(Power(p0,2) + Power(pf,2) + h13*pf*(4*v0 - 2*vf) - 2*p0*(pf + 2*h13*v0 - h13*vf) + Power(h13,2)*(-2*Power(v0,2) + 5*v0*vf - 3*Power(vf,2)))))/h17);
        profile.t[5] = 0;
        profile.t[6] = profile.t[4];

        // Set average as only sum of t2 and t4 needs to be >0
        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    // Solution 4
    // 2e-18 comes from sqrt(double precision)
    if (h14_c.real() > 0.0 && std::abs(h14_c.imag()) < 2e-8) {
        const double h14 = h14_c.real();

        profile.t[0] = h14;
        profile.t[1] = 0;
        profile.t[2] = (-4*PowerLong(a0,3) + 3.l*(h2*(h12 + h12_a) + 2*PowerLong(jMax,2)*(p0 - pf)) + 6l*a0*(h2 + jMax*(v0 - 2*vf))).real()/(6.l*h2*jMax);
        profile.t[3] = 0;
        profile.t[4] = -((Power(a0,7) + 13*Power(a0,6)*h14*jMax + 72*Power(jMax,4)*(-(h14*(jMax*Power(p0 - pf,2) - Power(v0 - vf,3))) + Power(h14,2)*jMax*(p0 - pf)*(v0 - vf) + 2*(p0 - pf)*v0*(v0 - vf) + Power(h14,3)*jMax*Power(v0 - vf,2)) + 6*Power(a0,5)*jMax*(7*Power(h14,2)*jMax + v0 + 3*vf) - 12*Power(a0,3)*Power(jMax,2)*(10*h14*jMax*(p0 - pf) - Power(v0,2) + Power(h14,2)*jMax*(13*v0 - 16*vf) - 2*v0*vf + 3*Power(vf,2)) + 6*Power(a0,4)*Power(jMax,2)*(3*Power(h14,3)*jMax - 8*p0 + 8*pf + h14*(v0 + 19*vf)) - 36*Power(a0,2)*Power(jMax,3)*(Power(h14,2)*jMax*(p0 - pf) + 2*(-p0 + pf)*v0 + 2*Power(h14,3)*jMax*(v0 - vf) + h14*(3*Power(v0,2) + 2*v0*vf - 3*Power(vf,2))) - 72*a0*Power(jMax,3)*(Power(v0,3) + Power(v0,2)*vf - 3*v0*Power(vf,2) + Power(vf,3) + jMax*(Power(p0,2) + Power(pf,2) + h14*pf*(4*v0 - 2*vf) - 2*p0*(pf + 2*h14*v0 - h14*vf) + Power(h14,2)*(-2*Power(v0,2) + 5*v0*vf - 3*Power(vf,2)))))/h17);
        profile.t[5] = 0;
        profile.t[6] = profile.t[4];

        // std::cout << std::setprecision(9) << profile.t[0] << std::endl;
        // std::cout << profile.t[1] << std::endl;
        // std::cout << profile.t[2] << std::endl;
        // std::cout << profile.t[3] << std::endl;
        // std::cout << profile.t[4] << std::endl;
        // std::cout << profile.t[5] << std::endl;
        // std::cout << profile.t[6] << std::endl << std::endl;

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    // Solution 2
    if (h15_c.real() > 0.0 && std::abs(h15_c.imag()) < 2e-8) {
        const double h15 = h15_c.real();

        profile.t[0] = h15;
        profile.t[1] = 0;
        profile.t[2] = (-4*PowerLong(a0,3) + 3.l*(h2*(-h12 + h12_b) + 2*PowerLong(jMax,2)*(p0 - pf)) + 6l*a0*(h2 + jMax*(v0 - 2*vf))).real()/(6.l*h2*jMax);
        profile.t[3] = 0;
        profile.t[4] = -((Power(a0,7) + 13*Power(a0,6)*h15*jMax + 72*Power(jMax,4)*(-(h15*(jMax*Power(p0 - pf,2) - Power(v0 - vf,3))) + Power(h15,2)*jMax*(p0 - pf)*(v0 - vf) + 2*(p0 - pf)*v0*(v0 - vf) + Power(h15,3)*jMax*Power(v0 - vf,2)) + 6*Power(a0,5)*jMax*(7*Power(h15,2)*jMax + v0 + 3*vf) - 12*Power(a0,3)*Power(jMax,2)*(10*h15*jMax*(p0 - pf) - Power(v0,2) + Power(h15,2)*jMax*(13*v0 - 16*vf) - 2*v0*vf + 3*Power(vf,2)) + 6*Power(a0,4)*Power(jMax,2)*(3*Power(h15,3)*jMax - 8*p0 + 8*pf + h15*(v0 + 19*vf)) - 36*Power(a0,2)*Power(jMax,3)*(Power(h15,2)*jMax*(p0 - pf) + 2*(-p0 + pf)*v0 + 2*Power(h15,3)*jMax*(v0 - vf) + h15*(3*Power(v0,2) + 2*v0*vf - 3*Power(vf,2))) - 72*a0*Power(jMax,3)*(Power(v0,3) + Power(v0,2)*vf - 3*v0*Power(vf,2) + Power(vf,3) + jMax*(Power(p0,2) + Power(pf,2) + h15*pf*(4*v0 - 2*vf) - 2*p0*(pf + 2*h15*v0 - h15*vf) + Power(h15,2)*(-2*Power(v0,2) + 5*v0*vf - 3*Power(vf,2)))))/h17);
        profile.t[5] = 0;
        profile.t[6] = profile.t[4];

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }

    // Solution 1
    if (h16_c.real() > 0.0 && std::abs(h16_c.imag()) < 2e-8) {
        const double h16 = h16_c.real();

        profile.t[0] = h16;
        profile.t[1] = 0;
        profile.t[2] = (-4*PowerLong(a0,3) + 3.l*(h2*(-h12 - h12_b) + 2*PowerLong(jMax,2)*(p0 - pf)) + 6l*a0*(h2 + jMax*(v0 - 2*vf))).real()/(6.l*h2*jMax);
        profile.t[3] = 0;
        profile.t[4] = -((PowerLong(a0,7) + 13*PowerLong(a0,6)*h16*jMax + 72*PowerLong(jMax,4)*(-(h16*(jMax*PowerLong(p0 - pf,2) - PowerLong(v0 - vf,3))) + PowerLong(h16,2)*jMax*(p0 - pf)*(v0 - vf) + 2*(p0 - pf)*v0*(v0 - vf) + PowerLong(h16,3)*jMax*PowerLong(v0 - vf,2)) + 6*PowerLong(a0,5)*jMax*(7*PowerLong(h16,2)*jMax + v0 + 3*vf) - 12*PowerLong(a0,3)*PowerLong(jMax,2)*(10*h16*jMax*(p0 - pf) - PowerLong(v0,2) + PowerLong(h16,2)*jMax*(13*v0 - 16*vf) - 2*v0*vf + 3*Power(vf,2)) + 6*PowerLong(a0,4)*PowerLong(jMax,2)*(3*PowerLong(h16,3)*jMax - 8*p0 + 8*pf + h16*(v0 + 19*vf)) - 36*PowerLong(a0,2)*PowerLong(jMax,3)*(PowerLong(h16,2)*jMax*(p0 - pf) + 2*(-p0 + pf)*v0 + 2*PowerLong(h16,3)*jMax*(v0 - vf) + h16*(3*PowerLong(v0,2) + 2*v0*vf - 3*PowerLong(vf,2))) - 72*a0*PowerLong(jMax,3)*(PowerLong(v0,3) + PowerLong(v0,2)*vf - 3*v0*PowerLong(vf,2) + PowerLong(vf,3) + jMax*(PowerLong(p0,2) + Power(pf,2) + h16*pf*(4*v0 - 2*vf) - 2*p0*(pf + 2*h16*v0 - h16*vf) + PowerLong(h16,2)*(-2*PowerLong(v0,2) + 5*v0*vf - 3*PowerLong(vf,2)))))/h17);
        profile.t[5] = 0;
        profile.t[6] = profile.t[4];

        profile.t[2] = (profile.t[2] + profile.t[4]) / 2;
        profile.t[4] = profile.t[2];

        profile.set(p0, v0, a0, {jMax, 0, -jMax, 0, -jMax, 0, jMax});
        if (profile.check(pf, vf, vMax, aMax)) {
            return true;
        }
    }
    return false;
}

bool RuckigStep2::time_up_none(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    if (std::abs(v0) < 1e-14 && std::abs(a0) < 1e-14 && std::abs(vf) < 1e-14) {
        profile.t[0] = tf/4;
        profile.t[1] = 0;
        profile.t[2] = profile.t[0];
        profile.t[3] = 0;
        profile.t[4] = profile.t[0];
        profile.t[5] = 0;
        profile.t[6] = profile.t[0];

        double jMaxNew = (-32*(p0 - pf))/Power(tf,3);

        profile.set(p0, v0, a0, {jMaxNew, 0, -jMaxNew, 0, -jMaxNew, 0, jMaxNew});
        return profile.check(pf, vf, vMax, aMax);
    }
    
    if (std::abs(v0) < 1e-14 && std::abs(a0) < 1e-14) {
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
        return profile.check(pf, vf, vMax, aMax);
    }

    if (std::abs(a0) < 1e-14 && std::abs(vf) < 1e-14) {
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
            if (profile.check(pf, vf, vMax, aMax)) {
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
                if (profile.check(pf, vf, vMax, aMax)) {
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
                if (profile.check(pf, vf, vMax, aMax)) {
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
                if (profile.check(pf, vf, vMax, aMax)) {
                    return true;
                }
            }
        }
    }
    
    double h1 = 8*p0 - 8*pf + a0*Power(tf,2) + 4*tf*v0 + 4*tf*vf;
    double h2 = -16*a0*p0 + 16*a0*pf + Power(a0,2)*Power(tf,2) + 8*Power(v0,2) - 16*a0*tf*vf - 16*v0*vf + 8*Power(vf,2);
    double h3 = 432*Power(h2,3) + 3888*Power(a0,4)*Power(h1,2)*Power(tf,2) + 2592*Power(a0,3)*h1*h2*Power(tf,2) + 1296*Power(a0,4)*h2*Power(tf,4) - 1296*Power(a0,6)*Power(tf,6);
    double h4 = 36*Power(h2,2) + 144*Power(a0,3)*h1*Power(tf,2) - 36*Power(a0,4)*Power(tf,4);
    auto h5 = PowerComplex(h3 + SqrtComplex(Power(h3,2) - 4*Power(h4,3)),1./3);
    auto h6 = SqrtComplex((4*Power(h1,2))/Power(tf,6) + (4*h2)/(3.*Power(tf,4)) - (Power(2,1./3)*h4)/(9.*h5*Power(tf,4)) - h5/(9.*Power(2,1./3)*Power(tf,4)));
    auto h7 = (8*Power(h1,2))/Power(tf,6) + (8*h2)/(3.*Power(tf,4)) + (Power(2,1./3)*h4)/(9.*h5*Power(tf,4)) + h5/(9.*Power(2,1./3)*Power(tf,4));
    auto h8 = ((-64*Power(h1,3))/Power(tf,9) - (32*h1*h2)/Power(tf,7) + (32*Power(a0,3))/(3.*Power(tf,3)))/(4.*h6);
    auto h9 = -h6/2. - SqrtComplex(h7 + h8)/2. - h1/Power(tf,3);
    auto h10 = -h6/2. + SqrtComplex(h7 + h8)/2. - h1/Power(tf,3);
    auto h11 = h6/2. - SqrtComplex(h7 + h8)/2. - h1/Power(tf,3);
    auto h12 = h6/2. + SqrtComplex(h7 + h8)/2. - h1/Power(tf,3);

    // std::cout << h9 << " " << h10 << " " << h11 << " " << h12 << std::endl;

    return false;
}

bool RuckigStep1::time_down_acc0_acc1_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc0_acc1_vel(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1_vel(profile, tf, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc1_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc1_vel(profile, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc1_vel(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc1_vel(profile, tf, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc0_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc0_vel(profile, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc0_vel(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc0_vel(profile, tf, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_vel(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_vel(profile, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_vel(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_vel(profile, tf, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc0_acc1(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc0_acc1(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc0_acc1(profile, tf, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc1(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc1(profile, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc1(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc1(profile, tf, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_acc0(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc0(profile, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_acc0(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_acc0(profile, tf, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep1::time_down_none(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_none(profile, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep2::time_down_none(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    return time_up_none(profile, tf, p0, v0, a0, pf, vf, -vMax, -aMax, -jMax);
}

bool RuckigStep1::get_profile(Profile& profile, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    // Test all cases to get ones that match
    if (pf > p0) {
        if (time_up_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1_VEL;

        } else if (time_down_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1_VEL;

        } else if (time_up_acc1_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1_VEL;

        } else if (time_down_acc1_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1_VEL;

        } else if (time_up_acc0_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_VEL;

        } else if (time_down_acc0_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_VEL;

        } else if (time_up_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_VEL;

        } else if (time_down_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_VEL;

        } else if (time_up_none(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_NONE;

        } else if (time_down_none(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_NONE;

        } else if (time_up_acc1(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1;

        } else if (time_down_acc1(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1;

        } else if (time_up_acc0(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0;

        } else if (time_down_acc0(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0;

        } else if (time_up_acc0_acc1(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1;

        } else if (time_down_acc0_acc1(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1;

        } else {
            return false;
        }

    } else {
        if (time_down_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1_VEL;

        } else if (time_up_acc0_acc1_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1_VEL;

        } else if (time_down_acc1_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1_VEL;

        } else if (time_up_acc1_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1_VEL;

        } else if (time_down_acc0_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_VEL;

        } else if (time_up_acc0_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_VEL;

        } else if (time_down_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_VEL;

        } else if (time_up_vel(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_VEL;

        } else if (time_down_none(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_NONE;

        } else if (time_up_none(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_NONE;

        } else if (time_down_acc1(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1;

        } else if (time_up_acc1(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1;

        } else if (time_down_acc0(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0;

        } else if (time_up_acc0(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0;

        } else if (time_down_acc0_acc1(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1;

        } else if (time_up_acc0_acc1(profile, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1;

        } else {
            return false;
        }
    }
    return true;
}

bool RuckigStep2::get_profile(Profile& profile, double tf, double p0, double v0, double a0, double pf, double vf, double vMax, double aMax, double jMax) {
    // Test all cases to get ones that match
    if (pf > p0) {
        if (time_up_none(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_NONE;

        } else if (time_down_none(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_NONE;

        } else if (time_up_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_VEL;

        } else if (time_down_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_VEL;

        } else if (time_up_acc0_acc1_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1_VEL;

        } else if (time_down_acc0_acc1_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1_VEL;

        } else if (time_up_acc0_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_VEL;

        } else if (time_down_acc0_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_VEL;

        } else if (time_up_acc1_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1_VEL;

        } else if (time_down_acc1_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1_VEL;

        } else if (time_up_acc0_acc1(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1;

        } else if (time_down_acc0_acc1(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1;

        } else if (time_up_acc1(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1;

        } else if (time_down_acc1(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1;

        } else if (time_up_acc0(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0;

        } else if (time_down_acc0(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0;

        } else {
            return false;
        }

    } else {
        if (time_down_none(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_NONE;

        } else if (time_up_none(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_NONE;

        } else if (time_down_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_VEL;

        } else if (time_up_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_VEL;

        } else if (time_down_acc0_acc1_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1_VEL;

        } else if (time_up_acc0_acc1_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1_VEL;

        } else if (time_down_acc0_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_VEL;

        } else if (time_up_acc0_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_VEL;

        } else if (time_down_acc1_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1_VEL;

        } else if (time_up_acc1_vel(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1_VEL;

        } else if (time_down_acc0_acc1(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0_ACC1;

        } else if (time_up_acc0_acc1(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0_ACC1;

        } else if (time_down_acc1(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC1;

        } else if (time_up_acc1(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC1;

        } else if (time_down_acc0(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::DOWN_ACC0;

        } else if (time_up_acc0(profile, tf, p0, v0, a0, pf, vf, vMax, aMax, jMax)) {
            profile.type = Profile::Type::UP_ACC0;

        } else {
            return false;
        }
    }
    return true;
}

inline double v_at_t(double v0, double a0, double j, double t) {
    return v0 + a0 * t + j * std::pow(t, 2) / 2;
}

void RuckigStep1::get_brake_trajectory(double v0, double a0, double vMax, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake) {
    t_brake[0] = 0.0;
    t_brake[1] = 0.0;
    j_brake[0] = 0.0;
    j_brake[1] = 0.0;

    const double a0_sq_jMax = std::pow(a0, 2)/(2 * jMax) - vMax;

    if (a0 > aMax) {
        j_brake[0] = -jMax;

        double t_to_a_max = (a0 - aMax) / jMax;
        double v_at_a_max = v_at_t(v0, a0, -jMax, t_to_a_max);

        t_brake[0] = t_to_a_max + 2e-15;

        double t_to_a_zero = a0 / jMax;
        double v_at_a_zero = v_at_t(v0, a0, -jMax, t_to_a_zero);

        if (v_at_a_zero > vMax) {
            t_brake[0] = t_to_a_zero + 2e-15;
            v_at_a_max = v_at_a_zero;
        }

        if (v_at_a_max < -vMax) {
            double t_to_v_max_while_a_max = -(v_at_a_max + vMax)/aMax;
            double t_to_v_max_in_reverse_j_direction = (-std::pow(aMax, 2) - 2 * jMax * v_at_a_max + 2 * jMax * vMax)/(2*aMax*jMax);
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
            
        } else if (v_at_a_max > vMax) {
            double t_to_other_a_max = (a0 + aMax) / jMax - 2e-15;
            double t_to_v_max = -(-a0 - std::sqrt(std::pow(a0,2) + 2 * jMax * v0 - 2 * jMax * vMax))/jMax;
            double t_to_v_max_brake_for_other = (2 * a0 + std::sqrt(2) * std::sqrt(std::pow(a0, 2) + 2 * jMax * (v0 + vMax)))/(2 * jMax);

            // std::cout << t_to_other_a_max << " " << t_to_v_max << " " << t_to_v_max_brake_for_other << std::endl;

            if (t_to_v_max_brake_for_other < t_to_other_a_max && t_to_v_max_brake_for_other < t_to_v_max) {
                t_brake[0] = t_to_v_max_brake_for_other - 2e-15;

            } else if (t_to_other_a_max < t_to_v_max) {
                double v_at_a_other_a_max = v_at_t(v0, a0, -jMax, t_to_other_a_max);
                double t_to_v_max_while_a_max = (v_at_a_other_a_max - vMax)/aMax;
                double t_to_v_max_brake_for_other_a_max = -(std::pow(aMax, 2) - 2 * jMax * (v_at_a_max + vMax))/(2 * aMax * jMax);

                t_brake[0] = t_to_other_a_max;
                t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_brake_for_other_a_max);
            
            } else {
                t_brake[0] = t_to_v_max + 2e-15;
            }
        }

    } else if (a0 < -aMax) {
        j_brake[0] = jMax;

        double t_to_a_max = -(a0 + aMax) / jMax;
        double v_at_a_max = v_at_t(v0, a0, jMax, t_to_a_max);

        t_brake[0] = t_to_a_max + 2e-15;

        double t_to_a_zero = -a0 / jMax;
        double v_at_a_zero = v_at_t(v0, a0, jMax, t_to_a_zero);
        if (v_at_a_zero < -vMax) {
            t_brake[0] = t_to_a_zero + 2e-15;
            v_at_a_max = v_at_a_zero;
        }

        if (v_at_a_max > vMax) {
            double t_to_v_max_while_a_max = (v_at_a_max - vMax)/aMax;
            double t_to_v_max_in_reverse_j_direction = (-std::pow(aMax, 2) + 2 * jMax * v_at_a_max + 2 * jMax * vMax)/(2*aMax*jMax);
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
            // {a1 == -aMax, v1 == v0 - aMax * t1, -vMax == v1 + a1 * t2 + jMax / 2 * t2^2, 0 == a1 + jMax t2}, {t1, t2, a1, v1}

        } else if (v_at_a_max < -vMax) {
            double t_to_other_a_max = -(a0 - aMax) / jMax - 2e-15;
            double t_to_v_max = (-a0 + std::sqrt(std::pow(a0,2) - 2 * jMax * (v0 + vMax)))/jMax;

            if (t_to_other_a_max < t_to_v_max) {
                double v_at_a_other_a_max = v_at_t(v0, a0, jMax, t_to_other_a_max);
                double t_to_v_max_while_a_max = -(v_at_a_other_a_max + vMax)/aMax;
                double t_to_v_max_brake_for_other_a_max = -(std::pow(aMax, 2) + 2 * jMax * (v_at_a_max - vMax))/(2 * aMax * jMax);
                
                t_brake[0] = t_to_other_a_max;
                t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_brake_for_other_a_max);
            } else {
                t_brake[0] = t_to_v_max + 2e-15;
            }
        }

    } else if ((v0 > vMax && (a0 > 0 || a0_sq_jMax < v0)) || (a0 > 0 && a0_sq_jMax > -v0)) {
        j_brake[0] = -jMax;
        double t_to_a_max = (a0 + aMax)/jMax;
        double t_to_v_max_in_j_direction = (a0 + std::sqrt(std::pow(a0, 2) - 2 * jMax * (vMax - v0)))/jMax;
        double t_to_v_max_in_reverse_j_direction = (2 * a0 + std::sqrt(2*(std::pow(a0, 2) + 2 * jMax * (v0 + vMax))))/(2 * jMax);

        if (t_to_a_max < t_to_v_max_in_j_direction && t_to_a_max < t_to_v_max_in_reverse_j_direction) {
            double v_at_a_max = v_at_t(v0, a0, -jMax, t_to_a_max);
            double t_to_v_max_while_a_max = (v_at_a_max - vMax)/aMax;
            double t_to_v_max_in_reverse_j_direction = -(std::pow(aMax, 2) - 2 * jMax * (v_at_a_max + vMax))/(2 * aMax * jMax);

            t_brake[0] = t_to_a_max;
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
        } else {
            t_brake[0] = std::min(t_to_v_max_in_j_direction, t_to_v_max_in_reverse_j_direction);
        }

        t_brake[0] = std::max(t_brake[0] - 2e-15, 0.0);

    } else if ((v0 < -vMax && (a0 < 0 || a0_sq_jMax > v0)) || (a0 < 0 && a0_sq_jMax > v0)) {
        j_brake[0] = jMax;
        double t_to_a_max = -(a0 - aMax)/jMax;
        double t_to_v_max_in_j_direction = (-a0 + std::sqrt(std::pow(a0, 2) + 2 * jMax * (-vMax - v0)))/jMax;
        double t_to_v_max_in_reverse_j_direction = (-2 * a0 + std::sqrt(2*(std::pow(a0, 2) - 2 * jMax * (v0 - vMax))))/(2 * jMax);

        if (t_to_a_max < t_to_v_max_in_j_direction && t_to_a_max < t_to_v_max_in_reverse_j_direction) {
            double v_at_a_max = v_at_t(v0, a0, jMax, t_to_a_max);
            double t_to_v_max_while_a_max = -(v_at_a_max + vMax)/aMax;
            double t_to_v_max_in_reverse_j_direction = -(std::pow(aMax, 2) + 2 * jMax * (v_at_a_max - vMax))/(2 * aMax * jMax);

            t_brake[0] = t_to_a_max;
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
        } else {
            t_brake[0] = std::min(t_to_v_max_in_j_direction, t_to_v_max_in_reverse_j_direction);
        }

        t_brake[0] = std::max(t_brake[0] - 2e-15, 0.0);
    }

    // std::cout << t_brake[0] << " " << t_brake[1] << std::endl;
}

} // namespace ruckig
