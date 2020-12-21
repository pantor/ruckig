#include <iomanip>

#include <ruckig/ruckig.hpp>


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
        v[i+1] = v[i] + t[i] * (a[i] + t[i] * j[i] / 2);
        p[i+1] = p[i] + t[i] * (v[i] + t[i] * (a[i] / 2 + t[i] * j[i] / 6));
    }
}

bool Profile::check(double pf, double vf, double af, double vMax, double aMax) const {
    // Velocity and acceleration limits can be broken in the beginnging if the initial velocity and acceleration are too high
    // std::cout << std::setprecision(15) << "target: " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << std::endl;
    return std::all_of(t.begin(), t.end(), [](double tm){ return tm >= 0; })
        && std::all_of(v.begin() + 3, v.end(), [vMax](double vm){ return std::abs(vm) < std::abs(vMax) + 1e-9; })
        && std::all_of(a.begin() + 2, a.end(), [aMax](double am){ return std::abs(am) < std::abs(aMax) + 1e-9; })
        && std::abs(p[7] - pf) < 1e-8 && std::abs(v[7] - vf) < 1e-8 && std::abs(a[7] - af) < 1e-8;
}

std::tuple<double, double, double> Profile::integrate(double t, double p0, double v0, double a0, double j) {
    const double p_new = p0 + t * (v0 + t * (a0 / 2 + t * j / 6));
    const double v_new = v0 + t * (a0 + t * j / 2);
    const double a_new = a0 + t * j;
    return {p_new, v_new, a_new};
}

} // namespace ruckig
