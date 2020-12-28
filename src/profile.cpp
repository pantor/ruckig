#include <iomanip>

#include <ruckig/ruckig.hpp>


namespace ruckig {

std::string Profile::to_string() const {
    std::string result;
    switch (direction) {
        case Direction::UP: result += "UP_"; break;
        case Direction::DOWN: result += "DOWN_"; break;
    }
    switch (limits) {
        case Limits::ACC0_ACC1_VEL: result += "ACC0_ACC1_VEL"; break;
        case Limits::VEL: result += "VEL"; break;
        case Limits::ACC0: result += "ACC0"; break;
        case Limits::ACC1: result += "ACC1"; break;
        case Limits::ACC0_ACC1: result += "ACC0_ACC1"; break;
        case Limits::ACC0_VEL: result += "ACC0_VEL"; break;
        case Limits::ACC1_VEL: result += "ACC1_VEL"; break;
        case Limits::NONE: result += "NONE"; break;
    }
    return result; 
}

bool Profile::check(double pf, double vf, double af, double vMax, double aMax, const std::array<double, 7>& j) {
    this->j = j;
    t_sum[0] = t[0];
    if (t[0] < 0) {
        return false;
    }

    for (size_t i = 0; i < 6; i += 1) {
        if (t[i+1] < 0) {
            return false;
        }

        t_sum[i+1] = t_sum[i] + t[i+1];
    }
    for (size_t i = 0; i < 7; i += 1) {
        a[i+1] = a[i] + t[i] * j[i];
        v[i+1] = v[i] + t[i] * (a[i] + t[i] * j[i] / 2);
        p[i+1] = p[i] + t[i] * (v[i] + t[i] * (a[i] / 2 + t[i] * j[i] / 6));
    }

    // Velocity and acceleration limits can be broken in the beginning if the initial velocity and acceleration are too high
    // std::cout << std::setprecision(15) << "target: " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << std::endl;
    return std::all_of(v.begin() + 3, v.end(), [vMax](double vm){ return std::abs(vm) < std::abs(vMax) + 1e-9; })
        && std::all_of(a.begin() + 2, a.end(), [aMax](double am){ return std::abs(am) < std::abs(aMax) + 1e-9; })
        && std::abs(p[7] - pf) < 1e-8 && std::abs(v[7] - vf) < 1e-8 && std::abs(a[7] - af) < 1e-8;
}

bool Profile::check(double tf, double pf, double vf, double af, double vMax, double aMax, const std::array<double, 7>& j) {
    // std::cout << std::setprecision(15) << "target: " << std::abs(t_sum[6]-tf) << " " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << std::endl;
    return check(pf, vf, af, vMax, aMax, j) && std::abs(t_sum[6] - tf) < 1e-8;
}

bool Profile::check(double tf, double pf, double vf, double af, double vMax, double aMax, double jMax, const std::array<double, 7>& j) {
    // std::cout << std::setprecision(15) << "target: " << std::abs(t_sum[6]-tf) << " " << std::abs(p[7]-pf) << " " << std::abs(v[7] - vf) << " " << std::abs(a[7] - af) << std::endl;
    // return check(pf, vf, af, vMax, aMax, j) && std::abs(t_sum[6] - tf) < 1e-8 && std::abs(j) < std::abs(jMax) + 1e-12;
}

std::tuple<double, double, double> Profile::integrate(double t, double p0, double v0, double a0, double j) {
    const double p_new = p0 + t * (v0 + t * (a0 / 2 + t * j / 6));
    const double v_new = v0 + t * (a0 + t * j / 2);
    const double a_new = a0 + t * j;
    return {p_new, v_new, a_new};
}

} // namespace ruckig
