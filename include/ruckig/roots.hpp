#pragma once

#define _USE_MATH_DEFINES
#include <cfloat>
#include <cmath>
#include <set>


namespace ruckig {

inline double Power(double v, int e) {
    return std::pow(v, e);
}

inline double Sqrt(double v) {
    return std::sqrt(v);
}

inline double Abs(double v) {
    return std::abs(v);
}

} // namespace ruckig


namespace Roots {

//! Calculate all roots of a*x^3 + b*x^2 + c*x + d = 0
inline std::set<double> solveCub(double a, double b, double c, double d) {
    std::set<double> roots;

    constexpr double cos120 = -0.50;
    constexpr double sin120 = 0.866025403784438646764;

    if (std::abs(d) < DBL_EPSILON) {
        // First solution is x = 0
        roots.insert(0.0);

        // Converting to a quadratic equation
        d = c;
        c = b;
        b = a;
        a = 0.0;
    }

    if (std::abs(a) < DBL_EPSILON) {
        if (std::abs(b) < DBL_EPSILON) {
            // Linear equation
            if (std::abs(c) > DBL_EPSILON) {
                roots.insert(-d / c);
            }
        
        } else {
            // Quadratic equation
            const double discriminant = c * c - 4 * b * d;
            if (discriminant >= 0) {
                const double inv2b = 1.0 / (2 * b);
                const double y = std::sqrt(discriminant);
                roots.insert((-c + y) * inv2b);
                roots.insert((-c - y) * inv2b);
            }
        }

    } else {
        // Cubic equation
        const double inva = 1.0 / a;
        const double invaa = inva * inva;
        const double bb = b * b;
        const double bover3a = b * (1.0 / 3.0) * inva;
        const double p = (3.0 * a * c - bb) * (1.0 / 3.0) * invaa;
        const double halfq = (2 * bb * b - 9.0 * a * b * c + 27.0 * a * a * d) * (0.5 / 27.0) * invaa * inva;
        const double yy = p * p * p / 27.0 + halfq * halfq;

        if (yy > DBL_EPSILON) {
            // Sqrt is positive: one real solution
            const double y = std::sqrt(yy);
            const double uuu = -halfq + y;
            const double vvv = -halfq - y;
            const double www = std::abs(uuu) > std::abs(vvv) ? uuu : vvv;
            const double w = std::cbrt(www);
            roots.insert(w - p / (3.0 * w) - bover3a);
        } else if (yy < -DBL_EPSILON) {
            // Sqrt is negative: three real solutions
            const double x = -halfq;
            const double y = std::sqrt(-yy);
            double theta;
            double r;
            double ux;
            double uyi;
            // Convert to polar form
            if (std::abs(x) > DBL_EPSILON) {
                theta = (x > 0.0) ? std::atan(y / x) : (std::atan(y / x) + M_PI);
                r = std::sqrt(x * x - yy);
            } else {
                // Vertical line
                theta = M_PI / 2;
                r = y;
            }
            // Calculate cube root
            theta /= 3.0;
            r = std::cbrt(r);
            // Convert to complex coordinate
            ux = std::cos(theta) * r;
            uyi = std::sin(theta) * r;
            // First solution
            roots.insert(ux + ux - bover3a);
            // Second solution, rotate +120 degrees
            roots.insert(2 * (ux * cos120 - uyi * sin120) - bover3a);
            // Third solution, rotate -120 degrees
            roots.insert(2 * (ux * cos120 + uyi * sin120) - bover3a);
        } else {
            // Sqrt is zero: two real solutions
            const double www = -halfq;
            const double w = std::cbrt(www);
            // First solution
            roots.insert(w + w - bover3a);
            // Second solution, rotate +120 degrees
            roots.insert(2 * w * cos120 - bover3a);
        }
    }
    return roots;
}

// Solve resolvent eqaution of corresponding Quartic equation
// The input x must be of length 3
// Number of zeros are returned
inline int solveResolvent(double *x, double a, double b, double c) {
    const double a2 = a * a;
    double q = (a2 - 3.0 * b) / 9.0;
    const double r = (a * (2 * a2 - 9.0 * b) + 27.0 * c) / 54.0;
    const double r2 = r * r;
    const double q3 = q * q * q;
    double A, B;
    if (r2 < q3) {
        double t = r / std::sqrt(q3);
        if (t < -1.0) {
            t = -1.0;
        }
        if (t > 1.0) {
            t = 1.0;
        }
        t = std::acos(t);
        a /= 3.0;
        q = -2 * std::sqrt(q);
        x[0] = q * std::cos(t / 3.0) - a;
        x[1] = q * std::cos((t + M_PI * 2) / 3.0) - a;
        x[2] = q * std::cos((t - M_PI * 2) / 3.0) - a;
        return 3;
    } else {
        A = -std::cbrt(std::abs(r) + std::sqrt(r2 - q3));
        if (r < 0.0) {
            A = -A;
        }
        B = (0.0 == A ? 0.0 : q / A);

        a /= 3.0;
        x[0] = (A + B) - a;
        x[1] = -0.5 * (A + B) - a;
        x[2] = 0.5 * std::sqrt(3.0) * (A - B);
        if (std::abs(x[2]) < DBL_EPSILON) {
            x[2] = x[1];
            return 2;
        }

        return 1;
    }
}

//! Calculate all roots of the monic quartic equation: x^4 + a*x^3 + b*x^2 + c*x + d = 0
inline std::set<double> solveQuartMonic(double a, double b, double c, double d) {
    std::set<double> roots;

    const double a3 = -b;
    const double b3 = a * c - 4 * d;
    const double c3 = -a * a * d - c * c + 4 * b * d;

    // Solve the resolvent: y^3 - b*y^2 + (ac - 4*d)*y - a^2*d - c^2 + 4*b*d = 0
    double x3[3];
    int iZeroes = solveResolvent(x3, a3, b3, c3);

    double q1, q2, p1, p2, D, sqrtD, y;

    y = x3[0];
    // Choosing Y with maximal absolute value.
    if (iZeroes != 1) {
        if (std::abs(x3[1]) > std::abs(y)) {
            y = x3[1];
        }
        if (std::abs(x3[2]) > std::abs(y)) {
            y = x3[2];
        }
    }

    // h1 + h2 = y && h1*h2 = d  <=>  h^2 - y*h + d = 0    (h === q)
    D = y * y - 4 * d;
    if (std::abs(D) < DBL_EPSILON) {
        q1 = q2 = y * 0.5;
        // g1 + g2 = a && g1 + g2 = b - y   <=>   g^2 - a*g + b - y = 0    (p === g)
        D = a * a - 4.0 * (b - y);
        if (std::abs(D) < DBL_EPSILON) {
            p1 = p2 = a * 0.5;
        } else {
            sqrtD = std::sqrt(D);
            p1 = (a + sqrtD) * 0.5;
            p2 = (a - sqrtD) * 0.5;
        }
    } else {
        sqrtD = std::sqrt(D);
        q1 = (y + sqrtD) * 0.5;
        q2 = (y - sqrtD) * 0.5;
        // g1 + g2 = a && g1*h2 + g2*h1 = c   ( && g === p )  Krammer
        p1 = (a * q1 - c) / (q1 - q2);
        p2 = (c - a * q2) / (q1 - q2);
    }

    constexpr double eps {2 * DBL_EPSILON};

    // Solve the quadratic equation: x^2 + p1*x + q1 = 0
    D = p1 * p1 - 4 * q1;
    if (std::abs(D) < eps) {
        roots.insert(-p1 * 0.5);
    } else if (D > 0.0) {
        sqrtD = std::sqrt(D);
        roots.insert((-p1 + sqrtD) * 0.5);
        roots.insert((-p1 - sqrtD) * 0.5);
    }

    // Solve the quadratic equation: x^2 + p2*x + q2 = 0
    D = p2 * p2 - 4 * q2;
    if (std::abs(D) < eps) {
        roots.insert(-p2 * 0.5);
    } else if (D > 0.0) {
        sqrtD = std::sqrt(D);
        roots.insert((-p2 + sqrtD) * 0.5);
        roots.insert((-p2 - sqrtD) * 0.5);
    }

    return roots;
}

//! Calculate the quartic equation: x^4 + b*x^3 + c*x^2 + d*x + e = 0
inline std::set<double> solveQuartMonic(const std::array<double, 5>& polynom) {
    return solveQuartMonic(polynom[1], polynom[2], polynom[3], polynom[4]);
}


//! Evaluate a polynomial of order N at x
template<size_t N>
inline double polyEval(std::array<double, N> p, double x) {
    double retVal = 0.0;

    if constexpr (N > 0) {
        if (std::abs(x) < DBL_EPSILON) {
            retVal = p[N - 1];
        } else if (x == 1.0) {
            for (int i = N - 1; i >= 0; i--) {
                retVal += p[i];
            }
        } else {
            double xn = 1.0;

            for (int i = N - 1; i >= 0; i--) {
                retVal += p[i] * xn;
                xn *= x;
            }
        }
    }

    return retVal;
}

// Calculate the derivative poly coefficients of a given poly
template<size_t N>
inline std::array<double, N-1> polyDeri(const std::array<double, N>& coeffs) {
    std::array<double, N-1> deriv;
    for (size_t i = 0; i < N - 1; ++i) {
        deriv[i] = (N - 1 - i) * coeffs[i];
    }
    return deriv;
}

template<size_t N>
inline std::array<double, N-1> polyMonicDeri(const std::array<double, N>& monic_coeffs) {
    std::array<double, N-1> deriv;
    deriv[0] = 1.0;
    for (size_t i = 1; i < N - 1; ++i) {
        deriv[i] = (N - 1 - i) * monic_coeffs[i] / (N - 1);
    }
    return deriv;
}

// Safe Newton Method
// Requirements: f(l)*f(h) <= 0
template <size_t maxIts, typename F, typename DF>
inline double safeNewton(const F& func, const DF& dfunc, double l, double h, const double tol) {
    const double fl = func(l);
    const double fh = func(h);
    if (fl == 0.0) {
        return l;
    }
    if (fh == 0.0) {
        return h;
    }
    if (fl > 0.0) {
        std::swap(l, h);
    }

    double rts = 0.5 * (l + h);
    double dxold = std::abs(h - l);
    double dx = dxold;
    double f = func(rts);
    double df = dfunc(rts);
    double temp;
    for (size_t j = 0; j < maxIts; j++) {
        if ((((rts - h) * df - f) * ((rts - l) * df - f) > 0.0) || (std::abs(2 * f) > std::abs(dxold * df))) {
            dxold = dx;
            dx = 0.5 * (h - l);
            rts = l + dx;
            if (l == rts) {
                break;
            }
        } else {
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts) {
                break;
            }
        }

        if (std::abs(dx) < tol) {
            break;
        }

        f = func(rts);
        df = dfunc(rts);
        if (f < 0.0) {
            l = rts;
        } else {
            h = rts;
        }
    }

    return rts;
}

// Calculate a single zero of polynom p(x) inside [lbound, ubound]
// Requirements: p(lbound)*p(ubound) < 0, lbound < ubound

constexpr double tolerance {1e-14};

template<size_t N, size_t maxDblIts = 128>
inline double shrinkInterval(const std::array<double, N>& p, double lbound, double ubound) {
    auto deriv = polyDeri(p);
    auto func = [&p](double x) { return polyEval(p, x); };
    auto dfunc = [&deriv](double x) { return polyEval(deriv, x); };
    return safeNewton<maxDblIts>(func, dfunc, lbound, ubound, tolerance);
}

} // namespace Roots
