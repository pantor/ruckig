#pragma once

#include <array>
#include <random>
#include <type_traits>


template<size_t DOFs, class T>
class Randomizer {
    template<class U> using Vector = typename std::conditional<DOFs >= 1, std::array<U, DOFs>, std::vector<U>>::type;

    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform_dist;
    T dist;

public:
    explicit Randomizer() { }
    explicit Randomizer(T dist, int seed): dist(dist), uniform_dist(std::uniform_real_distribution<double>(0.0, 1.0)) {
        gen.seed(seed);
    }

    void fill(Vector<double>& input) {
        for (size_t dof = 0; dof < input.size(); ++dof) {
            input[dof] = dist(gen);
        }
    }

    void fill_or_zero(Vector<double>& input, double p) {
        for (size_t dof = 0; dof < input.size(); ++dof) {
            input[dof] = uniform_dist(gen) < p ? dist(gen) : 0.0;
        }
    }

    void fill(Vector<double>& input, const Vector<double>& offset) {
        for (size_t dof = 0; dof < input.size(); ++dof) {
            input[dof] = dist(gen) + std::abs(offset[dof]);
        }
    }

    void fill_min(Vector<double>& input, const Vector<double>& offset) {
        for (size_t dof = 0; dof < input.size(); ++dof) {
            input[dof] = dist(gen) - std::abs(offset[dof]);
        }
    }
};
