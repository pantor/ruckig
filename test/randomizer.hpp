#pragma once

#include <array>
#include <random>

template<size_t DOFs, class T>
class Randomizer {
    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform_dist;
    T dist;

public:
    explicit Randomizer(T dist): dist(dist), uniform_dist(std::uniform_real_distribution<double>(0.0, 1.0)) { }

    void fill(std::array<double, DOFs>& input) {
        for (size_t dof = 0; dof < DOFs; ++dof) {
            input[dof] = dist(gen);
        }
    }

    void fill_or_zero(std::array<double, DOFs>& input, double p) {
        for (size_t dof = 0; dof < DOFs; ++dof) {
            input[dof] = uniform_dist(gen) < p ? dist(gen) : 0.0;
        }
    }

    void fill(std::array<double, DOFs>& input, const std::array<double, DOFs>& offset) {
        for (size_t dof = 0; dof < DOFs; ++dof) {
            input[dof] = dist(gen) + std::abs(offset[dof]);
        }
    }
};
