#pragma once

#include <limits>
#include <optional>

#include <ruckig/profile.hpp>


namespace ruckig {

//! Which times are possible for synchronization?
struct Block {
    struct Interval {
        double left, right; // [s]
        Profile profile; // Profile corresponding to right (end) time

        explicit Interval(double left, double right, const Profile& profile): left(left), right(right), profile(profile) { };
    };

    Profile p_min; // Save min profile so that it doesn't need to be recalculated in Step2
    double t_min; // [s]

    // Max. 2 intervals can be blocked: called a and b with corresponding profiles, order does not matter
    std::optional<Interval> a, b;

    explicit Block() { }
    explicit Block(const Profile& p_min): p_min(p_min), t_min(p_min.t_sum[6] + p_min.t_brake.value_or(0.0)) { }

    bool is_blocked(double t) const {
        return (t < t_min) || (a && a->left < t && t < a->right) || (b && b->left < t && t < b->right);
    }

    inline static void add_interval(std::optional<Block::Interval>& interval, const Profile& left, const Profile& right) {
        const double left_duration = left.t_sum[6] + left.t_brake.value_or(0.0);
        const double right_duraction = right.t_sum[6] + right.t_brake.value_or(0.0);
        if (left_duration < right_duraction) {
            interval = Block::Interval(left_duration, right_duraction, right);
        } else {
            interval = Block::Interval(right_duraction, left_duration, left);
        }
    }

    template<size_t N>
    static bool calculate_block(Block& block, const std::array<Profile, N>& valid_profiles, const size_t valid_profile_counter) {
        // if (valid_profile_counter > 0)
        // {
        //     std::cout << "---\n " << valid_profile_counter << std::endl;
        //     for (size_t i = 0; i < valid_profile_counter; ++i) {
        //         std::cout << valid_profiles[i].t_sum[6] << " " << valid_profiles[i].to_string() << std::endl;
        //     }
        // }

        if (
            (valid_profile_counter == 1)
            || (valid_profile_counter == 2 && std::abs(valid_profiles[0].t_sum[6] - valid_profiles[1].t_sum[6]) < std::numeric_limits<double>::epsilon())
        ) {
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

            Block::add_interval(block.a, valid_profiles[idx_else_1], valid_profiles[idx_else_2]);
            return true;

        } else if (valid_profile_counter == 5) {
            size_t idx_else_1 = (idx_min + 1) % 5;
            size_t idx_else_2 = (idx_min + 2) % 5;
            size_t idx_else_3 = (idx_min + 3) % 5;
            size_t idx_else_4 = (idx_min + 4) % 5;

            if (valid_profiles[idx_else_1].direction == valid_profiles[idx_else_2].direction) {
                Block::add_interval(block.a, valid_profiles[idx_else_1], valid_profiles[idx_else_2]);
                Block::add_interval(block.b, valid_profiles[idx_else_3], valid_profiles[idx_else_4]);
            } else {
                Block::add_interval(block.a, valid_profiles[idx_else_1], valid_profiles[idx_else_4]);
                Block::add_interval(block.b, valid_profiles[idx_else_2], valid_profiles[idx_else_3]);
            }
            return true;
        }

        return false;
    }

    template<size_t DOFs>
    static bool synchronize(const std::array<Block, DOFs>& blocks, std::optional<double> t_min, double& t_sync, int& limiting_dof, std::array<Profile, DOFs>& profiles, bool discrete_duration, double delta_time) {
        if (DOFs == 1 && !t_min && !discrete_duration) {
            limiting_dof = 0;
            t_sync = blocks[0].t_min;
            profiles[0] = blocks[0].p_min;
            return true;
        }

        // Possible t_syncs are the start times of the intervals and optional t_min
        std::array<double, 3*DOFs+1> possible_t_syncs;
        std::array<int, 3*DOFs+1> idx;
        for (size_t dof = 0; dof < DOFs; ++dof) {
            possible_t_syncs[3 * dof] = blocks[dof].t_min;
            possible_t_syncs[3 * dof + 1] = blocks[dof].a ? blocks[dof].a->right : std::numeric_limits<double>::infinity();
            possible_t_syncs[3 * dof + 2] = blocks[dof].b ? blocks[dof].b->right : std::numeric_limits<double>::infinity();
        }
        possible_t_syncs[3 * DOFs] = t_min.value_or(std::numeric_limits<double>::infinity());

        if (discrete_duration) {
            for (size_t i = 0; i < 3*DOFs+1; ++i) {
                possible_t_syncs[i] = std::ceil(possible_t_syncs[i] / delta_time) * delta_time;
            }
        }

        // Test them in sorted order
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&possible_t_syncs](size_t i, size_t j) { return possible_t_syncs[i] < possible_t_syncs[j]; });

        // Start at last tmin (or worse)
        for (auto i = idx.begin() + DOFs - 1; i != idx.end(); ++i) {
            const double possible_t_sync = possible_t_syncs[*i];
            if (std::any_of(blocks.begin(), blocks.end(), [possible_t_sync](auto block){ return block.is_blocked(possible_t_sync); }) || possible_t_sync < t_min.value_or(0.0)) {
                continue;
            }

            t_sync = possible_t_sync;
            if (*i == 3*DOFs) { // Optional t_min
                limiting_dof = -1;
                return true;
            }

            const auto div = std::div(*i, 3);
            limiting_dof = div.quot;
            switch (div.rem) {
                case 0: {
                    profiles[limiting_dof] = blocks[limiting_dof].p_min;
                } break;
                case 1: {
                    profiles[limiting_dof] = blocks[limiting_dof].a->profile;
                } break;
                case 2: {
                    profiles[limiting_dof] = blocks[limiting_dof].b->profile;
                } break;
            }
            return true;
        }

        return false;
    }
};

} // namespace ruckig
