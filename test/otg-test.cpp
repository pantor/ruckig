#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <random>
#include "randomizer.hpp"

#include <ruckig/ruckig.hpp>
#include <ruckig/alternative/quintic.hpp>

#ifdef WITH_REFLEXXES
#include <ruckig/alternative/reflexxes.hpp>
#endif


using namespace ruckig;


namespace ruckig {
    template<size_t DOFs>
    std::ostream& operator<< (std::ostream& os, const InputParameter<DOFs>& value) {
        os << value.to_string();
        return os;
    }
}


template<size_t DOFs, class OTGType>
inline void check_duration(OTGType& otg, InputParameter<DOFs>& input, double duration) {
    OutputParameter<DOFs> output;

    while (otg.update(input, output) == Result::Working) {
        input.current_position = output.new_position;
        input.current_velocity = output.new_velocity;
        input.current_acceleration = output.new_acceleration;
    }

    CHECK( output.trajectory.duration == doctest::Approx(duration) );
}


template<size_t DOFs, class OTGType>
inline void check_calculation(OTGType& otg, InputParameter<DOFs>& input) {
    OutputParameter<DOFs> output;

    CAPTURE( input );

    auto result = otg.update(input, output);
    if (result == Result::ErrorTrajectoryDuration) {
        return;
    }

    CHECK( (result == Result::Working || (result == Result::Finished && output.trajectory.duration < 0.005)) );
    CHECK( output.trajectory.duration >= 0.0 );

    for (size_t dof = 0; dof < DOFs; ++dof) {
        CHECK_FALSE( (std::isnan(output.new_position[dof]) || std::isnan(output.new_velocity[dof]) || std::isnan(output.new_acceleration[dof])) );
    }
}


template<size_t DOFs, class OTGType, class OTGCompType>
inline void check_comparison(OTGType& otg, InputParameter<DOFs>& input, OTGCompType& otg_comparison) {
    OutputParameter<DOFs> output;

    CAPTURE( input );

    auto result = otg.update(input, output);
    if (result == Result::ErrorTrajectoryDuration) {
        return;
    }

    CHECK( (result == Result::Working || (result == Result::Finished && output.trajectory.duration < 0.005)) );

    OutputParameter<DOFs> output_comparison;
    auto result_comparison = otg_comparison.update(input, output_comparison);
    CHECK( output.trajectory.duration <= doctest::Approx(output_comparison.trajectory.duration) );
}


int seed {42};
size_t number_trajectories {120000}; // Some user variable you want to be able to set
size_t random_1, random_3, random_discrete_3, random_direction_3, comparison_1, comparison_3, velocity_random_3;

std::normal_distribution<double> position_dist {0.0, 4.0};
std::normal_distribution<double> dynamic_dist {0.0, 0.8};
std::uniform_real_distribution<double> limit_dist {0.08, 12.0};
std::uniform_real_distribution<double> min_limit_dist {-12.0, -0.08};


TEST_CASE("quintic") {
    InputParameter<3> input;
    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {1.0, 1.0, 1.0};
    input.target_velocity = {0.0, 0.0, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {1.0, 1.0, 1.0};
    input.max_acceleration = {1.0, 1.0, 1.0};
    input.max_jerk = {1.0, 1.0, 1.0};

    Quintic<3> otg {0.005};
    check_duration(otg, input, 3.9148676412);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.max_jerk = {2.0, 2.0, 2.0};
    check_duration(otg, input, 3.107232506);
}

TEST_CASE("known" * doctest::description("Known examples")) {
    Ruckig<3> otg {0.005};

    InputParameter<3> input;
    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {0.0, 0.0, 0.0};
    input.target_velocity = {0.0, 0.0, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {1.0, 1.0, 1.0};
    input.max_acceleration = {1.0, 1.0, 1.0};
    input.max_jerk = {1.0, 1.0, 1.0};
    check_duration(otg, input, 0.0);

    input.target_position = {1.0, 1.0, 1.0};
    check_duration(otg, input, 3.1748021039);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.max_jerk = {2.0, 2.0, 2.0};
    check_duration(otg, input, 2.5615528128);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {0.6, 0.6, 0.6};
    check_duration(otg, input, 2.7666666667);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {0.4, 0.4, 0.4};
    check_duration(otg, input, 3.394427191);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.3, 0.3, 0.3};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {1.0, 1.0, 1.0};
    check_duration(otg, input, 2.2319602829);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.3, 0.3, 0.3};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {0.6, 0.6, 0.6};
    check_duration(otg, input, 2.410315834);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {-1.0, -1.0, -1.0};
    check_duration(otg, input, 2.7666666667);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.2, 0.2, 0.2};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {-1.0, -1.0, -1.0};
    input.max_velocity = {10.0, 10.0, 10.0};
    input.max_acceleration = {10.0, 10.0, 10.0};
    check_duration(otg, input, 2.7338531701);

    input.current_position = {-1.0, -1.0, -1.0};
    input.current_velocity = {0.2, 0.2, 0.2};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {1.0, 1.0, 1.0};
    input.max_velocity = {0.4, 0.4, 0.4};
    input.max_acceleration = {1.0, 1.0, 1.0};
    check_duration(otg, input, 5.6053274785);
}

TEST_CASE("random_discrete_3" * doctest::description("Random discrete input with 3 DoF and target velocity, acceleration")) {
    constexpr size_t DOFs {3};
    Ruckig<DOFs, true> otg {0.005};
    InputParameter<DOFs> input;

    std::uniform_int_distribution<int> position_discrete_dist(-1, 1);
    std::uniform_int_distribution<int> dynamic_discrete_dist(-1, 1);
    std::uniform_int_distribution<int> limit_discrete_dist(1, 2);

    Randomizer<DOFs, decltype(position_discrete_dist)> p { position_discrete_dist, seed };
    Randomizer<DOFs, decltype(dynamic_discrete_dist)> d { dynamic_discrete_dist, seed + 1 };
    Randomizer<DOFs, decltype(limit_discrete_dist)> l { limit_discrete_dist, seed + 2 };

    for (size_t i = 0; i < random_discrete_3; ++i) {
        p.fill(input.current_position);
        d.fill(input.current_velocity);
        d.fill(input.current_acceleration);
        p.fill(input.target_position);
        d.fill(input.target_velocity);
        d.fill(input.target_acceleration);
        l.fill(input.max_velocity, input.target_velocity);
        l.fill(input.max_acceleration, input.target_acceleration);
        l.fill(input.max_jerk);

        if (!otg.validate_input(input)) {
            --i;
            continue;
        }

        check_calculation(otg, input);
    }
}

TEST_CASE("random_1" * doctest::description("Random input with 1 DoF and target velocity, acceleration")) {
    constexpr size_t DOFs {1};
    Ruckig<DOFs, true> otg {0.005};
    InputParameter<DOFs> input;

    Randomizer<DOFs, decltype(position_dist)> p { position_dist, seed };
    Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist, seed + 1 };
    Randomizer<DOFs, decltype(limit_dist)> l { limit_dist, seed + 2 };

    for (size_t i = 0; i < random_1; ++i) {
        p.fill(input.current_position);
        d.fill_or_zero(input.current_velocity, 0.9);
        d.fill_or_zero(input.current_acceleration, 0.8);
        p.fill(input.target_position);
        d.fill_or_zero(input.target_velocity, 0.7);
        d.fill_or_zero(input.target_acceleration, 0.6);
        l.fill(input.max_velocity, input.target_velocity);
        l.fill(input.max_acceleration, input.target_acceleration);
        l.fill(input.max_jerk);

        if (!otg.validate_input(input)) {
            --i;
            continue;
        }

        check_calculation(otg, input);
    }
}

TEST_CASE("velocity_random_3" * doctest::description("Random input with 3 DoF and target velocity, acceleration in velocity interface")) {
    constexpr size_t DOFs {3};
    Ruckig<DOFs, true> otg {0.005};
    InputParameter<DOFs> input;
    input.interface = Interface::Velocity;

    Randomizer<DOFs, decltype(position_dist)> p { position_dist, seed + 3 };
    Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist, seed + 4 };
    Randomizer<DOFs, decltype(limit_dist)> l { limit_dist, seed + 5 };

    for (size_t i = 0; i < velocity_random_3; ++i) {
        d.fill_or_zero(input.current_velocity, 0.9);
        d.fill_or_zero(input.current_acceleration, 0.8);
        d.fill_or_zero(input.target_velocity, 0.7);
        d.fill_or_zero(input.target_acceleration, 0.6);
        l.fill(input.max_acceleration, input.target_acceleration);
        l.fill(input.max_jerk);

        check_calculation(otg, input);
    }
}

TEST_CASE("random_3" * doctest::description("Random input with 3 DoF and target velocity, acceleration")) {
    constexpr size_t DOFs {3};
    Ruckig<DOFs, true> otg {0.005};
    InputParameter<DOFs> input;

    Randomizer<DOFs, decltype(position_dist)> p { position_dist, seed + 3 };
    Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist, seed + 4 };
    Randomizer<DOFs, decltype(limit_dist)> l { limit_dist, seed + 5 };

    for (size_t i = 0; i < random_3; ++i) {
        p.fill(input.current_position);
        d.fill_or_zero(input.current_velocity, 0.9);
        d.fill_or_zero(input.current_acceleration, 0.8);
        p.fill(input.target_position);
        d.fill_or_zero(input.target_velocity, 0.7);
        d.fill_or_zero(input.target_acceleration, 0.6);
        l.fill(input.max_velocity, input.target_velocity);
        l.fill(input.max_acceleration, input.target_acceleration);
        l.fill(input.max_jerk);

        if (!otg.validate_input(input)) {
            --i;
            continue;
        }

        check_calculation(otg, input);
    }
}

TEST_CASE("random_direction_3" * doctest::description("Random input with 3 DoF and target velocity, acceleration and min velocity, acceleration")) {
    constexpr size_t DOFs {3};
    Ruckig<DOFs, true> otg {0.005};
    InputParameter<DOFs> input;

    Randomizer<DOFs, decltype(position_dist)> p { position_dist, seed + 3 };
    Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist, seed + 4 };
    Randomizer<DOFs, decltype(limit_dist)> l { limit_dist, seed + 5 };
    Randomizer<DOFs, decltype(min_limit_dist)> min_l { min_limit_dist, seed + 5 };

    for (size_t i = 0; i < random_direction_3; ++i) {
        p.fill(input.current_position);
        d.fill_or_zero(input.current_velocity, 0.9);
        d.fill_or_zero(input.current_acceleration, 0.8);
        p.fill(input.target_position);
        d.fill_or_zero(input.target_velocity, 0.7);
        d.fill_or_zero(input.target_acceleration, 0.6);
        l.fill(input.max_velocity, input.target_velocity);
        l.fill(input.max_acceleration, input.target_acceleration);
        l.fill(input.max_jerk);
        min_l.fill_min(*input.min_velocity, input.target_velocity);
        min_l.fill_min(*input.min_acceleration, input.target_acceleration);

        if (!otg.validate_input(input)) {
            --i;
            continue;
        }

        check_calculation(otg, input);
    }
}

#ifdef WITH_REFLEXXES
TEST_CASE("comparison_1" * doctest::description("Comparison with Reflexxes with 1 DoF")) {
    constexpr size_t DOFs {1};
    Ruckig<DOFs, true> otg {0.005};
    Reflexxes<DOFs> rflx {0.005};
    InputParameter<DOFs> input;

    Randomizer<DOFs, decltype(position_dist)> p { position_dist, seed + 6 };
    Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist, seed + 7 };
    Randomizer<DOFs, decltype(limit_dist)> l { limit_dist, seed + 8 };

    for (size_t i = 0; i < comparison_1; ++i) {
        p.fill(input.current_position);
        d.fill_or_zero(input.current_velocity, 0.9);
        d.fill_or_zero(input.current_acceleration, 0.8);
        p.fill(input.target_position);
        d.fill_or_zero(input.target_velocity, 0.7);
        l.fill(input.max_velocity, input.target_velocity);
        l.fill(input.max_acceleration);
        l.fill(input.max_jerk);

        if (!otg.validate_input(input)) {
            --i;
            continue;
        }

        check_comparison(otg, input, rflx);
    }
    // WARN(counter_faster << " / " << 128*1024 << " trajectories are faster.");
}

TEST_CASE("comparison_3" * doctest::description("Comparison with Reflexxes with 3 DoF")) {
    constexpr size_t DOFs {3};
    Ruckig<DOFs, true> otg {0.005};
    Reflexxes<DOFs> rflx {0.005};
    InputParameter<DOFs> input;

    Randomizer<DOFs, decltype(position_dist)> p { position_dist, seed + 9 };
    Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist, seed + 10 };
    Randomizer<DOFs, decltype(limit_dist)> l { limit_dist, seed + 11 };

    for (size_t i = 0; i < comparison_3; ++i) {
        p.fill(input.current_position);
        d.fill_or_zero(input.current_velocity, 0.9);
        d.fill_or_zero(input.current_acceleration, 0.8);
        p.fill(input.target_position);
        d.fill_or_zero(input.target_velocity, 0.7);
        l.fill(input.max_velocity);
        l.fill(input.max_acceleration);
        l.fill(input.max_jerk);

        if (!otg.validate_input(input)) {
            --i;
            continue;
        }

        check_comparison(otg, input, rflx);
    }
}
#endif


int main(int argc, char** argv) {
    doctest::Context context;

    if (argc > 1 && std::isdigit(argv[1][0])) {
        number_trajectories = std::stoi(argv[1]);
    }
    if (argc > 2 && std::isdigit(argv[2][0])) {
        seed = std::stoi(argv[2]);
    }

    context.applyCommandLine(argc, argv);

    comparison_1 = std::min<size_t>(250000, number_trajectories / 10);
    comparison_3 = std::min<size_t>(250000, number_trajectories / 10);
    random_discrete_3 = std::min<size_t>(250000, number_trajectories / 10);
    random_1 = number_trajectories / 10;
    random_direction_3 = number_trajectories / 50;
    velocity_random_3 = number_trajectories / 10;
    random_3 = number_trajectories - (random_1 + random_direction_3 + comparison_1 + comparison_3 + velocity_random_3 + random_discrete_3);
    std::cout << "<number_trajectories> Random 1 DoF: " << random_1 << " 3 DoF: " << random_3 << "  Comparison 1 DoF: " << comparison_1 << " 3 DoF: " << comparison_3 << " Total: " << number_trajectories << std::endl;

    return context.run();
}
