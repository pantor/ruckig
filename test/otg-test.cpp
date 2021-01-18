#define CATCH_CONFIG_MAIN

#include <random>

#include <catch2/catch.hpp>
#include "randomizer.hpp"

#include <ruckig/parameter.hpp>
#include <ruckig/ruckig.hpp>
#include <ruckig/alternative/quintic.hpp>

#ifdef WITH_REFLEXXES
#include <ruckig/alternative/reflexxes.hpp>
#endif


using namespace ruckig;

template<size_t DOFs, class OTGType>
void check_duration(OTGType& otg, InputParameter<DOFs>& input, double duration) {
    OutputParameter<DOFs> output;

    while (otg.update(input, output) == Result::Working) {
        input.current_position = output.new_position;
        input.current_velocity = output.new_velocity;
        input.current_acceleration = output.new_acceleration;
    }

    CHECK( output.duration == Approx(duration) );
}


template<size_t DOFs, class OTGType>
void check_calculation(OTGType& otg, InputParameter<DOFs>& input) {
    OutputParameter<DOFs> output;

    CAPTURE( input.current_position, input.current_velocity, input.current_acceleration );
    CAPTURE( input.target_position, input.target_velocity, input.target_acceleration );
    CAPTURE( input.max_velocity, input.max_acceleration, input.max_jerk );

    auto result = otg.update(input, output);

    CHECK( result == Result::Working );
    CHECK( output.duration > 0.0 );

    for (size_t dof = 0; dof < DOFs; dof += 1) {
        CHECK_FALSE( (std::isnan(output.new_position[dof]) || std::isnan(output.new_velocity[dof]) || std::isnan(output.new_acceleration[dof])) );
    }
}


template<size_t DOFs, class OTGType, class OTGCompType>
void check_comparison(OTGType& otg, InputParameter<DOFs>& input, OTGCompType& otg_comparison) {
    OutputParameter<DOFs> output;

    CAPTURE( input.current_position, input.current_velocity, input.current_acceleration );
    CAPTURE( input.target_position, input.target_velocity, input.target_acceleration );
    CAPTURE( input.max_velocity, input.max_acceleration, input.max_jerk );

    auto result = otg.update(input, output);
    CHECK( result == Result::Working );

    OutputParameter<DOFs> output_comparison;
    auto result_comparison = otg_comparison.update(input, output_comparison);
    CHECK( output.duration <= Approx(output_comparison.duration) );

    // if (output.duration == Approx(output_comparison.duration)) {
    //     double half_duration = output.duration / 2;
    //     otg.at_time(half_duration, output);
    //     otg_comparison.at_time(half_duration, output_comparison);

    //     for (size_t dof = 0; dof < DOFs; dof += 1) {
    //         CHECK( output.new_position[dof] == Approx(output_comparison.new_position[dof]).margin(1e-8) );
    //     }
    // } else {
    //     WARN("Ruckig and Reflexxes differ! Maybe Reflexxes error...");
    // }
}


TEST_CASE("Quintic") {
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

TEST_CASE("Ruckig") {
    constexpr bool full {true};

    std::normal_distribution<double> position_dist {0.0, 4.0};
    std::normal_distribution<double> dynamic_dist {0.0, 0.8};
    std::uniform_real_distribution<double> limit_dist {0.04, 12.0};

    SECTION("Known examples") {
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

    SECTION("Random input with 3 DoF and target velocity") {
        constexpr size_t DOFs {3};
        Ruckig<DOFs, true> otg {0.005};
        InputParameter<DOFs> input;

        srand(39);
        Randomizer<DOFs, decltype(position_dist)> p { position_dist };
        Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist };
        Randomizer<DOFs, decltype(limit_dist)> l { limit_dist };

        for (size_t i = 0; i < (full ? 256 : 1) * 1024; ++i) {
            p.fill(input.current_position);
            d.fill_or_zero(input.current_velocity, 0.9);
            d.fill_or_zero(input.current_acceleration, 0.8);
            p.fill(input.target_position);
            d.fill_or_zero(input.target_velocity, 0.7);
            l.fill(input.max_velocity, input.target_velocity);
            l.fill(input.max_acceleration);
            l.fill(input.max_jerk);
            check_calculation(otg, input);
        }
    }

    SECTION("Random input with 1 DoF and target velocity, acceleration") {
        constexpr size_t DOFs {1};
        Ruckig<DOFs, true> otg {0.005};
        InputParameter<DOFs> input;

        srand(47);
        Randomizer<DOFs, decltype(position_dist)> p { position_dist };
        Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist };
        Randomizer<DOFs, decltype(limit_dist)> l { limit_dist };

        for (size_t i = 0; i < (full ? 256 : 1) * 1024; ++i) {
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
                continue;
            }

            check_calculation(otg, input);
        }
    }

    SECTION("Random input with 3 DoF and target velocity, acceleration") {
        constexpr size_t DOFs {3};
        Ruckig<DOFs, true> otg {0.005};
        InputParameter<DOFs> input;
        
        srand(44);
        Randomizer<DOFs, decltype(position_dist)> p { position_dist };
        Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist };
        Randomizer<DOFs, decltype(limit_dist)> l { limit_dist };

        for (size_t i = 0; i < (full ? 320 : 1) * 1024; ++i) {
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
                continue;
            }

            check_calculation(otg, input);
        }
    }

#ifdef WITH_REFLEXXES
    SECTION("Comparison with Reflexxes with 1 DoF") {
        constexpr size_t DOFs {1};
        Ruckig<DOFs, true> otg {0.005};
        Reflexxes<DOFs> rflx {0.005};
        InputParameter<DOFs> input;

        srand(43);
        Randomizer<DOFs, decltype(position_dist)> p { position_dist };
        Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist };
        Randomizer<DOFs, decltype(limit_dist)> l { limit_dist };

        for (size_t i = 0; i < (full ? 154 : 1) * 1024; ++i) {
            p.fill(input.current_position);
            d.fill_or_zero(input.current_velocity, 0.9);
            d.fill_or_zero(input.current_acceleration, 0.8);
            p.fill(input.target_position);
            d.fill_or_zero(input.target_velocity, 0.6);
            l.fill(input.max_velocity, input.target_velocity);
            l.fill(input.max_acceleration);
            l.fill(input.max_jerk);
            check_comparison(otg, input, rflx);
        }
    }

    /* SECTION("Comparison with Reflexxes with 2 DoF") {
        constexpr size_t DOFs {2};
        using Vec = Eigen::Matrix<double, DOFs, 1>;

        Ruckig<DOFs, true> otg {0.005};
        Reflexxes<DOFs> rflx {0.005};
        InputParameter<DOFs> input;

        srand(48);
        Randomizer<DOFs, decltype(position_dist)> p { position_dist };
        Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist };
        Randomizer<DOFs, decltype(limit_dist)> l { limit_dist };

        for (size_t i = 0; i < 1*1024; ++i) {
            p.fill(input.current_position);
            d.fill_or_zero(input.current_velocity, 0.9);
            d.fill_or_zero(input.current_acceleration, 0.8);
            p.fill(input.target_position);
            l.fill(input.max_velocity);
            l.fill(input.max_acceleration);
            l.fill(input.max_jerk);
            check_comparison(otg, input, rflx);
        }
    } */
#endif
}
