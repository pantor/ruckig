#define CATCH_CONFIG_MAIN
#include <random>

#include <catch2/catch.hpp>
#include <Eigen/Core>

#include <ruckig/parameter.hpp>
#include <ruckig/ruckig.hpp>

#include <ruckig/alternative/quintic.hpp>

#ifdef WITH_REFLEXXES
#include <ruckig/alternative/reflexxes.hpp>
#endif


using namespace ruckig;
using Vec1 = InputParameter<1>::Vector;
using Vec2 = InputParameter<2>::Vector;
using Vec = InputParameter<3>::Vector;


template<size_t DOFs, class OTGType>
void check(OTGType& otg, InputParameter<DOFs>& input, double time) {
    OutputParameter<DOFs> output;

    while (otg.update(input, output) == Result::Working) {
        input.current_position = output.new_position;
        input.current_velocity = output.new_velocity;
        input.current_acceleration = output.new_acceleration;
    }

    CHECK( output.duration == Approx(time).margin(0.005) );
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

    if (output.duration == Approx(output_comparison.duration)) {
        double half_duration = output.duration / 2;
        otg.at_time(half_duration, output);
        otg_comparison.at_time(half_duration, output_comparison);

        for (size_t dof = 0; dof < DOFs; dof += 1) {
            CHECK( output.new_position[dof] == Approx(output_comparison.new_position[dof]).margin(1e-8) );
        }
    } else {
        WARN("Ruckig and Reflexxes differ! Maybe Reflexxes error...");
    }
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
    check(otg, input, 3.915);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.max_jerk = {2.0, 2.0, 2.0};
    check(otg, input, 3.110);
}

TEST_CASE("Ruckig") {
    bool full {true};

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
        check(otg, input, 0.0);

        input.target_position = {1.0, 1.0, 1.0};
        check(otg, input, 3.170);

        input.current_position = {0.0, 0.0, 0.0};
        input.current_velocity = {0.0, 0.0, 0.0};
        input.current_acceleration = {0.0, 0.0, 0.0};
        input.max_jerk = {2.0, 2.0, 2.0};
        check(otg, input, 2.560);

        input.current_position = {0.0, 0.0, 0.0};
        input.current_velocity = {0.0, 0.0, 0.0};
        input.current_acceleration = {0.0, 0.0, 0.0};
        input.max_velocity = {0.6, 0.6, 0.6};
        check(otg, input, 2.765);

        input.current_position = {0.0, 0.0, 0.0};
        input.current_velocity = {0.0, 0.0, 0.0};
        input.current_acceleration = {0.0, 0.0, 0.0};
        input.max_velocity = {0.4, 0.4, 0.4};
        check(otg, input, 3.390);

        input.current_position = {0.0, 0.0, 0.0};
        input.current_velocity = {0.3, 0.3, 0.3};
        input.current_acceleration = {0.0, 0.0, 0.0};
        input.max_velocity = {1.0, 1.0, 1.0};
        check(otg, input, 2.230);

        input.current_position = {0.0, 0.0, 0.0};
        input.current_velocity = {0.3, 0.3, 0.3};
        input.current_acceleration = {0.0, 0.0, 0.0};
        input.max_velocity = {0.6, 0.6, 0.6};
        check(otg, input, 2.410);

        input.current_position = {0.0, 0.0, 0.0};
        input.current_velocity = {0.0, 0.0, 0.0};
        input.current_acceleration = {0.0, 0.0, 0.0};
        input.target_position = {-1.0, -1.0, -1.0};
        check(otg, input, 2.765);

        input.current_position = {0.0, 0.0, 0.0};
        input.current_velocity = {0.2, 0.2, 0.2};
        input.current_acceleration = {0.0, 0.0, 0.0};
        input.target_position = {-1.0, -1.0, -1.0};
        input.max_velocity = {10.0, 10.0, 10.0};
        input.max_acceleration = {10.0, 10.0, 10.0};
        check(otg, input, 2.730);

        input.current_position = {-1.0, -1.0, -1.0};
        input.current_velocity = {0.2, 0.2, 0.2};
        input.current_acceleration = {0.0, 0.0, 0.0};
        input.target_position = {1.0, 1.0, 1.0};
        input.max_velocity = {0.4, 0.4, 0.4};
        input.max_acceleration = {1.0, 1.0, 1.0};
        check(otg, input, 5.605);
    }

    SECTION("Random input with 3 DoF") {
        Ruckig<3> otg {0.005};
        InputParameter<3> input;

        // Eigen returns uniform random floats between -1 and 1
        srand(42);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < (full ? 64*1024 : 1*1024); i += 1) {
            input.current_position = Vec::Random();
            input.current_velocity = dist(gen) < 0.9 ? (Vec)Vec::Random() : (Vec)Vec::Zero();
            input.current_acceleration = dist(gen) < 0.8 ? (Vec)Vec::Random() : (Vec)Vec::Zero();
            input.target_position = Vec::Random();
            input.max_velocity = 10 * Vec::Random().array().abs() + input.target_velocity.array().abs() + 0.1;
            input.max_acceleration = 10 * Vec::Random().array().abs() + 0.1;
            input.max_jerk = 10 * Vec::Random().array().abs() + 0.1;

            check_calculation(otg, input);
        }
    }

    SECTION("Random input with 3 DoF, target velocity") {
        Ruckig<3> otg {0.005};
        InputParameter<3> input;

        // Eigen returns uniform random floats between -1 and 1
        srand(39);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < (full ? 32*1024 : 1*1024); i += 1) {
            input.current_position = Vec::Random();
            input.current_velocity = dist(gen) < 0.9 ? (Vec)Vec::Random() : (Vec)Vec::Zero();
            input.current_acceleration = dist(gen) < 0.8 ? (Vec)Vec::Random() : (Vec)Vec::Zero();
            input.target_position = Vec::Random();
            input.target_velocity = dist(gen) < 0.6 ? (Vec)Vec::Random() : (Vec)Vec::Zero();
            input.max_velocity = 10 * Vec::Random().array().abs() + input.target_velocity.array().abs() + 0.1;
            input.max_acceleration = 10 * Vec::Random().array().abs() + 0.1;
            input.max_jerk = 10 * Vec::Random().array().abs() + 0.1;

            check_calculation(otg, input);
        }
    }

    SECTION("Random input with 1 DoF with target velocity, acceleration") {
        Ruckig<1> otg {0.005};
        InputParameter<1> input;

        srand(47);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < (full ? 32*1024 : 1*1024); i += 1) {
            input.current_position = Vec1::Random();
            input.current_velocity = dist(gen) < 0.9 ? (Vec1)Vec1::Random() : (Vec1)Vec1::Zero();
            input.current_acceleration = dist(gen) < 0.8 ? (Vec1)Vec1::Random() : (Vec1)Vec1::Zero();
            input.target_position = Vec1::Random();
            input.target_velocity = dist(gen) < 0.6 ? (Vec1)Vec1::Random() : (Vec1)Vec1::Zero();
            input.target_acceleration = dist(gen) < 0.4 ? (Vec1)Vec1::Random() : (Vec1)Vec1::Zero();
            input.max_velocity = 10 * Vec1::Random().array().abs() + input.target_velocity.array().abs(); // Target velocity needs to be smaller than max velocity
            input.max_acceleration = 10 * Vec1::Random().array().abs() + input.target_acceleration.array().abs() + 0.1;
            input.max_jerk = 10 * Vec1::Random().array().abs() + 0.1;

            auto max_target_acceleration = (2 * input.max_jerk.array() * (input.max_velocity.array() - input.target_velocity.array().abs())).sqrt();
            if ((input.target_acceleration.array().abs() > max_target_acceleration.array()).any()) {
                continue;
            }

            check_calculation(otg, input);
        }
    }

#ifdef WITH_REFLEXXES
    SECTION("Comparison with Reflexxes with 1 DoF") {
        Ruckig<1> otg {0.005};
        Reflexxes<1> rflx {0.005};
        InputParameter<1> input;

        srand(43);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < (full ? 64*1024 : 1024); i += 1) {
            input.current_position = Vec1::Random();
            input.current_velocity = dist(gen) < 0.9 ? (Vec1)Vec1::Random() : (Vec1)Vec1::Zero();
            input.current_acceleration = dist(gen) < 0.8 ? (Vec1)Vec1::Random() : (Vec1)Vec1::Zero();
            input.target_position = Vec1::Random();
            input.target_velocity = dist(gen) < 0.6 ? (Vec1)Vec1::Random() : (Vec1)Vec1::Zero();
            input.max_velocity = 10 * Vec1::Random().array().abs() + input.target_velocity.array().abs(); // Target velocity needs to be smaller than max velocity
            input.max_acceleration = 10 * Vec1::Random().array().abs() + 0.1;
            input.max_jerk = 10 * Vec1::Random().array().abs() + 0.1;

            check_comparison(otg, input, rflx);
        }
    }

    SECTION("Comparison with Reflexxes with 2 DoF") {
        Ruckig<2> otg {0.005};
        Reflexxes<2> rflx {0.005};
        InputParameter<2> input;

        srand(48);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < 1*1024; i += 1) {
            input.current_position = Vec2::Random();
            input.current_velocity = dist(gen) < 0.9 ? (Vec2)Vec2::Random() : (Vec2)Vec2::Zero();
            input.current_acceleration = dist(gen) < 0.8 ? (Vec2)Vec2::Random() : (Vec2)Vec2::Zero();
            input.target_position = Vec2::Random();
            input.max_velocity = 10 * Vec2::Random().array().abs() + 0.1;
            input.max_acceleration = 10 * Vec2::Random().array().abs() + 0.1;
            input.max_jerk = 10 * Vec2::Random().array().abs() + 0.1;

            check_comparison(otg, input, rflx);
        }
    }

    /* SECTION("Comparison with Reflexxes with 3 DoF") {
        Ruckig<3> otg {0.005};
        Reflexxes<3> rflx {0.005};
        InputParameter<3> input;

        srand(44);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < 16; i += 1) {
            input.current_position = Vec::Random();
            input.current_velocity = dist(gen) < 0.9 ? (Vec)Vec::Random() : (Vec)Vec::Zero();
            input.current_acceleration = dist(gen) < 0.8 ? (Vec)Vec::Random() : (Vec)Vec::Zero();
            input.target_position = Vec::Random();
            input.max_velocity = 10 * Vec::Random().array().abs() + 0.1;
            input.max_acceleration = 10 * Vec::Random().array().abs() + 0.1;
            input.max_jerk = 10 * Vec::Random().array().abs() + 0.1;

            check_comparison(otg, input, rflx);
        }
    } */
#endif
}
