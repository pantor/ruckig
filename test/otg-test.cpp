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

template<size_t DOFs>
inline std::array<double, DOFs> Random() {
    using EigenVector = Eigen::Matrix<double, DOFs, 1>;
    EigenVector a = EigenVector::Random(); // Eigen returns uniform random floats between -1 and 1
    std::array<double, DOFs> result;
    std::copy_n(a.data(), DOFs, result.begin());
    return result;
}

template<size_t DOFs>
constexpr inline std::array<double, DOFs> Zero() {
    return std::array<double, DOFs> {};
}

template<size_t DOFs>
inline std::array<double, DOFs> RandomOrZero(double r, double p_random) {
    return r < p_random ? Random<DOFs>() : Zero<DOFs>();
}


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
    check_duration(otg, input, 3.9148676412);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.max_jerk = {2.0, 2.0, 2.0};
    check_duration(otg, input, 3.107232506);
}

TEST_CASE("Ruckig") {
    constexpr bool full {true};
    constexpr double limit_min {0.1};
    constexpr double limit_max = {10.0};

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
        using Vec = Eigen::Matrix<double, DOFs, 1>;

        Ruckig<DOFs, true> otg {0.005};
        InputParameter<DOFs> input;

        srand(39);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < (full ? 256 : 1) * 1024; ++i) {
            input.current_position = Random<DOFs>();
            input.current_velocity = RandomOrZero<DOFs>(dist(gen), 0.9);
            input.current_acceleration = RandomOrZero<DOFs>(dist(gen), 0.8);
            input.target_position = Random<DOFs>();
            input.target_velocity = RandomOrZero<DOFs>(dist(gen), 0.7);

            Vec max_v = limit_max * Vec::Random().array().abs() + Eigen::Map<Vec>(input.target_velocity.data(), 3, 1).array().abs() + limit_min;
            std::copy_n(max_v.data(), DOFs, input.max_velocity.begin());
            Vec max_a = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_a.data(), DOFs, input.max_acceleration.begin());
            Vec max_j = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_j.data(), DOFs, input.max_jerk.begin());

            check_calculation(otg, input);
        }
    }

    SECTION("Random input with 1 DoF and target velocity, acceleration") {
        constexpr size_t DOFs {1};
        using Vec = Eigen::Matrix<double, DOFs, 1>;

        Ruckig<DOFs, true> otg {0.005};
        InputParameter<DOFs> input;

        srand(47);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < (full ? 256 : 1) * 1024; ++i) {
            input.current_position = Random<DOFs>();
            input.current_velocity = RandomOrZero<DOFs>(dist(gen), 0.9);
            input.current_acceleration = RandomOrZero<DOFs>(dist(gen), 0.8);
            input.target_position = Random<DOFs>();
            input.target_velocity = RandomOrZero<DOFs>(dist(gen), 0.7);
            input.target_acceleration = RandomOrZero<DOFs>(dist(gen), 0.6);

            Vec max_v = limit_max * Vec::Random().array().abs() + Eigen::Map<Vec>(input.target_velocity.data()).array().abs() + limit_min;
            std::copy_n(max_v.data(), DOFs, input.max_velocity.begin());
            Vec max_a = limit_max * Vec::Random().array().abs() + Eigen::Map<Vec>(input.target_acceleration.data()).array().abs() + limit_min;
            std::copy_n(max_a.data(), DOFs, input.max_acceleration.begin());
            Vec max_j = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_j.data(), DOFs, input.max_jerk.begin());

            Vec t_v = Eigen::Map<Vec>(input.target_velocity.data());
            Vec t_a = Eigen::Map<Vec>(input.target_acceleration.data());
            auto max_target_acceleration = (2 * max_j.array() * (max_v.array() - t_v.array().abs())).sqrt();
            if ((t_a.array().abs() > max_target_acceleration.array()).any()) {
                continue;
            }

            check_calculation(otg, input);
        }
    }

    SECTION("Random input with 3 DoF and target velocity, acceleration") {
        constexpr size_t DOFs {3};
        using Vec = Eigen::Matrix<double, DOFs, 1>;

        Ruckig<DOFs, true> otg {0.005};
        InputParameter<DOFs> input;
        
        srand(39);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < (full ? 272 : 1) * 1024; ++i) {
            input.current_position = Random<DOFs>();
            input.current_velocity = RandomOrZero<DOFs>(dist(gen), 0.9);
            input.current_acceleration = RandomOrZero<DOFs>(dist(gen), 0.8);
            input.target_position = Random<DOFs>();
            input.target_velocity = RandomOrZero<DOFs>(dist(gen), 0.7);
            input.target_acceleration = RandomOrZero<DOFs>(dist(gen), 0.6);

            Vec max_v = limit_max * Vec::Random().array().abs() + Eigen::Map<Vec>(input.target_velocity.data()).array().abs() + limit_min;
            std::copy_n(max_v.data(), DOFs, input.max_velocity.begin());
            Vec max_a = limit_max * Vec::Random().array().abs() + Eigen::Map<Vec>(input.target_acceleration.data()).array().abs() + limit_min;
            std::copy_n(max_a.data(), DOFs, input.max_acceleration.begin());
            Vec max_j = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_j.data(), DOFs, input.max_jerk.begin());

            Vec t_v = Eigen::Map<Vec>(input.target_velocity.data());
            Vec t_a = Eigen::Map<Vec>(input.target_acceleration.data());
            auto max_target_acceleration = (2 * max_j.array() * (max_v.array() - t_v.array().abs())).sqrt();
            if ((t_a.array().abs() > max_target_acceleration.array()).any()) {
                continue;
            }

            check_calculation(otg, input);
        }
    }

#ifdef WITH_REFLEXXES
    SECTION("Comparison with Reflexxes with 1 DoF") {
        constexpr size_t DOFs {1};
        using Vec = Eigen::Matrix<double, DOFs, 1>;

        Ruckig<DOFs, true> otg {0.005};
        Reflexxes<DOFs> rflx {0.005};
        InputParameter<DOFs> input;

        srand(43);
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < (full ? 128 : 1) * 1024; ++i) {
            input.current_position = Random<DOFs>();
            input.current_velocity = RandomOrZero<DOFs>(dist(gen), 0.9);
            input.current_acceleration = RandomOrZero<DOFs>(dist(gen), 0.8);
            input.target_position = Random<DOFs>();
            input.target_velocity = RandomOrZero<DOFs>(dist(gen), 0.6);

            Vec max_v = limit_max * Vec::Random().array().abs() + Eigen::Map<Vec>(input.target_velocity.data()).array().abs() + limit_min;
            std::copy_n(max_v.data(), DOFs, input.max_velocity.begin());
            Vec max_a = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_a.data(), DOFs, input.max_acceleration.begin());
            Vec max_j = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_j.data(), DOFs, input.max_jerk.begin());

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
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < 1*1024; ++i) {
            input.current_position = Random<DOFs>();
            input.current_velocity = RandomOrZero<DOFs>(dist(gen), 0.9);
            input.current_acceleration = RandomOrZero<DOFs>(dist(gen), 0.8);
            input.target_position = Random<DOFs>();

            Vec max_v = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_v.data(), DOFs, input.max_velocity.begin());
            Vec max_a = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_a.data(), DOFs, input.max_acceleration.begin());
            Vec max_j = limit_max * Vec::Random().array().abs() + limit_min;
            std::copy_n(max_j.data(), DOFs, input.max_jerk.begin());

            check_comparison(otg, input, rflx);
        }
    } */
#endif
}
