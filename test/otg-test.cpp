#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <random>
#include "randomizer.hpp"

#include <ruckig/ruckig.hpp>

#ifdef WITH_REFLEXXES
#include <ruckig/reflexxes_comparison.hpp>
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

    CHECK( output.trajectory.get_duration() == doctest::Approx(duration) );
}


template<size_t DOFs, class OTGType>
inline void check_calculation(OTGType& otg, InputParameter<DOFs>& input) {
    OutputParameter<DOFs> output;

    CAPTURE( input );

    auto result = otg.update(input, output);
    if (result == Result::ErrorTrajectoryDuration) {
        return;
    }

    CHECK( (result == Result::Working || (result == Result::Finished && output.trajectory.get_duration() < 0.005)) );
    CHECK( output.trajectory.get_duration() >= 0.0 );

    for (size_t dof = 0; dof < otg.degrees_of_freedom; ++dof) {
        CHECK_FALSE( (std::isnan(output.new_position[dof]) || std::isnan(output.new_velocity[dof]) || std::isnan(output.new_acceleration[dof])) );
    }
}


template<size_t DOFs, class OTGType>
inline size_t step_through_and_check_calculation(OTGType& otg, InputParameter<DOFs>& input, size_t max_number_checks) {
    OutputParameter<DOFs> output;
    OTGType otg_second {0.001};

    check_calculation<DOFs, OTGType>(otg, input);

    size_t number_checks {1};
    while (otg.update(input, output) == Result::Working) {
        input.current_position = output.new_position;
        input.current_velocity = output.new_velocity;
        input.current_acceleration = output.new_acceleration;

        // Or randomize input with small noise here?
        check_calculation<DOFs, OTGType>(otg_second, input);

        number_checks += 1;
        if (number_checks == max_number_checks) {
            break;
        }
    }

    return number_checks;
}


template<size_t DOFs, class OTGType, class OTGCompType>
inline void check_comparison(OTGType& otg, InputParameter<DOFs>& input, OTGCompType& otg_comparison) {
    OutputParameter<DOFs> output;

    CAPTURE( input );

    auto result = otg.update(input, output);
    if (result == Result::ErrorTrajectoryDuration) {
        return;
    }

    CHECK( (result == Result::Working || (result == Result::Finished && output.trajectory.get_duration() < 0.005)) );

    OutputParameter<DOFs> output_comparison;
    auto result_comparison = otg_comparison.update(input, output_comparison);
    CHECK( output.trajectory.get_duration() <= doctest::Approx(output_comparison.trajectory.get_duration()) );
}


int seed {42};
size_t number_trajectories {120000}; // Some user variable you want to be able to set
size_t random_1, random_3, step_through_3, random_discrete_3, random_direction_3, comparison_1, comparison_3, velocity_random_3;

std::normal_distribution<double> position_dist {0.0, 4.0};
std::normal_distribution<double> dynamic_dist {0.0, 0.8};
std::uniform_real_distribution<double> limit_dist {0.08, 16.0};
// std::uniform_real_distribution<double> limit_dist {800.0, 160000.0};
std::uniform_real_distribution<double> min_limit_dist {-16.0, -0.08};


template<class T>
void check_array(const T& first, const T& second) {
    for (size_t dof = 0; dof < first.size(); ++dof) {
        CHECK( first[dof] == doctest::Approx(second[dof]) );
    }
}


TEST_CASE("secondary" * doctest::description("Secondary Features")) {
    Ruckig<3, true> otg {0.005};
    InputParameter<3> input;
    OutputParameter<3> output;

    input.current_position = {0.0, -2.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {1.0, -3.0, 2.0};
    input.target_velocity = {0.0, 0.3, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {1.0, 1.0, 1.0};
    input.max_acceleration = {1.0, 1.0, 1.0};
    input.max_jerk = {1.0, 1.0, 1.0};

    auto result = otg.update(input, output);

    CHECK( result == Result::Working );
    CHECK( output.trajectory.get_duration() == doctest::Approx(4.0) );

    std::array<double, 3> new_position, new_velocity, new_acceleration;
    output.trajectory.at_time(0.0, new_position, new_velocity, new_acceleration);
    check_array(new_position, input.current_position);
    check_array(new_velocity, input.current_velocity);
    check_array(new_acceleration, input.current_acceleration);

    output.trajectory.at_time(output.trajectory.get_duration(), new_position, new_velocity, new_acceleration);
    check_array(new_position, input.target_position);
    check_array(new_velocity, input.target_velocity);
    check_array(new_acceleration, input.target_acceleration);

    output.trajectory.at_time(2.0, new_position, new_velocity, new_acceleration);
    check_array(new_position, {0.5, -2.6871268303, 1.0});

    auto independent_min_durations = output.trajectory.get_independent_min_durations();
    CHECK( independent_min_durations[0] == doctest::Approx(3.1748021039) );
    CHECK( independent_min_durations[1] == doctest::Approx(3.6860977315) );
    CHECK( independent_min_durations[2] == doctest::Approx(output.trajectory.get_duration()) );

    auto position_extrema = output.trajectory.get_position_extrema();
    CHECK( position_extrema[0].t_max == doctest::Approx(4.0) );
    CHECK( position_extrema[0].max == doctest::Approx(1.0) );
    CHECK( position_extrema[0].t_min == doctest::Approx(0.0) );
    CHECK( position_extrema[0].min == doctest::Approx(0.0) );

    CHECK( position_extrema[1].t_max == doctest::Approx(0.0) );
    CHECK( position_extrema[1].max == doctest::Approx(-2.0) );
    CHECK( position_extrema[1].t_min == doctest::Approx(3.2254033308) );
    CHECK( position_extrema[1].min == doctest::Approx(-3.1549193338) );

    CHECK( position_extrema[2].t_max == doctest::Approx(4.0) );
    CHECK( position_extrema[2].max == doctest::Approx(2.0) );
    CHECK( position_extrema[2].t_min == doctest::Approx(0.0) );
    CHECK( position_extrema[2].min == doctest::Approx(0.0) );


    double time;
    CHECK( output.trajectory.get_first_time_at_position(0, 0.0, time) );
    CHECK( time == doctest::Approx(0.0) );

    CHECK( output.trajectory.get_first_time_at_position(0, 0.5, time) );
    CHECK( time == doctest::Approx(2.0) );

    CHECK( output.trajectory.get_first_time_at_position(0, 1.0, time) );
    CHECK( time == doctest::Approx(4.0) );

    CHECK( output.trajectory.get_first_time_at_position(1, -3.0, time) );
    CHECK( time == doctest::Approx(2.6004877902) );

    CHECK( output.trajectory.get_first_time_at_position(1, -3.1, time) );
    CHECK( time == doctest::Approx(2.8644154489) );

    CHECK( output.trajectory.get_first_time_at_position(2, 0.05, time) );
    CHECK( time == doctest::Approx(0.6694329501) );

    CHECK_FALSE( output.trajectory.get_first_time_at_position(0, -1.0, time) );
    CHECK_FALSE( output.trajectory.get_first_time_at_position(1, -3.4, time) );
    CHECK_FALSE( output.trajectory.get_first_time_at_position(6, 0.0, time) );


    input.current_position = {0.0, -2.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {1.0, -3.0, 2.0};
    input.target_velocity = {2.0, 0.3, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {1.0, 1.0, 1.0};
    input.max_acceleration = {1.0, 1.0, 1.0};
    input.max_jerk = {1.0, 1.0, 1.0};

    result = otg.update(input, output);

    CHECK( result == Result::ErrorInvalidInput );
    CHECK_FALSE( output.new_calculation );

    input.target_velocity = {0.2, -0.3, 0.8};
    result = otg.update(input, output);

    CHECK( result == Result::Working );
    CHECK( output.new_calculation );
}

TEST_CASE("dynamic-dofs" * doctest::description("Dynamic DoFs")) {
    Ruckig<0, true> otg {3, 0.005};
    InputParameter<0> input {3};
    OutputParameter<0> output {3};

    input.current_position = {0.0, -2.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {1.0, -3.0, 2.0};
    input.target_velocity = {0.0, 0.3, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {1.0, 1.0, 1.0};
    input.max_acceleration = {1.0, 1.0, 1.0};
    input.max_jerk = {1.0, 1.0, 1.0};

    auto result = otg.update(input, output);

    CHECK( result == Result::Working );
    CHECK( output.trajectory.get_duration() == doctest::Approx(4.0) );

    std::vector<double> new_position(3), new_velocity(3), new_acceleration(3);
    output.trajectory.at_time(0.0, new_position, new_velocity, new_acceleration);
    check_array(new_position, input.current_position);
    check_array(new_velocity, input.current_velocity);
    check_array(new_acceleration, input.current_acceleration);
}

TEST_CASE("known" * doctest::description("Known examples")) {
    Ruckig<3, true> otg {0.005};

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

    input.target_position = {0.0, 1e-17, -1e-17};
    check_duration(otg, input, 1e-18);

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

    input.current_position = {0.3888899206957 - 10e-14, 0.0, 0.0};
    input.current_velocity = {0.2231429352410215, 0.0, 0.0};
    input.current_acceleration = {-0.2987593916455, 0.0, 0.0};
    input.target_position = {0.5, 0.0, 0.0};
    input.target_velocity = {0.0, 0.0, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {3.0, 3.0, 3.0};
    input.max_acceleration = {0.5, 0.5, 0.5};
    input.max_jerk = {0.2, 0.2, 0.2};
    check_duration(otg, input, 1.4939456041);

    input.current_position = {-5.54640573838539, -2.34195463203842, 5.10070661762967};
    input.current_velocity = {0.824843228617216, -1.03863337183304, -0.749451523227729};
    input.current_acceleration = {-0.119403564898501, 0.923861820607788, 3.04022341347259};
    input.target_position = {-1.58293112753888, 0.383405919465141, 5.79349604610299};
    input.target_velocity = {-1.59453676324393, 0.0, -0.0693173526513803};
    input.target_acceleration = {-0.664429703711622, 0.0, 0.0};
    input.max_velocity = {12.9892953062198, 3.74169932927481, 1.42398447457303};
    input.max_acceleration = {4.2162106624246, 10.2906731766853, 2.1869079548297};
    input.max_jerk = {8.03496976453435, 0.200684346397485 - 1.0e-14, 0.0848503482861296};
    check_duration(otg, input, 1921.0797627836);

    input.current_position = {-6.49539540831446, 6.14883133273172, -2.02636240900911};
    input.current_velocity = {-1.14327601654428, 0.00991019970085593, -1.00932863927626};
    input.current_acceleration = {-1.73501068960131, -0.584885092422228, 0};
    input.target_position = {4.4676187540058, 2.93367894961155, -0.646008452514058};
    input.target_velocity = {-0.544559915133859, 0.298517792372943, 1.6058847848484};
    input.target_acceleration = {-1.31832055647831, 0, 0};
    input.max_velocity = {8.65978706670502, 5.94921088330542, 10.7652253566829};
    input.max_acceleration = {3.40137210377608, 4.04166318018487, 10.8617860610581};
    input.max_jerk = {10.9542353113865, 3.11056302676629, 0.798055744482636 + 9e-12};
    check_duration(otg, input, 4.6277455678);

    input.current_position = {7.06378251402596, -2.4834697862831, -0.843847405371359};
    input.current_velocity = {0.436985859305842, 0.0708113515655622, -0.751266816040307};
    input.current_acceleration = {-0.80835350359544, 0, -0.355284934641626};
    input.target_position = {4.40606827118048, -2.84629921001043, -2.91829890043522};
    input.target_velocity = {0.555084596169823, 0, -1.24631524923535};
    input.target_acceleration = {0.463000173872542, 0, 0};
    input.max_velocity = {7.97399137456765, 2.68591430972239, 9.54987666746364};
    input.max_acceleration = {5.44679859206862, 7.61752909348119, 0.473482772614085};
    input.max_jerk = {7.88958080921515, 5.26855927512131, 0.764061581326592 - 1e-14};
    check_duration(otg, input, 8.8739464323);

    input.current_position = {-7.962737259350095, 0, 0};
    input.current_velocity = {-0.8844863500141733, 0, 0};
    input.current_acceleration = {2.252932547031004, 0, 0};
    input.target_position = {-3.547368989678775, 0, 0};
    input.target_velocity = {0, 0, 0};
    input.target_acceleration = {0.217242176687843, 0, 0};
    input.max_velocity = {0.1241065584614779, 1, 1};
    input.max_acceleration = {1.808598147153279, 1, 1};
    input.max_jerk = {2.516849090900998 - 3e-15, 1, 1};
    check_duration(otg, input, 38.3409477609);

    input.current_position = {-4.180150148354134, 1.030371049895473, -2.660154279239869};
    input.current_velocity = {1.673805463302308, -1.435796222257198, 0.9711306630275642};
    input.current_acceleration = {1.412175048500792, 1.892262449040863, -1.128847905860926};
    input.target_position = {2.079913937916431, 1.839862681333277, 2.341421542126605};
    input.target_velocity = {0.7537566830764975, 0, 0.02507782261105568};
    input.target_acceleration = {-0.8610296259045267, -0.07876324073516261, 0};
    input.max_velocity = {1.863775561344568, 0.4357836109021987, 6.260907804906162};
    input.max_acceleration = {9.49223908896113, 9.002562577262177, 1.119142029086944};
    input.max_jerk = {8.689575453772798, 0.09322235504216797, 0.1594452521517275 + 3e-15};
    check_duration(otg, input, 1135.0135089249);

    input.current_position = {-4.490717417930574, 3.467236624628543, -0.7545929089757601};
    input.current_velocity = {0.1839756723363622, -0.4356283320280516, 0.7490399525818022};
    input.current_acceleration = {-1.057769973808928, 0, -2.368645439140517};
    input.target_position = {-4.928244836531066, -4.821780824003112, -8.20567952461017};
    input.target_velocity = {0.1097319156272965, -0.9272874846270881, 0};
    input.target_acceleration = {0.03089046366221739, -0.9744054582899561, 0};
    input.max_velocity = {6.144314006624488, 2.93258338415229, 0.1820021269527196};
    input.max_acceleration = {5.199401036221791, 1.848176490768948, 11.11168017805234};
    input.max_jerk = {9.940940357283978, 10.46997753899755, 0.08166297169205029};
    check_duration(otg, input, 7295.4375633935);

    input.current_position = {0.01073568005271233, -0.7002627264230739, 0};
    input.current_velocity = {0.05656281587106524, 1.011281770884991, 0};
    input.current_acceleration = {-5.348847133445708, -3.400994300842285, 0};
    input.target_position = {0.0698, 0.6283, 0};
    input.target_velocity = {0, 0, 0};
    input.target_acceleration = {0, 0, 0};
    input.max_velocity = {1, 1, 1};
    input.max_acceleration = {7, 7, 7};
    input.max_jerk = {1000, 1000, 1000};
    check_duration(otg, input, 1.403613276);

    input.current_position = {0.0001215, 0, 0};
    input.current_velocity = {0.00405, 0, 0};
    input.current_acceleration = {0.09, 0, 0};
    input.target_position = {0.1421083333333333087, 0, 0};
    input.target_velocity = {0.37, 0, 0};
    input.target_acceleration = {0.5, 0, 0};
    input.max_velocity = {1, 1, 1};
    input.max_acceleration = {0.5, 0.5, 0.5};
    input.max_jerk = {1, 1, 1};
    check_duration(otg, input, 0.9);

    input.current_position = {0, 0, 0};
    input.current_velocity = {0, 0, 0};
    input.current_acceleration = {0, 0, 0};
    input.target_position = {400, 4000, 40000};
    input.target_velocity = {0, 0, 0};
    input.target_acceleration = {0, 0, 0};
    input.max_velocity = {1800, 18000, 180000};
    input.max_acceleration = {20000, 200000, 2000000};
    input.max_jerk = {200000, 2000000, 20000000};
    check_duration(otg, input, 0.4119588818);
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
        if (i < random_discrete_3 / 2) {
            input.synchronization = Synchronization::Phase;
        } else {
            input.synchronization = Synchronization::Time;
        }

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

    input.current_position = {0.0, 0.0, 0.0};

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
        if (i < random_3 / 2) {
            input.synchronization = Synchronization::Phase;
        } else {
            input.synchronization = Synchronization::Time;
        }

        if (i < random_3 / 20) {
            input.duration_discretization = DurationDiscretization::Discrete;
        } else {
            input.duration_discretization = DurationDiscretization::Continuous;
        }

        p.fill(input.current_position);
        d.fill_or_zero(input.current_velocity, 0.9);
        d.fill_or_zero(input.current_acceleration, 0.8);
        p.fill(input.target_position);
        d.fill_or_zero(input.target_velocity, 0.7);
        d.fill_or_zero(input.target_acceleration, 0.6);
        l.fill(input.max_velocity, input.target_velocity);
        l.fill(input.max_acceleration, input.target_acceleration);
        l.fill(input.max_jerk);

        // if (i % 1000 == 0) {
        //     const double factor = 1e3;
        //     for (size_t d = 0; d < DOFs; ++d) {
        //         input.current_position[d] *= factor;
        //         input.current_velocity[d] *= factor;
        //         input.current_acceleration[d] *= factor;
        //         input.target_position[d] *= factor;
        //         input.target_velocity[d] *= factor;
        //         input.target_acceleration[d] *= factor;
        //         input.max_velocity[d] *= factor;
        //         input.max_acceleration[d] *= factor;
        //         input.max_jerk[d] *= factor;
        //     }
        // }

        if (!otg.validate_input(input)) {
            --i;
            continue;
        }

        check_calculation(otg, input);
    }
}

TEST_CASE("step_through_3" * doctest::description("Step through random input with 3 DoF and target velocity, acceleration")) {
    constexpr size_t DOFs {3};
    Ruckig<DOFs, true> otg {0.01};
    InputParameter<DOFs> input;

    Randomizer<DOFs, decltype(position_dist)> p { position_dist, seed + 3 };
    Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist, seed + 4 };
    Randomizer<DOFs, decltype(limit_dist)> l { limit_dist, seed + 5 };

    for (size_t i = 0; i < step_through_3; ++i) {
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

        i += step_through_and_check_calculation(otg, input, 1000);
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

    // To "activate" std::optional
    input.min_velocity = {0.0, 0.0, 0.0};
    input.min_acceleration = {0.0, 0.0, 0.0};

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
    step_through_3 = 0; // number_trajectories / 20;
    random_direction_3 = number_trajectories / 50;
    velocity_random_3 = number_trajectories / 10;
    random_3 = number_trajectories - (random_1 + step_through_3 + random_direction_3 + comparison_1 + comparison_3 + velocity_random_3 + random_discrete_3);
    std::cout << "<number_trajectories> Random (1 DoF): " << random_1 << " (3 DoF): " << random_3 << "  Step Through (3 Dof): " << step_through_3 << "  Comparison (1 DoF): " << comparison_1 << " (3 DoF): " << comparison_3 << " Total: " << number_trajectories << std::endl;

    return context.run();
}
