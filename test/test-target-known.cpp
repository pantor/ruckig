#include <doctest/doctest.h>

#include <optional>

#include <ruckig/ruckig.hpp>


using namespace ruckig;


template<size_t DOFs, class OTGType>
void check_duration(OTGType& otg, InputParameter<DOFs>& input, double duration) {
    OutputParameter<DOFs> output;

    while (otg.update(input, output) == Result::Working) {
        input.current_position = output.new_position;
        input.current_velocity = output.new_velocity;
        input.current_acceleration = output.new_acceleration;
    }

    CHECK( output.trajectory.get_duration() == doctest::Approx(duration) );
}


TEST_CASE("known" * doctest::description("Known examples")) {
    RuckigThrow<3> otg {0.004};
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
    check_duration(otg, input, 1.493805);

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

    input.current_position = {0.02853333333333339, 0.0285, 0.0285};
    input.current_velocity = {0.6800000000000006, 0.68, 0.68};
    input.current_acceleration = {7.999999999999993, 8.0, 8.0};
    input.target_position = {0, 0, 0};
    input.target_velocity = {0, 0, 0};
    input.target_acceleration = {0, 0, 0};
    input.max_velocity = {1, 1 ,1};
    input.max_acceleration = {10, 10, 10};
    input.max_jerk = {100 + 1e-14, 100 + 1e-14, 100 + 1e-14};
    check_duration(otg, input, 0.58);

    input.current_position = {-0.05598571695553641, -0.534847776106059, 0.0978130731424748};
    input.current_velocity = {-0.03425673149926184, -0.8169926404190487, -0.004506245841081729};
    input.current_acceleration = {-2.720000000000001, 1.440254448401435, 0};
    input.target_position = {-0.0534691276550293, -0.6224863891601563, 0.09690575408935546};
    input.target_velocity = {0.0, 0.0, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {0.8500000000000001, 0.8500000000000001, 0.8500000000000001};
    input.max_acceleration = {4.25, 4.25, 4.25};
    input.max_jerk = {85.00000000000001, 85.00000000000001, 85.00000000000001};
    check_duration(otg, input, 0.2281604414);

    input.current_position = {0.0, 0.0, 0.3736320740840176};
    input.current_velocity = {0.0, 0.0, -0.60486324450823};
    input.current_acceleration = {0.0, 0.0, -0.4953501898933239};
    input.target_position = {0.0, 0.0, 0.233562911156468};
    input.target_velocity = {0.0, 0.0, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {1.0, 1.0, 10.01369296498101};
    input.max_acceleration = {1.0, 1.0, 14.72621077848741};
    input.max_jerk = {1.0, 1.0, 7.770133554060553};
    input.min_velocity = {-1.0, -1.0, -1.94898305867544};
    input.min_acceleration = {-1.0, -1.0, -0.6829625196960336};
    check_duration(otg, input, 1.08732372);

    input.current_position = {-0.01919986582215404, 0.0, 0.0};
    input.current_velocity = {-0.3858205249368821, 0.0, 0.0};
    input.current_acceleration = {0.1889847091893647, 0.0, 0.0};
    input.target_position = {1.297187158009963, 0.0, 0.0};
    input.target_velocity = {1.160424379732321, 0.0, 0.0};
    input.target_acceleration = {0.4552736879206988, 0.0, 0.0};
    input.max_velocity = {14.64378197125325, 3.0, 3.0};
    input.max_acceleration = {0.4552736879216988, 2.5, 2.5};
    input.max_jerk = {12.15045820314999, 2.2, 2.2};
    input.minimum_duration = 3.408914;
    check_duration(otg, input, 3.408914);

    input.current_position = {0.2473592757796861, 0.2921606775204735, 0.7758663276711127};
    input.current_velocity = {-0.2426115138900957, 0.2200706500820608, -0.01891492763905089};
    input.current_acceleration = {0.01877538437863763, -0.6573642866096158, 0.5};
    input.target_position = {0.2308075286416321, 0.3066442218484541, 0.7733155040940536};
    input.target_velocity = {0.0, 0.0, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {0.25, 0.25, 0.25};
    input.max_acceleration = {2.5, 2.5, 2.5};
    input.max_jerk = {62.5, 62.5, 62.5};
    input.min_velocity = std::nullopt;
    input.min_acceleration = std::nullopt;
    input.minimum_duration = std::nullopt;
    input.duration_discretization = DurationDiscretization::Discrete;
    check_duration(otg, input, 0.14);

    input.current_position = {0.5289912019692077, -0.2461593579591288, -0.2728396804501142};
    input.current_velocity = {0.0287779218983349, -0.005980397028399779, 0.04763314105835294};
    input.current_acceleration = {-1.374657509445135, 0.2564786317266955, -2.497210506808601};
    input.target_position = {0, 0, 0};
    input.target_velocity = {0, 0, 0};
    input.target_acceleration = {0, 0, 0};
    input.max_velocity = {6.457718232379019, 5.410520681182422, 7.155849933176751};
    input.max_acceleration = {26.89552377323262, 22.5322006432468, 29.81022362406315};
    input.max_jerk = {224.2224490037115, 187.8497873921497, 248.4650723139127};
    input.min_velocity = std::nullopt;
    input.min_acceleration = std::nullopt;
    input.minimum_duration = std::nullopt;
    input.duration_discretization = DurationDiscretization::Discrete;
    input.control_interface = ControlInterface::Velocity;
    check_duration(otg, input, 0.024);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.1119477497536703, -0.005706738140158095, 0.0};
    input.current_acceleration = {-2.943871184141059, 0.1638588832925878, 0.0};
    input.target_position = {0, 0, 0};
    input.target_velocity = {0, 0, 0};
    input.target_acceleration = {0, 0, 0};
    input.max_velocity = {10.0, 10.0, 10.0};
    input.max_acceleration = {10.0, 10.0, 10.0};
    input.max_jerk = {224.2224490037115, 187.8497873921497, 10.0};
    input.min_velocity = std::nullopt;
    input.min_acceleration = std::nullopt;
    input.minimum_duration = std::nullopt;
    input.duration_discretization = DurationDiscretization::Continuous;
    input.control_interface = ControlInterface::Velocity;
    check_duration(otg, input, 0.0352632);

    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {1.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};
    input.target_position = {0.0, 0.0, 0.0};
    input.target_velocity = {1.0, 0.0, 1.0};
    input.target_acceleration = {0.0, 0.0, 0.0};
    input.max_velocity = {1.0, 1.0, 1.0};
    input.max_acceleration = {1.0, 1.0, 1.0};
    input.max_jerk = {1.0, 1.0, 1.0};
    input.control_interface = ControlInterface::Velocity;
    check_duration(otg, input, 2.0);

    input.current_position = {-0.54516231864478149, -1.4206629551476477, 0.32727720821160067};
    input.current_velocity = {6.1348323207600686, -0.0040477518609957240, 0.0088178054233364854};
    input.current_acceleration = {0.0000000000000000, 0.19725113486681556, -0.23382835447542066};
    input.target_velocity = {6.1348323207600686, 0.014259490365703332, -0.012798219033181931};
    input.max_velocity = {6.4577182323790190, 5.4105206811824216, 7.1558499331767509};
    input.max_acceleration = {26.895523773232618, 22.532200643246796, 29.810223624063148};
    input.max_jerk = {224.22244900371152, 187.84978739214969, 248.46507231391274};
    input.control_interface = ControlInterface::Velocity;
    check_duration(otg, input, 0.0187497625);

    input.current_position = {-1.2563817016634644, -0.51079124473334669, 2.4607439315667303};
    input.current_velocity = {3.4194572623820010, 5.7548168684509085, -1.5675191369761612};
    input.current_acceleration = {27.722541632719086, -1.7945991957745204, -22.971380686605098};
    input.target_velocity = {5.1559110036108979, 5.5636136585722848, -3.6173923596604807};
    input.max_velocity = {9.5993108859688121, 9.4247779607693793, 17.453292519943297};
    input.max_acceleration = {39.985493163190093, 39.618974020271281, 72.710416638083771};
    input.max_jerk = {333.30552725335713, 330.26865435488696, 606.01322287747109};
    input.control_interface = ControlInterface::Velocity;
    check_duration(otg, input, 0.1030382161);

    input.current_position = {0, 0, 0};
    input.current_velocity = {0, -0.0577630321017372216, 0};
    input.current_acceleration = {0, -1.9198621771937627, 0};
    input.target_position = {0, 0, 0};
    input.target_velocity = {0, -0.1, 0};
    input.target_acceleration = {0, 0, 0};
    input.max_velocity = {5.23, 3.926990816987241394, 3.92};
    input.max_acceleration = {4.36, 4.363323129985824, 4.36};
    input.max_jerk = {43.63, 43.633231299858238116, 43.63};
    input.control_interface = ControlInterface::Velocity;
    input.duration_discretization = ruckig::DurationDiscretization::Discrete;
    check_duration(otg, input, 0.048);

    input.current_position = {-19.93333333333424, -0.4983333333333563, 0};
    input.current_velocity = {10, 0.25, 0};
    input.current_acceleration = {0, 0, 0};
    input.target_position = {20, 0.5, 0.2};
    input.target_velocity = {10, 0.25, 0};
    input.target_acceleration = {0, 0, 0};
    input.max_velocity = {10, 10, 10};
    input.max_acceleration = {15, 15, 15};
    input.max_jerk = {15, 5, 2};
    input.control_interface = ControlInterface::Position;
    input.duration_discretization = ruckig::DurationDiscretization::Continuous;

    Trajectory<3> traj;
    otg.calculate(input, traj);
    check_duration(otg, input, 3.99333);

    std::array<double, 3> new_position, new_velocity, new_acceleration;
    traj.at_time(traj.get_duration() / 2, new_position, new_velocity, new_acceleration);
    CHECK( new_position[0] == doctest::Approx(0.033333) );
    CHECK( new_position[1] == doctest::Approx(0.000833) );


    RuckigThrow<38> otg38 {0.004};
    InputParameter<38> input38;
    input38.current_position = {0.5, -0.7, -0.7, -1.5, -0, -0.2, 0.5, 0.7, 0.7, -1.5, 0, -0.2, -0, 0, -0, 0, 0, 0, 0, -0, 0, -0, 0, 0, 0, -0, 0, -0, 0, 0, -0, 0, -0, 0, 0, 0, -0, 0};
    input38.current_velocity = {2.689534009085704e-06, -2.689534009088848e-06, -2.689534009086246e-06, -2.689534009088414e-06, 2.689534009087655e-06, -2.689534009089241e-06, 2.68953400908874e-06, 2.689534009087113e-06, 2.689534009086246e-06, -2.689534009088848e-06, -2.689534009089065e-06, -2.689534009089228e-06, 2.689534009088198e-06, -2.68953400908668e-06, 2.689534009089715e-06, -2.689534009087113e-06, -2.4911206284336e-12, 0, -2.689534009089065e-06, 2.689534009089065e-06, -2.689534009087764e-06, 2.689534009059358e-06, -2.4911206284336e-12, 0, -2.689534009088957e-06, 2.689534009088631e-06, -2.689534009088198e-06, 2.689534009031602e-06, -2.4911206284336e-12, 0, 2.689534009088848e-06, -2.689534009088414e-06, 2.689534009086463e-06, -2.689534009142625e-06, -2.4911206284336e-12, 0, -2.689534009086869e-06, -2.689534009086951e-06};
    input38.current_acceleration = {-0.001639979880695402, 0.001639979880696266, 0.001639979880695815, 0.00163997988069603, -0.00163997988069612, 0.001639979880696489, -0.001639979880696373, -0.001639979880695849, -0.001639979880695794, 0.001639979880696522, 0.001639979880696495, 0.001639979880696486, -0.00163997988069612, 0.001639979880695946, -0.001639979880696397, 0.001639979880696085, 0, 0, 0.001639979880696328, -0.001639979880696418, 0.001639979880696141, -0.001639979880696363, 0, 0, 0.001639979880696443, -0.001639979880696245, 0.001639979880696175, -0.001639979880695752, 0, 0, -0.001639979880696356, 0.001639979880696342, -0.001639979880695766, 0.001639979880696307, 0, 0, 0.001639979880695763, 0.001639979880695815};
    input38.target_position = {1.5, -0.7, -0.7, -1.5, 0, -0.2, 1.5, 0.7, 0.7, -1.5, 0, -0.2, 0, 1.7, 2.3, 0.6, 0, 0, 0, -1.7, -2.3, -0.6, 0, 0, 0, -1.7, -2.3, -0.6, 0, 0, 0, 1.7, 2.3, 0.6, 0, 0, 0, 0};
    input38.target_velocity = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    input38.target_acceleration = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    input38.max_velocity = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    input38.max_acceleration = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    input38.max_jerk = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    check_duration(otg38, input38, 6.60328);


    RuckigThrow<1> otg1 {0.01};
    InputParameter<1> input1;
    input1.current_position = {0.0};
    input1.current_velocity = {1.0};
    input1.current_acceleration = {-2.0};
    input1.target_position = {1.0};
    input1.target_velocity = {0.5};
    input1.target_acceleration = {0.0};
    input1.max_velocity = {0.8};
    input1.max_acceleration = {2};
    input1.max_jerk = {5};

    check_duration(otg1, input1, 1.6041);


    RuckigThrow<6> otg6 {0.004};
    InputParameter<6> input6;
    input6.current_position = {0.1720894642466409, 3.62255615705579, 1.675754445058354, 2.181127634157466, 0.7829547024446563, 0};
    input6.current_velocity = {0, -0.06422550797049381, -0.006910864193129071, 0.1207111377965507, 0, 0};
    input6.current_acceleration = {0, 0, 0, -0.3999999999978001, 0, 0};
    input6.target_position = {0.1703441349946466, 3.587689540855463, 1.671344092861895, 2.271916653384632, 0.7826056365942573, -0.003141592653589793};
    input6.target_velocity = {0, -0.06423623454575711, -0.007006078877926925, 0.1187909578397501, 0, 0};
    input6.target_acceleration = {0, 0, 0, 0, 0, 0};
    input6.max_velocity = {0.06981316999999999, 0.06423623454575711, 0.007006078877926925, 0.1187909578397501, 6.981317008, 34.90658504};
    input6.max_acceleration = {0.034906585, 0.261799388, 0.043633231, 0.4, 4.01425728, 34.90658504};
    input6.max_jerk = {100, 100, 100, 100, 100, 100};
    check_duration(otg6, input6, 0.764274);

    input6.duration_discretization = DurationDiscretization::Discrete;
    input6.current_position = { 0.1174161453909641,-0.4222080077114712,0.3003863125814782,2.960444926889819,-0.005791262234199017,-2.364847659719411 };
    input6.current_velocity = { -3.662950284795513e-09,-5.585273134233322e-09,0.01325125623820627,1.800912308381978e-08,-1.307199431356758e-09,-1.908875705031016e-08 };
    input6.current_acceleration = { 1.350015487635199e-07,2.058506025917496e-07,-0.4883877683641483,-6.637435179724953e-07,4.81780898060446e-08,7.035344641338199e-07 };
    input6.target_position = { 0.117416145324707,-0.4222080078125,0.3006260070800781,2.960444927215576,-0.00579126225784421,-2.364847660064697 };
    input6.target_velocity = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    input6.target_acceleration = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    input6.max_velocity = { 0.002,0.002,0.002,0.003141592629253865,0.003141592629253865,0.003141592629253865 };
    input6.max_acceleration = { 0.01,0.01,0.01,0.01570796314626932,0.01570796314626932,0.01570796314626932 };
    input6.max_jerk = { 0.18,0.18,0.18,1.570796314626932,1.570796314626932,1.570796314626932 };
    check_duration(otg6, input6, 11198.4);
}
