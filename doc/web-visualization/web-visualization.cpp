#include <iostream>
#include <sstream>
#include <emscripten.h>
#include <ruckig/ruckig.hpp>
#include <json/json.hpp>

using namespace ruckig;

EM_JS(void, NoReturnValueWithStringParameter, (const char* string_pointer), {
  updateGraph(Module.UTF8ToString(string_pointer));
})

extern "C" {
const char* computeCurve(float current_position, float current_velocity, float current_acceleration,
                  float target_position, float target_velocity, float target_acceleration,
                  float max_velocity, float max_jerk, float max_acceleration) {
    // Create instances: the Ruckig OTG as well as input and output parameters
    Ruckig<1> otg{0.01};  // control cycle
    InputParameter<1> input;
    OutputParameter<1> output;
    unsigned int i = 0;
    nlohmann::json json;
    nlohmann::json dataArray = nlohmann::json::array();
    nlohmann::json formatArray = nlohmann::json::array();

    input.control_interface = ControlInterface::Position;

    input.current_position = {current_position};
    input.current_velocity = {current_velocity};
    input.current_acceleration = {current_acceleration};

    input.target_position = {target_position};
    input.target_velocity = {target_velocity};
    input.target_acceleration = {target_acceleration};

    input.max_velocity = {max_velocity};
    input.max_jerk = {max_jerk};
    input.max_acceleration = {max_acceleration};

//    std::stringstream ss;
    static std::string errorMessage;

    std::cout << "Validating input" << std::endl;
    try {
        bool ret = otg.validate_input(input, false, true);
        if (ret == false) {
            std::stringstream error;
            error << "input validation error, can't continue!" << std::endl;
            std::cout << error.str();
            errorMessage = error.str();
            return errorMessage.c_str();
        }
    } catch (const ruckig::RuckigError& e) {
        std::stringstream error;
        error << "input validation RuckigError exception occurred: " << e.what() << std::endl;
        std::cout << error.str();
        errorMessage = error.str();
        return errorMessage.c_str();
    } catch (...) {
        std::stringstream error;
        error << "input validation unknown exception occurred: " << std::endl;
        std::cout << error.str();
        errorMessage = error.str();
        return errorMessage.c_str();
    }
    std::cout << "Validating complete, all good! Computing otg.update" << std::endl;
    while (true) {
        Result result;
        try {
            result = otg.update(input, output);
        } catch (const ruckig::RuckigError& e) {
            std::stringstream error;
            error << "input validation RuckigError exception occurred: " << e.what() << std::endl;
            std::cout << error.str();
            errorMessage = error.str();
            return errorMessage.c_str();
        } catch (...) {
            std::stringstream error;
            error << "input validation unknown exception occurred: " << std::endl;
            std::cout << error.str();
            errorMessage = error.str();
            return errorMessage.c_str();
        }
        json["trajectory_duration"] = output.trajectory.get_duration();
        if (i == 0) {
            std::cout << "Trajectory duration: " << output.trajectory.get_duration() << " [s]." << std::endl;
            std::cout << "degrees_of_freedom: " << output.trajectory.degrees_of_freedom << std::endl;
            std::cout << "output.trajectory.get_independent_min_durations(): " << std::endl;
            for (int x=0; x < output.trajectory.get_independent_min_durations().size(); x++) {
                std::cout << output.trajectory.get_independent_min_durations().at(x) << std::endl;
            }
            std::cout << "output.trajectory.get_intermediate_durations(): " << std::endl;
            for (int x=0; x < output.trajectory.get_intermediate_durations().size(); x++) {
                std::cout << output.trajectory.get_intermediate_durations().at(x) << std::endl;
            }
            std::cout << "output.trajectory.get_position_extrema(): " << std::endl;
            for (int x=0; x < output.trajectory.get_position_extrema().size(); x++) {
                std::cout << "min: " << output.trajectory.get_position_extrema().at(x).min << std::endl;
                std::cout << "max: " << output.trajectory.get_position_extrema().at(x).max << std::endl;
                std::cout << "t_max: " << output.trajectory.get_position_extrema().at(x).t_max << std::endl;
                std::cout << "t_min: " << output.trajectory.get_position_extrema().at(x).t_min << std::endl;
            }
            std::cout << "output.trajectory.get_profiles(): " << std::endl;
            for (int x=0; x < output.trajectory.get_profiles().size(); x++) {
                for (int y=0; y < output.trajectory.get_profiles().size(); y++) {
                    for (int z=0; z < output.trajectory.get_profiles().at(x).at(y).j.size(); z++) {
                        std::cout << j0(output.trajectory.get_profiles().at(x).at(y).j.at(z)) << std::endl;
                    }
                }
            }
        }
        if (result == Result::Working) {
            auto &p = output.new_position;
            auto &a = output.new_acceleration;
            auto &v = output.new_velocity;
            auto &j = output.new_jerk;
            if (i == 0) {
                formatArray.push_back("time");
                formatArray.push_back("position");
                formatArray.push_back("velocity");
                formatArray.push_back("acceleration");
                formatArray.push_back( "jerk");
            }
            nlohmann::json entry = nlohmann::json::array();
            entry.push_back(output.time);
            entry.push_back(p[0]);
            entry.push_back(v[0]);
            entry.push_back(a[0]);
            entry.push_back(j[0]);
            dataArray.push_back(entry);
            //std::cout << output.trajectory << std::endl;
            output.pass_to_input(input);
        } else if (result == Result::Finished) {
            std::cout << "Finished calculation, js callback!" << std::endl;
            json["curve"]["data"]=dataArray;
            json["curve"]["format"]=formatArray;
//            std::cout << json.dump() << std::endl;
            NoReturnValueWithStringParameter(json.dump().c_str());
            break;
        } else {
            std::stringstream error;
            switch(result) {
                case 0:
                    error << "Working";
                    break;
                case 1:
                    error << "Finished";;
                    break;
                case -1:
                    error << "Error";
                    break;
                case -100:
                    error << "ErrorInvalidInput";
                    break;
                case -101:
                    error << "ErrorTrajectoryDuration";
                    break;
                case -102:
                    error << "ErrorPositionalLimits";
                    break;
                case -104:
                    error << "ErrorZeroLimits";
                    break;
                case -110:
                    error << "ErrorExecutionTimeCalculation";
                    break;
                case -111:
                    error << "ErrorSynchronizationCalculation";
                    break;
            }
            std::cout << error.str();
            errorMessage = error.str();
            return errorMessage.c_str();
        }
        i += 1;
    }
    return NULL;
}
}