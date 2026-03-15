// Only with Ruckig Pro

#include <ruckig/ruckig.hpp>

#include "plotter.hpp"


using namespace ruckig;

int main() {
    // Create instances: the Ruckig trajectory generator as well as input and output parameters
    Ruckig<3> ruckig(0.01);  // control cycle
    InputParameter<3> input;
    OutputParameter<3> output;

    // Set input parameters
    input.current_position = {0.0, 0.0, 0.5};
    input.current_velocity = {0.0, -2.2, -0.5};
    input.current_acceleration = {0.0, 2.5, -0.5};

    input.target_position = {5.0, -2.0, -3.5};
    input.target_velocity = {0.0, -0.5, -2.0};
    input.target_acceleration = {0.0, 0.0, 0.5};

    input.max_velocity = {3.0, 1.0, 3.0};
    input.max_acceleration = {3.0, 2.0, 1.0};
    input.max_jerk = {4.0, 3.0, 2.0};

    // We want to break to a paused state and re-accelerate to the normal trajectory
    // using Ruckig's speed control feature.
    // In this example, we have the four phases: 1. start, 2. brake, 3. accelerate, 4. end.
    std::string phase = "start";
    const double speed_change_duration = 1.0;  // [s]
    double time = 0.0;

    // Generate the trajectory within the control loop
    std::cout << "t | t_traj | position" << std::endl;
    std::cout << std::fixed << std::setprecision(2);

    Result result = Result::Working;
    while (result == Result::Working || result == Result::Paused) {
        result = ruckig.update(input, output);

        // The out.time parameter denotes the time on the trajectory,
        // which is not the same as the time in the control loop as soon as the speed is not 1.0.
        time += ruckig.delta_time;

        std::cout << time << " | " << output.time << " | " << pretty_print(output.new_position) << std::endl;

        if (output.time > 1.8 && phase == "start") {
            phase = "brake";
        }
        if (result == Result::Paused && phase == "brake") {
            phase = "accel";
        }
        if (phase == "accel" && ruckig.speed >= 1.0) {
            phase = "end";
        }

        if (phase == "brake") {
            ruckig.speed = std::max(ruckig.speed - ruckig.delta_time / speed_change_duration, 0.0);
        } else if (phase == "accel") {
            ruckig.speed = std::min(ruckig.speed + ruckig.delta_time / speed_change_duration, 1.0);
        }

        output.pass_to_input(input);
    }

    std::cout << "Trajectory duration: " << output.trajectory.get_duration() << " [s]." << std::endl;
}
