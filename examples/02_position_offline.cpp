#include <ruckig/ruckig.hpp>

#include "plotter.hpp"


using namespace ruckig;

int main() {
    const size_t degrees_of_freedom = 3;

    // Create input parameters
    InputParameter input(degrees_of_freedom);
    input.current_position = {0.0, 0.0, 0.5};
    input.current_velocity = {0.0, -2.2, -0.5};
    input.current_acceleration = {0.0, 2.5, -0.5};

    input.target_position = {5.0, -2.0, -3.5};
    input.target_velocity = {0.0, -0.5, -2.0};
    input.target_acceleration = {0.0, 0.0, 0.5};

    input.max_velocity = {3.0, 1.0, 3.0};
    input.max_acceleration = {3.0, 2.0, 1.0};
    input.max_jerk = {4.0, 3.0, 2.0};

    // Set different constraints for negative direction
    input.min_velocity = {-2.0, -0.5, -3.0};
    input.min_acceleration = {-2.0, -2.0, -2.0};

    // We don't need to pass the control rate (cycle time) when using only offline features
    Ruckig ruckig(degrees_of_freedom);
    Trajectory trajectory(degrees_of_freedom);

    // Calculate the trajectory in an offline manner (outside of the control loop)
    Result result = ruckig.calculate(input, trajectory);
    if (result == Result::ErrorInvalidInput) {
        std::cout << "Invalid input!" << std::endl;
        return -1;
    }

    // Get duration of the trajectory
    std::cout << "Trajectory duration: " << trajectory.get_duration() << " [s]." << std::endl;

    double new_time = 1.0;

    // Then, we can calculate the kinematic state at a given time
    std::vector<double> new_position(degrees_of_freedom), new_velocity(degrees_of_freedom), new_acceleration(degrees_of_freedom);
    trajectory.at_time(new_time, new_position, new_velocity, new_acceleration);

    std::cout << "Position at time " << new_time << " [s]: " << pretty_print(new_position) << std::endl;

    // Get some info about the position extrema of the trajectory
    std::vector<Bound> position_extrema(degrees_of_freedom);
    trajectory.get_position_extrema(position_extrema);
    std::cout << "Position extremas for DoF 4 are " << position_extrema[2].min << " (min) to " << position_extrema[2].max << " (max)" << std::endl;
}
