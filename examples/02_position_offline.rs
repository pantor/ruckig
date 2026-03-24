use ruckig::{InputParameter, Result, Ruckig, Trajectory};

fn main() {
    // Create input parameters
    let mut input = InputParameter::new(3);
    input.current_position = vec![0.0, 0.0, 0.5];
    input.current_velocity = vec![0.0, -2.2, -0.5];
    input.current_acceleration = vec![0.0, 2.5, -0.5];

    input.target_position = vec![5.0, -2.0, -3.5];
    input.target_velocity = vec![0.0, -0.5, -2.0];
    input.target_acceleration = vec![0.0, 0.0, 0.5];

    input.max_velocity = vec![3.0, 1.0, 3.0];
    input.max_acceleration = vec![3.0, 2.0, 1.0];
    input.max_jerk = vec![4.0, 3.0, 2.0];

    // We don't need to pass the control rate (cycle time) when using only offline features
    let mut ruckig = Ruckig::new_direct_and_offline(3);
    let mut trajectory = Trajectory::new_direct(3);

    // Calculate the trajectory in an offline manner (outside of the control loop)
    let result = ruckig.calculate(&mut input, &mut trajectory);
    if result == Result::ErrorInvalidInput {
        println!("Invalid input!");
        return;
    }

    // Get duration of the trajectory
    println!("Trajectory duration: {} [s].", trajectory.get_duration());
}

