import time
from pathlib import Path
from sys import path

# Path to the build directory including a file similar to 'ruckig.cpython-37m-x86_64-linux-gnu'.
build_path = Path(__file__).parent.absolute().parent / 'build'
path.insert(0, str(build_path))

from ruckig import InputParameter, Ruckig, Trajectory, Result


if __name__ == '__main__':
    inp = InputParameter(3)
    inp.current_position = [0.2, 0, -1]
    inp.current_velocity = [0, 0.2, 0]
    inp.current_acceleration = [0, 1, 0]
    inp.target_position = [0, -1, -1]
    inp.target_velocity = [0.2, 0, 0]
    inp.target_acceleration = [0, 0.1, -0.1]
    inp.max_velocity = [2, 1, 1]
    inp.max_acceleration = [0.2, 2, 2]
    inp.max_jerk = [3, 4, 5]

    # We don't need to pass the control rate (cycle time) when using only offline features
    otg = Ruckig(3)
    trajectory = Trajectory(3)

    # We now calculate the trajectory in an offline manner
    result = otg.calculate(inp, trajectory)
    if result == Result.ErrorInvalidInput:
        raise Exception('Invalid input!')

    print(f'Trajectory duration: {trajectory.duration:0.4f} [s]')

    new_time = 1.0

    # Then, we can calculate the kinematic state at a given time
    new_position, new_velocity, new_acceleration = trajectory.at_time(new_time)

    print(f'Position at time {new_time:0.4f} [s]: {new_position}')
    print(f'Position extremas are {trajectory.position_extrema}')
