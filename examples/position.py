from copy import copy
from pathlib import Path
from sys import path

# Path to the build directory including a file similar to 'ruckig.cpython-37m-x86_64-linux-gnu'.
build_path = Path(__file__).parent.parent / 'build'
path.insert(0, str(build_path / 'Debug'))  # For windows
path.insert(0, str(build_path))

from ruckig import InputParameter, OutputParameter, Result, Ruckig


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

    otg = Ruckig(3, 0.05)

    out = OutputParameter(3)
    first_output = None

    print('\t'.join(['t'] + [str(i) for i in range(otg.degrees_of_freedom)]))

    res = Result.Working
    while res == Result.Working:
        res = otg.update(inp, out)

        print('\t'.join([f'{out.time:0.3f}'] + [f'{p:0.3f}' for p in out.new_position]))

        inp.current_position = out.new_position
        inp.current_velocity = out.new_velocity
        inp.current_acceleration = out.new_acceleration

        if not first_output:
            first_output = copy(out)

    print(f'Calculation duration: {first_output.calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {first_output.trajectory.duration:0.4f} [s]')
