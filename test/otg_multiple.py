from copy import copy
from pathlib import Path
from sys import path

from plotter import Plotter

path.insert(0, str(Path(__file__).parent.parent / 'build'))

from ruckig import InputParameter, OutputParameter, Result, Ruckig, Synchronization, Interface, DurationDiscretization
from ruckig import Reflexxes


def walk_through_trajectory(otg, inp, intermediate_targets):
    t_list, out_list = [], []
    out = OutputParameter(inp.degrees_of_freedom)

    old_target = inp.target_position, inp.target_velocity, inp.target_acceleration
    intermediate_targets.append(old_target)

    time_offset = 0.0
    for inp.target_position, inp.target_velocity, inp.target_acceleration in intermediate_targets:
        res = Result.Working
        while res == Result.Working:
            res = otg.update(inp, out)

            inp.current_position = out.new_position
            inp.current_velocity = out.new_velocity
            inp.current_acceleration = out.new_acceleration

            t_list.append(time_offset + out.time)
            out_list.append(copy(out))

        time_offset += out.trajectory.duration

    return t_list, out_list, time_offset


if __name__ == '__main__':
    inp = InputParameter(3)
    inp.current_position = [0, 0, 0]
    inp.current_velocity = [0, 0, 0]
    inp.current_acceleration = [0, 0, 0]
    inp.target_position = [1, 0, 0]
    inp.target_velocity = [0, 0, 0]
    inp.target_acceleration = [0, 0, 0]
    inp.max_velocity = [100, 100, 100]
    inp.max_acceleration = [100, 100, 100]
    inp.max_jerk = [1, 1, 1]

    intermediate_targets = [(
        [0.7, 0.0, 0.0],
        [0.56, 0.0, 0.0],
        [-0.3, 0, 0],
    )]

    otg = Ruckig(inp.degrees_of_freedom, 0.005)

    t_list, out_list, duration = walk_through_trajectory(otg, inp, intermediate_targets)


    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {duration:0.4f} [s]')

    Plotter.plot_trajectory('otg_trajectory.png', otg, inp, t_list, out_list)
