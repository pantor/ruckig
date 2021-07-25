from copy import copy
from pathlib import Path
from sys import path

from plotter import Plotter

path.insert(0, str(Path(__file__).parent.parent / 'build'))

from ruckig import InputParameter, OutputParameter, Result, Ruckig, Synchronization, Interface, DurationDiscretization
from ruckig import Reflexxes


def walk_through_trajectory(otg, inp):
    t_list, out_list = [], []
    out = OutputParameter(inp.degrees_of_freedom)

    res = Result.Working
    while res == Result.Working:
        res = otg.update(inp, out)

        inp.current_position = out.new_position
        inp.current_velocity = out.new_velocity
        inp.current_acceleration = out.new_acceleration

        t_list.append(out.time)
        out_list.append(copy(out))

    return t_list, out_list


if __name__ == '__main__':
    inp = InputParameter(3)
    # inp.interface = Interface.Velocity
    inp.synchronization = Synchronization.Phase
    # inp.duration_discretization = DurationDiscretization.Discrete

    inp.current_position = [0, 0, 0]
    inp.current_velocity = [0.2, -0.4, 0.6]
    inp.current_acceleration = [0, 0, 0]
    inp.target_position = [1, -2, 3]
    inp.target_velocity = [0, 0, 0]
    inp.target_acceleration = [0, 0, 0]
    inp.max_velocity = [2, 2, 2]
    inp.max_acceleration = [1, 1, 1]
    inp.max_jerk = [2, 2, 2]

    # otg = Reflexxes(inp.degrees_of_freedom, 0.005)
    otg = Ruckig(inp.degrees_of_freedom, 0.005)

    t_list, out_list = walk_through_trajectory(otg, inp)

    # print(out_list[0].trajectory.position_extrema)

    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {out_list[0].trajectory.duration:0.4f} [s]')

    Plotter.plot_trajectory('otg_trajectory.png', otg, inp, t_list, out_list)
