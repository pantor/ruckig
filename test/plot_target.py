from copy import copy
from pathlib import Path
from sys import path

from plotter import Plotter

path.insert(0, str(Path(__file__).parent.absolute().parent / 'build'))

from ruckig import InputParameter, OutputParameter, Result, Ruckig, Synchronization, ControlInterface, DurationDiscretization


def walk_through_trajectory(otg, inp):
    out_list = []
    out = OutputParameter(inp.degrees_of_freedom)

    res = Result.Working
    while res == Result.Working:
        res = otg.update(inp, out)

        out.pass_to_input(inp)
        out_list.append(copy(out))

    return out_list


if __name__ == '__main__':
    inp = InputParameter(3)
    # inp.control_interface = ControlInterface.Velocity
    # inp.synchronization = Synchronization.No
    # inp.duration_discretization = DurationDiscretization.Discrete

    # inp.per_dof_control_interface = [ControlInterface.Position, ControlInterface.Velocity, ControlInterface.Position]
    # inp.per_dof_synchronization = [Synchronization.Phase, Synchronization.Time, Synchronization.Phase]

    inp.current_position = [0.0, -2.0, 0.0]
    inp.current_velocity = [0.0, 0.0, 0.0]
    inp.current_acceleration = [0.0, 0.0, 0.0]

    inp.target_position = [1.0, -3.0, 2.0]
    inp.target_velocity = [0.0, 0.0, 0.0]
    inp.target_acceleration = [0.0, 0.0, 0.0]

    inp.max_velocity = [1.0, 1.0, 1.0]
    inp.max_acceleration = [1.0, 1.0, 1.0]
    inp.max_jerk = [1.0, 1.0, 1.0]

    # inp.minimum_duration = 5.0


    otg = Ruckig(inp.degrees_of_freedom, 0.001)

    out_list = walk_through_trajectory(otg, inp)

    # print(out_list[0].trajectory.position_extrema)
    # print(out_list[0].trajectory.independent_min_durations)

    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {out_list[0].trajectory.duration:0.4f} [s]')

    Plotter.plot_trajectory('otg_trajectory.png', otg, inp, out_list, plot_jerk=False)
