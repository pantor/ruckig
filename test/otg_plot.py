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
    # inp.synchronization = Synchronization.Phase
    # inp.duration_discretization = DurationDiscretization.Discrete

    inp.current_position = [6.22856405025901, -3.690144255821358, 1.664557711692404]
    inp.current_velocity = [-0.434815099852722, 0.5760942342059757, 0.459449163887151]
    inp.current_acceleration = [0, -2.072115307337411, 0.689486369647744]
    inp.target_position = [0.6959453258241575, 0.3766650922417908, -2.986784389954874]
    inp.target_velocity = [0, 0.05033093271396352, -0.1453937359671692]
    inp.target_acceleration = [1.305000072573007, 0, -0.09904524070763976]
    inp.max_velocity = [0.434815099852722, 1.864268119408352, 2.002840733466842]
    inp.max_acceleration = [9.112772043597223, 4.574478625381595, 10.10043773003713]
    inp.max_jerk = [7.654313459653499, 8.328157077423564, 2.012501506433343]

    # otg = Reflexxes(3, 0.005)
    otg = Ruckig(inp.degrees_of_freedom, 0.005)

    t_list, out_list = walk_through_trajectory(otg, inp)

    # print(out_list[0].trajectory.position_extrema)

    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {out_list[0].trajectory.duration:0.4f} [s]')

    Plotter.plot_trajectory('otg_trajectory.png', otg, inp, t_list, out_list)
