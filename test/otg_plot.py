from copy import copy
from pathlib import Path
from sys import path

import matplotlib.pyplot as plt
import numpy as np

path.insert(0, str(Path(__file__).parent.parent / 'build'))

from ruckig import InputParameter, OutputParameter, Result, Ruckig, Synchronization, Interface, DurationDiscretization
from ruckig import Reflexxes


def walk_through_trajectory(otg, inp, print_table=True):
    t_list, out_list = [], []
    out = OutputParameter(3)

    res = Result.Working
    old_acc = 0
    print_dof = 0
    while res == Result.Working:
        res = otg.update(inp, out)

        inp.current_position = out.new_position
        inp.current_velocity = out.new_velocity
        inp.current_acceleration = out.new_acceleration

        if print_table:
            jerk = (old_acc - out.new_acceleration[print_dof]) / otg.delta_time
            old_acc = out.new_acceleration[print_dof]
            # print(str(out.time) + '\t' + str(inp.current_position[print_dof]) + '\t' + str(inp.current_velocity[print_dof]) + '\t' + str(inp.current_acceleration[print_dof]) + '\t' + str(jerk))
            # print(str(inp.current_position[0]) + '\t' + str(inp.current_position[1]))

        t_list.append(out.time)
        out_list.append(copy(out))

    return t_list, out_list


def plot_trajectory(t_list, out_list):
    qaxis = np.array(list(map(lambda x: x.new_position, out_list)))
    dqaxis = np.array(list(map(lambda x: x.new_velocity, out_list)))
    ddqaxis = np.array(list(map(lambda x: x.new_acceleration, out_list)))
    dddqaxis = np.diff(ddqaxis, axis=0, prepend=ddqaxis[0, 0]) / otg.delta_time
    dddqaxis[0, :] = 0.0
    dddqaxis[-1, :] = 0.0

    plt.figure(figsize=(8.0, 2.0 + 3.0 * inp.degrees_of_freedom), dpi=120)

    for dof in range(inp.degrees_of_freedom):
        global_max = np.max([qaxis[:, dof], dqaxis[:, dof], ddqaxis[:, dof], dddqaxis[:, dof]])
        global_min = np.min([qaxis[:, dof], dqaxis[:, dof], ddqaxis[:, dof], dddqaxis[:, dof]])

        plt.subplot(inp.degrees_of_freedom, 1, dof + 1)
        plt.plot(t_list, qaxis[:, dof], label=f'Position {dof+1}')
        plt.plot(t_list, dqaxis[:, dof], label=f'Velocity {dof+1}')
        plt.plot(t_list, ddqaxis[:, dof], label=f'Acceleration {dof+1}')
        plt.plot(t_list, dddqaxis[:, dof], label=f'Jerk {dof+1}')

        # Plot limit lines
        if inp.max_velocity[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_velocity[dof], color='orange', linestyle='--', linewidth=1.1)

        min_velocity = inp.min_velocity[dof] if inp.min_velocity else -inp.max_velocity[dof]
        if min_velocity > 1.4 * global_min:
            plt.axhline(y=min_velocity, color='orange', linestyle='--', linewidth=1.1)

        if inp.max_acceleration[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_acceleration[dof], color='g', linestyle='--', linewidth=1.1)

        min_acceleration = inp.min_acceleration[dof] if inp.min_acceleration else -inp.max_acceleration[dof]
        if min_acceleration > 1.4 * global_min:
            plt.axhline(y=min_acceleration, color='g', linestyle='--', linewidth=1.1)

        if inp.max_jerk[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_jerk[dof], color='red', linestyle='--', linewidth=1.1)

        if -inp.max_jerk[dof] > 1.4 * global_min:
            plt.axhline(y=-inp.max_jerk[dof], color='red', linestyle='--', linewidth=1.1)

        plt.legend()
        plt.grid(True)

    plt.xlabel('t')
    plt.savefig(Path(__file__).parent.parent / 'build' / 'otg_trajectory.png')
    # plt.show()


if __name__ == '__main__':
    inp = InputParameter(3)
    # inp.interface = Interface.Velocity
    # inp.synchronization = Synchronization.Phase
    # inp.duration_discretization = DurationDiscretization.Discrete


    for i in [-299]:
    # for i in range(200):
        eps = (i - 100) * 1e-16
        print(i)

        inp.current_position = [6.22856405025901, -3.690144255821358, 1.664557711692404]
        inp.current_velocity = [-0.434815099852722, 0.5760942342059757, 0.459449163887151]
        inp.current_acceleration = [0, -2.072115307337411, 0.689486369647744]
        inp.target_position = [0.6959453258241575, 0.3766650922417908, -2.986784389954874]
        inp.target_velocity = [0, 0.05033093271396352, -0.1453937359671692]
        inp.target_acceleration = [1.305000072573007, 0, -0.09904524070763976]
        inp.max_velocity = [0.434815099852722, 1.864268119408352, 2.002840733466842]
        inp.max_acceleration = [9.112772043597223, 4.574478625381595, 10.10043773003713]
        inp.max_jerk = [7.654313459653499, 8.328157077423564  +eps, 2.012501506433343]

        # otg = Reflexxes(3, 0.005)
        otg = Ruckig(3, 0.005)

        t_list, out_list = walk_through_trajectory(otg, inp)

    # print(out_list[0].trajectory.get_position_extrema())

    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {out_list[0].trajectory.duration:0.4f} [s]')

    plot_trajectory(t_list, out_list)
