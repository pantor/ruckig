import copy
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / 'build'))

from _ruckig import Quintic, InputParameter, OutputParameter, Result, Ruckig, Smoothie
from _ruckig import Reflexxes


def walk_through_trajectory(otg, inp):
    t = 0.0
    t_list, out_list = [], []
    out = OutputParameter()

    res = Result.Working
    while res == Result.Working:
        res = otg.update(inp, out)

        inp.current_position = out.new_position
        inp.current_velocity = out.new_velocity
        inp.current_acceleration = out.new_acceleration

        t_list.append(t)
        out_list.append(copy.copy(out))
        t += otg.delta_time

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
        plt.plot(t_list, qaxis[:, dof], label=f'r_{dof+1}')
        plt.plot(t_list, dqaxis[:, dof], label=f'v_{dof+1}')
        plt.plot(t_list, ddqaxis[:, dof], label=f'a_{dof+1}')
        plt.plot(t_list, dddqaxis[:, dof], label=f'j_{dof+1}')

        # Plot limit lines
        if inp.max_velocity[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_velocity[dof], color='orange', linestyle='--', linewidth=1.1)

        if -inp.max_velocity[dof] > 1.4 * global_min:
            plt.axhline(y=-inp.max_velocity[dof], color='orange', linestyle='--', linewidth=1.1)

        if inp.max_acceleration[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_acceleration[dof], color='g', linestyle='--', linewidth=1.1)

        if -inp.max_acceleration[dof] > 1.4 * global_min:
            plt.axhline(y=-inp.max_acceleration[dof], color='g', linestyle='--', linewidth=1.1)

        if inp.max_jerk[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_jerk[dof], color='red', linestyle='--', linewidth=1.1)

        if -inp.max_jerk[dof] > 1.4 * global_min:
            plt.axhline(y=-inp.max_jerk[dof], color='red', linestyle='--', linewidth=1.1)

        plt.legend()
        plt.grid(True)

    plt.xlabel('t')
    plt.savefig(Path(__file__).parent.parent / 'build' / 'otg_trajectory.png')
    # plt.show()

def print_input_for_mathematica(inp, dof, tf=None):
    result = f"""p0->{inp.current_position[dof]}, \
v0->{inp.current_velocity[dof]}, \
a0->{inp.current_acceleration[dof]}, \
pf->{inp.target_position[dof]}, \
vf->{inp.target_velocity[dof]}, \
af->{inp.target_acceleration[dof]}, \
vMax->{inp.max_velocity[dof]}, \
aMax->{inp.max_acceleration[dof]}, \
jMax->{inp.max_jerk[dof]}"""

    if tf:
        result += f', tf->{tf}'
    print('{ ' + result + ' }')


if __name__ == '__main__':
    inp = InputParameter()
    inp.current_position = [0.08031159597, -0.2030066108, 0.06789170442]
    inp.current_velocity = [-0.9441238828,   0.109902212, -0.8735232222]
    inp.current_acceleration = [ 0.6952040874,  0.2950967826, -0.3083741829]
    inp.target_position = [-0.8448925111,  -0.108433588, -0.4433133944]
    inp.target_velocity = [-0.7682199151,  0.5278867341,  0.1923404258]
    inp.target_acceleration = [0, 0, 0]
    inp.max_velocity = [7.523579228,  4.00390983,  1.11251056]
    inp.max_acceleration = [4.699438046, 2.855245582, 7.512490303]
    inp.max_jerk = [ 1.82452848, 4.250166602,  8.24992113]
    inp.minimum_duration = None

    print_input_for_mathematica(inp, 1)

    # otg = Quintic(0.005)
    # otg = Smoothie(0.005)
    # otg = Reflexxes(0.005)
    otg = Ruckig(0.005)

    t_list, out_list = walk_through_trajectory(otg, inp)

    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {out_list[0].duration:0.3f} [s]')

    plot_trajectory(t_list, out_list)
