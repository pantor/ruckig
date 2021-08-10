from pathlib import Path
from sys import path
import time

import numpy as np
import toppra as ta
import toppra.constraint as constraint
import toppra.algorithm as algo

from plotter import Plotter

path.insert(0, str(Path(__file__).parent.absolute().parent / 'build'))

from ruckig import InputParameter, Ruckig


class SimpleOut:
    time = None
    new_position = []
    new_velocity = []
    new_acceleration = []


def generate_new_problem(i, seed=9):
    # way_pts = np.random.uniform(-2, 2, size=(4, 3))
    # way_pts = np.array([
    #     np.linspace(0.0, 1.0, 20),
    #     np.linspace(0.01, 1.01, 20),
    #     np.linspace(0.02, 1.02, 20)
    # ]).T
    way_pts = np.array([
        [0.2, 0.1, 0.3],
        [0.3, 0.2, 0.4],
        [0.4, 0.3, 0.5],
        [0.4, 0.4, 0.6],
        [0.5, 0.5, 0.7],
    ])

    way_pts = np.concatenate([[[0, 0, 0]], way_pts, [[1, 1, 1]]])

    return (
        np.linspace(0, 1, way_pts.shape[0]),
        way_pts,
        [2, 2, 2],
        [2, 2, 2],
    )


if __name__ == '__main__':
    np.random.seed(42)
    ta.setup_logging("INFO")

    # durations = 0.0
    # for i in range(250):
    #     ss, way_pts, vlims, alims = generate_new_problem(i)
    #     # print(way_pts)

    #     path = ta.SplineInterpolator(ss, way_pts)
    #     pc_vel = constraint.JointVelocityConstraint(vlims)
    #     pc_acc = constraint.JointAccelerationConstraint(alims)

    #     instance = algo.TOPPRA([pc_vel, pc_acc], path, parametrizer="ParametrizeConstAccel")
    #     s = time.time()
    #     jnt_traj = instance.compute_trajectory()
    #     e = time.time()
    #     durations += (e - s) * 1000
    #     # durations += jnt_traj.duration
    # print(durations/250)

    ss, way_pts, vlims, alims = generate_new_problem(None)
    path = ta.SplineInterpolator(ss, way_pts)
    pc_vel = constraint.JointVelocityConstraint(vlims)
    pc_acc = constraint.JointAccelerationConstraint(alims)

    instance = algo.TOPPRA([pc_vel, pc_acc], path, parametrizer="ParametrizeConstAccel")
    s = time.time()
    jnt_traj = instance.compute_trajectory()


    otg = Ruckig(3, 0.01)
    inp = InputParameter(3)
    inp.max_jerk = [1000, 1000, 1000]
    inp.max_acceleration = alims
    inp.max_velocity = vlims

    out_list = []
    ts_sample = np.linspace(0, jnt_traj.duration, 100)
    for t in ts_sample:
        out = SimpleOut()
        out.time = t
        out.new_position = jnt_traj(t)
        out.new_velocity = jnt_traj(t, 1)
        out.new_acceleration = jnt_traj(t, 2)
        out_list.append(out)

    Plotter.plot_trajectory('otg_toppra.png', otg, inp, out_list, plot_jerk=False)
