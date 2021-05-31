<div align="center">
  <h1 align="center">Ruckig</h1>
  <h3 align="center">
    Online Trajectory Generation. Real-time. Time-optimal. Jerk-constrained.
  </h3>
</div>
<p align="center">
  <a href="https://github.com/pantor/ruckig/actions">
    <img src="https://github.com/pantor/ruckig/workflows/CI/badge.svg" alt="CI">
  </a>

  <a href="https://github.com/pantor/ruckig/issues">
    <img src="https://img.shields.io/github/issues/pantor/ruckig.svg" alt="Issues">
  </a>

  <a href="https://github.com/pantor/ruckig/releases">
    <img src="https://img.shields.io/github/v/release/pantor/ruckig.svg?include_prereleases&sort=semver" alt="Releases">
  </a>

  <a href="https://github.com/pantor/ruckig/blob/master/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green.svg" alt="LGPL">
  </a>
</p>

Ruckig calculates a time-optimal trajectory to a *target* waypoint with position, velocity, and acceleration starting from *any* initial state limited by velocity, acceleration, and jerk constraints. Ruckig is a more powerful and open-source alternative to the [Reflexxes Type IV](http://reflexxes.ws/) library. In fact, Ruckig is the first Type V trajectory generator for arbitrary target states and even supports directional velocity and acceleration limits, while also being faster on top. For robotics and machining applications, Ruckig allows both instant reactions to unforeseen events as well as simple offline trajectory planning.

More information can be found in the corresponding paper [Jerk-limited Real-time Trajectory Generation with Arbitrary Target States](https://arxiv.org/abs/2105.04830), accepted for the *Robotics: Science and Systems (RSS), 2021* conference.


## Installation

Ruckig has no dependencies (except for testing). To build Ruckig using CMake, just run

```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

To install Ruckig in a system-wide directory, use `(sudo) make install`. An example of using Ruckig in your CMake project is given by `examples/CMakeLists.txt`. However, you can also include Ruckig as a directory within your project and call `add_subdirectory(ruckig)` in your parent `CMakeLists.txt`.

Ruckig is also available as a Python module, in particular for development or debugging purposes. It can be installed from [PyPI](https://pypi.org/project/ruckig/) via
```bash
pip install ruckig
```
When using CMake, the Python module can be built using the `BUILD_PYTHON_MODULE` flag. If you're only interested in the Python module (and not in the C++ library), you can build and install Ruckig via `pip install .`.


## Tutorial

Furthermore, a tutorial will explain the basics to include online generated trajectories within your application. A working example can be found in the `examples` directory. A time-optimal trajectory for a single degree of freedom is shown in the figure below.

![Trajectory Profile](https://github.com/pantor/ruckig/raw/master/doc/example_profile.png?raw=true)


### Waypoint-based Trajectory Generation

Ruckig provides three main interface classes: the *Ruckig*, the *InputParameter*, and the *OutputParameter* class.

First, you'll need to create a Ruckig instance with the number of DoFs as a template parameter, and the control cycle (e.g. in seconds) in the constructor.

```.cpp
Ruckig<6> ruckig {0.001}; // Number DoFs; control cycle in [s]
```

The input type has 3 blocks of data: the *current* state, the *target* state and the corresponding kinematic *limits*.

```.cpp
InputParameter<6> input; // Number DoFs
input.current_position = {0.2, ...};
input.current_velocity = {0.1, ...};
input.current_acceleration = {0.1, ...};
input.target_position = {0.5, ...};
input.target_velocity = {-0.1, ...};
input.target_acceleration = {0.2, ...};
input.max_velocity = {0.4, ...};
input.max_acceleration = {1.0, ...};
input.max_jerk = {4.0, ...};

OutputParameter<6> output; // Number DoFs
```

Given all input and output resources, we can iterate over the trajectory at each discrete time step. For most applications, this loop must run within a real-time thread and controls the actual hardware.

```.cpp
while (ruckig.update(input, output) == Result::Working) {
  // Make use of the new state here!

  input.current_position = output.new_position;
  input.current_velocity = output.new_velocity;
  input.current_acceleration = output.new_acceleration;
}
```

During your update step, you'll need to copy the new kinematic state into the current state. If the current state is not the expected, pre-calculated trajectory, ruckig will calculate a new trajectory with the novel input. When the trajectory has reached the target state, the `update` function will return `Result::Finished`.


### Input Parameter

To go into more detail, the *InputParameter* type has following members:

```.cpp
using Vector = std::array<double, DOFs>;

Vector current_position;
Vector current_velocity; // Initialized to zero
Vector current_acceleration; // Initialized to zero

Vector target_position;
Vector target_velocity; // Initialized to zero
Vector target_acceleration; // Initialized to zero

Vector max_velocity;
Vector max_acceleration;
Vector max_jerk;

std::optional<Vector> min_velocity; // If not given, the negative maximum velocity will be used.
std::optional<Vector> min_acceleration; // If not given, the negative maximum acceleration will be used.

std::array<bool, DOFs> enabled; // Initialized to true
std::optional<double> minimum_duration;

Interface interface; // The default position interface controls the full kinematic state.
Synchronization synchronization; // Synchronization behavior of multiple DoFs
DurationDiscretization duration_discretization; // Whether the duration should be a discrete multiple of the control cycle (off by default)
```

Members are implemented using the C++ standard `array` and `optional` type. Note that there are range constraints due to numerical reasons, see below for more details. To check the input before a calculation step, the `ruckig.validate_input(input)` method returns `false` if an input is not valid. Of course, the target state needs to be within the given kinematic limits. Additionally, the target acceleration needs to fulfil
```
Abs(target_acceleration) <= Sqrt(2 * max_jerk * (max_velocity - Abs(target_velocity)))
```
If a DoF is not enabled, it will be ignored in the calculation. A minimum duration can be optionally given. Furthermore, the minimum velocity and acceleration can be specified. If it is not given, the negative maximum velocity or acceleration will be used (similar to the jerk limit). For example, this might be useful in human robot collaboration settings with a different velocity limit towards a human. Or, the dynamic limits at a given configuration of the robot can be approximated much better with different acceleration limits.

Furthermore, there are some options for advanced functionality, e.g.:
- for different synchronization behavior (i.a. phase or time synchonization).
- for the control interface (position or velocity control).
- for discrete trajectory durations.

We refer to the [API documentation](https://pantor.github.io/ruckig/namespaceruckig.html) of the enumerations within the `ruckig` namespace for all available options.


### Result Type

The `update` function of the Ruckig class returns a Result type that indicates the current state of the algorithm. Currently, this can either be **working**, **finished** if the trajectory has finished, or an **error** type if something went wrong during calculation. The result type can be compared as a standard integer.

State                           | Error Code
------------------------------- | ----------
Working                         | 0
Finished                        | 1
Error                           | -1
ErrorInvalidInput               | -100
ErrorTrajectoryDuration         | -101
ErrorExecutionTimeCalculation   | -110
ErrorSynchronizationCalculation | -111


### Output Parameter

The output class gives the new kinematic state of the trajectory.

```.cpp
Vector new_position;
Vector new_velocity;
Vector new_acceleration;

bool new_calculation; // Whether a new calculation was performed in the last cycle
double calculation_duration; // Duration of the calculation in the last cycle [µs]

Trajectory trajectory; // The current trajectory
double time; // The current, auto-incremented time. Reset to 0 at a new calculation.
```
Moreover, the trajectory class has a range of useful parameters and methods.

```.cpp
double duration; // Duration of the trajectory
std::array<double, DOFs> independent_min_durations; // Time-optimal profile for each independent DoF

<...> at_time(double time); // Get the kinematic state of the trajectory at a given time
<...> get_position_extrema(); // Returns information about the position extrema and their times
```
Again, we refer to the [API documentation](https://pantor.github.io/ruckig/) for the exact signatures.


## Tests and Numerical Stability

The current test suite validates over 5.000.000.000 random trajectories. The numerical exactness is tested for the final position and final velocity to be within `1e-8`, for the velocity, acceleration and jerk limit to be within `1e-12`, and for the final acceleration as well to be within a numerical error of `1e-12`. The maximal supported trajectory duration is `7e3`, which sounds short but should suffice for most applications seeking for time-optimality. Note that Ruckig will also output values outside of this range, there is however no guarantee for correctness.


## Benchmark

We find that Ruckig is around twice as fast as Reflexxes Type IV and well-suited for control cycles as low as half a millisecond.

![Benchmark](https://github.com/pantor/ruckig/raw/master/doc/benchmark.png?raw=true)


## Development

Ruckig is written in C++17. It is continuously tested on `ubuntu-latest`, `macos-latest`, and `windows-latest` against following versions

- Doctest v2.4 (only for testing)
- Pybind11 v2.6 (only for python wrapper)


## Citation

```
@article{berscheid2021jerk,
  title={Jerk-limited Real-time Trajectory Generation with Arbitrary Target States},
  author={Berscheid, Lars and Kr{\"o}ger, Torsten},
  journal={Robotics: Science and Systems XVII},
  year={2021}
}
```
