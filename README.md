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

Ruckig calculates a time-optimal trajectory given a *target* waypoint with position, velocity, and acceleration, starting from *any* initial state limited by velocity, acceleration, and jerk constraints. Robotics. Machine control. Ruckig is a more powerful and open-source alternative to the [Reflexxes Type IV](http://reflexxes.ws/) library. In fact, Ruckig is a Type V trajectory generator. In general, Ruckig allows for instant reactions to unforeseen events.


## Installation

For normal usage, Ruckig has no dependencies. To build Ruckig using CMake, just 

```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

To install Ruckig in a system-wide directory, use `(sudo) make install`. We recommend to include Ruckig as a directory within your project and call `add_subdirectory(ruckig)` in your main `CMakeLists.txt`. A python module can be built using the `BUILD_PYTHON_MODULE` CMake flag.


## Tutorial

Currently only the (more-complex) *position* interface is implemented. A time-optimal trajectory is shown in the figure below.

![Trajectory Profile](/doc/example_profile.png?raw=true)

Furthermore, a tutorial will explain the basics to include online generated trajectories within your robotics or machining application. A working example can be found in the `examples` directory.


### Real-time trajectory generation

Ruckig provides three interface classes: the *Ruckig*, the *InputParameter*, and the *OutputParameter* class. 

First, you'll need to create a Ruckig instance with the number of DoFs as a template parameter, and the control cycle in seconds in the constructor.

```c++
Ruckig<6> ruckig {0.001}; // Number DoFs; control cycle in [s]
```

The input type has 3 blocks of data: the *current* state, the *target* state and the corresponding dynamical *limits*.

```c++
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

while (otg.update(input, output) == Result::Working) {
  // Make use of the new dynamic state here!

  input.current_position = output.new_position;
  input.current_velocity = output.new_velocity;
  input.current_acceleration = output.new_acceleration;
}
```

During your update step, you'll need to copy the new dynamic state into the current state. If the current state is not the expected, pre-calculated trajectory, ruckig will calculate a new trajectory with the new input. The output state at a given time using the last calculated trajectory is returned by the `ruckig.at_time(t)` function. 


### Input Parameter

The *InputParameter* type has following members: 

```c++
std::array<double, DOFs> current_position;
std::array<double, DOFs> current_velocity; // Initialized to zero
std::array<double, DOFs> current_acceleration; // Initialized to zero

std::array<double, DOFs> target_position;
std::array<double, DOFs> target_velocity; // Initialized to zero
std::array<double, DOFs> target_acceleration; // Initialized to zero

std::array<double, DOFs> max_velocity;
std::array<double, DOFs> max_acceleration;
std::array<double, DOFs> max_jerk;

std::array<bool, DOFs> enabled; // Initialized to true
std::optional<double> minimum_duration;
```

To check the input in front of a calculation step, the `ruckig.validate_input(input)` method returns `false` if an input is not valid. Of course, the target state needs to be within the given dynamic limits. Additionally, the target acceleration needs to fulfill
```
target_acceleration <= Sqrt(2 * max_jerk * (max_velocity - Abs(target_velocity)))
``` 
If a DoF is not enabled, it will be ignored in the calculation. A minimum duration can be optionally given.


### Result Type

The `update` function of the Ruckig class returns a Result type that indicates the current state of the algorithm. Currently, this can either be **working**, **finished** if the trajectory has finished, or **error** if something went wrong during calculation. In this case, an exception (see below for more details) is thrown.

State                           | Error Code
------------------------------- | ----------
Working                         | 0
Finished                        | 1
Error                           | -1
ErrorInvalidInput               | -100
ErrorExecutionTimeCalculation   | -110
ErrorSynchronizationCalculation | -111


### Output Parameter

The output class gives the new dynamical state of the trajectory.

```c++
std::array<double, DOFs> new_position;
std::array<double, DOFs> new_velocity;
std::array<double, DOFs> new_acceleration;

double duration; // Duration of the trajectory [s]
bool new_calculation; // Whether a new calactuion was performed in the last cycle
double calculation_duration; // Duration of the calculation in the last cycle [µs]

std::array<double, DOFs> independent_min_durations; // [s]
```
Moreover, a range of additional parameter about the duration of the trajectory are included.


## Tests

The current test suite validates over 918.000 (random) trajectories. The numerical exactness is tested for the position, velocity, acceleration, and time target to be within `1e-8`, for the velocity and acceleration limit to be withing `1e-9`, and for the jerk limit to be within a numerical error of `1e-12`.


## Development

Ruckig is written in C++17. It is currently tested against following versions

- Eigen v3.3.9 (only for testing)
- Catch2 v2.13.3 (only for testing)
- Reflexxes v1.2.7 (only for testing)
- Pybind11 v2.6.0 (only for prototyping)


## Citation

A paper will follow...
