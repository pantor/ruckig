#[cxx::bridge]
mod ffi {
    unsafe extern "C++" {
        include!("rust.hpp");

        type Trajectory;
        type InputParameter;
        type OutputParameter;
        type Ruckig;

        // Trajectory
        fn traj_new_direct(dofs: usize) -> UniquePtr<Trajectory>;
        fn traj_new(dofs: usize, waypoints: usize) -> UniquePtr<Trajectory>;
        fn traj_get_dofs(traj: &Trajectory) -> usize;
        fn traj_get_duration(traj: &Trajectory) -> f64;
        unsafe fn traj_at_time(traj: &Trajectory, time: f64, pos: *mut f64, vel: *mut f64, acc: *mut f64, len: usize);

        // InputParameter
        fn input_new(dofs: usize) -> UniquePtr<InputParameter>;
        fn input_get_control_interface(input: &InputParameter) -> i32;
        fn input_set_control_interface(input: Pin<&mut InputParameter>, v: i32);
        fn input_get_synchronization(input: &InputParameter) -> i32;
        fn input_set_synchronization(input: Pin<&mut InputParameter>, v: i32);
        unsafe fn input_set_current_position(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_current_position(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_current_velocity(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_current_velocity(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_current_acceleration(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_current_acceleration(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_target_position(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_target_position(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_target_velocity(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_target_velocity(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_target_acceleration(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_target_acceleration(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_max_velocity(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_max_velocity(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_max_acceleration(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_max_acceleration(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_max_jerk(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_max_jerk(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_max_position(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_max_position(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_min_position(input: Pin<&mut InputParameter>, data: *const f64, len: usize);
        fn input_get_min_position(input: &InputParameter) -> &CxxVector<f64>;
        unsafe fn input_set_intermediate_positions(input: Pin<&mut InputParameter>, data: *const f64, dofs: usize, count: usize);
        fn input_clear_intermediate_positions(input: Pin<&mut InputParameter>);
        unsafe fn input_get_intermediate_positions(input: &InputParameter, buf: *mut f64, capacity: usize) -> usize;

        // OutputParameter
        fn output_new_direct(dofs: usize) -> UniquePtr<OutputParameter>;
        fn output_new(dofs: usize, waypoints: usize) -> UniquePtr<OutputParameter>;
        fn output_get_trajectory(output: &OutputParameter) -> &Trajectory;
        fn output_get_new_position(output: &OutputParameter) -> &CxxVector<f64>;
        fn output_get_new_velocity(output: &OutputParameter) -> &CxxVector<f64>;
        fn output_get_new_acceleration(output: &OutputParameter) -> &CxxVector<f64>;
        fn output_get_time(output: &OutputParameter) -> f64;
        fn output_get_new_section(output: &OutputParameter) -> usize;
        fn output_get_new_calculation(output: &OutputParameter) -> bool;
        fn output_get_calculation_duration(output: &OutputParameter) -> f64;
        fn output_pass_to_input(output: &OutputParameter, input: Pin<&mut InputParameter>);

        // Ruckig
        fn ruckig_new_direct_and_offline(dofs: usize) -> UniquePtr<Ruckig>;
        fn ruckig_new_direct(dofs: usize, delta_time: f64) -> UniquePtr<Ruckig>;
        fn ruckig_new_offline(dofs: usize, waypoints: usize) -> UniquePtr<Ruckig>;
        fn ruckig_new(dofs: usize, delta_time: f64, waypoints: usize) -> UniquePtr<Ruckig>;
        fn ruckig_update(ruckig: Pin<&mut Ruckig>, input: &InputParameter, output: Pin<&mut OutputParameter>) -> i32;
        fn ruckig_calculate(ruckig: Pin<&mut Ruckig>, input: &InputParameter, trajectory: Pin<&mut Trajectory>) -> i32;
        fn ruckig_reset(ruckig: Pin<&mut Ruckig>);
    }
}

use cxx::UniquePtr;

pub struct Trajectory(UniquePtr<ffi::Trajectory>);

pub struct TrajectoryRef<'a> {
    inner: &'a ffi::Trajectory,
}

fn traj_at_time_impl(inner: &ffi::Trajectory, time: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let dofs = ffi::traj_get_dofs(inner);
    let mut pos = vec![0.0f64; dofs];
    let mut vel = vec![0.0f64; dofs];
    let mut acc = vec![0.0f64; dofs];
    unsafe {
        ffi::traj_at_time(inner, time, pos.as_mut_ptr(), vel.as_mut_ptr(), acc.as_mut_ptr(), dofs);
    }
    (pos, vel, acc)
}

impl Trajectory {
    pub fn new_direct(dofs: usize) -> Self {
        Self(ffi::traj_new_direct(dofs))
    }

    pub fn new(dofs: usize, waypoints: usize) -> Self {
        Self(ffi::traj_new(dofs, waypoints))
    }

    pub fn get_duration(&self) -> f64 {
        ffi::traj_get_duration(&self.0)
    }

    pub fn at_time(&self, time: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        traj_at_time_impl(&self.0, time)
    }
}

impl<'a> TrajectoryRef<'a> {
    pub fn get_duration(&self) -> f64 {
        ffi::traj_get_duration(self.inner)
    }

    pub fn at_time(&self, time: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        traj_at_time_impl(self.inner, time)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum ControlInterface {
    #[default]
    Position = 0,
    Velocity = 1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum Synchronization {
    #[default]
    Time = 0,
    TimeIfNecessary = 1,
    Phase = 2,
    None = 3,
}

#[repr(i32)]
#[derive(Debug, PartialEq, Eq)]
pub enum Result {
    Working = 0,
    Finished = 1,
    Error = -1,
    ErrorInvalidInput = -100,
    ErrorTrajectoryDuration = -101,
    ErrorPositionalLimits = -102,
    ErrorZeroLimits = -104,
    ErrorExecutionTimeCalculation = -110,
    ErrorSynchronizationCalculation = -111,
}

impl From<i32> for Result {
    fn from(v: i32) -> Self {
        match v {
            0    => Result::Working,
            1    => Result::Finished,
            -100 => Result::ErrorInvalidInput,
            -101 => Result::ErrorTrajectoryDuration,
            -102 => Result::ErrorPositionalLimits,
            -104 => Result::ErrorZeroLimits,
            -110 => Result::ErrorExecutionTimeCalculation,
            -111 => Result::ErrorSynchronizationCalculation,
            _    => Result::Error,
        }
    }
}

pub struct InputParameter {
    pub degrees_of_freedom: usize,
    pub control_interface: ControlInterface,
    pub synchronization: Synchronization,
    pub current_position: Vec<f64>,
    pub current_velocity: Vec<f64>,
    pub current_acceleration: Vec<f64>,
    pub target_position: Vec<f64>,
    pub target_velocity: Vec<f64>,
    pub target_acceleration: Vec<f64>,
    pub max_velocity: Vec<f64>,
    pub max_acceleration: Vec<f64>,
    pub max_jerk: Vec<f64>,
    pub max_position: Vec<f64>,
    pub min_position: Vec<f64>,
    pub intermediate_positions: Vec<Vec<f64>>,
    inner: UniquePtr<ffi::InputParameter>,
}

impl InputParameter {
    pub fn new(dofs: usize) -> Self {
        Self {
            degrees_of_freedom: dofs,
            control_interface: ControlInterface::Position,
            synchronization: Synchronization::Time,
            current_position: vec![0.0; dofs],
            current_velocity: vec![0.0; dofs],
            current_acceleration: vec![0.0; dofs],
            target_position: vec![0.0; dofs],
            target_velocity: vec![0.0; dofs],
            target_acceleration: vec![0.0; dofs],
            max_velocity: vec![0.0; dofs],
            max_acceleration: vec![f64::INFINITY; dofs],
            max_jerk: vec![f64::INFINITY; dofs],
            max_position: vec![f64::INFINITY; dofs],
            min_position: vec![f64::NEG_INFINITY; dofs],
            intermediate_positions: Vec::new(),
            inner: ffi::input_new(dofs),
        }
    }

    fn sync(&mut self) {
        self.control_interface = unsafe { std::mem::transmute(ffi::input_get_control_interface(&self.inner)) };
        self.synchronization   = unsafe { std::mem::transmute(ffi::input_get_synchronization(&self.inner)) };

        self.current_position     = ffi::input_get_current_position(&self.inner).iter().cloned().collect();
        self.current_velocity     = ffi::input_get_current_velocity(&self.inner).iter().cloned().collect();
        self.current_acceleration = ffi::input_get_current_acceleration(&self.inner).iter().cloned().collect();
        self.target_position      = ffi::input_get_target_position(&self.inner).iter().cloned().collect();
        self.target_velocity      = ffi::input_get_target_velocity(&self.inner).iter().cloned().collect();
        self.target_acceleration  = ffi::input_get_target_acceleration(&self.inner).iter().cloned().collect();
        self.max_velocity         = ffi::input_get_max_velocity(&self.inner).iter().cloned().collect();
        self.max_acceleration     = ffi::input_get_max_acceleration(&self.inner).iter().cloned().collect();
        self.max_jerk             = ffi::input_get_max_jerk(&self.inner).iter().cloned().collect();
        self.max_position         = ffi::input_get_max_position(&self.inner).iter().cloned().collect();
        self.min_position         = ffi::input_get_min_position(&self.inner).iter().cloned().collect();

        let cap = self.degrees_of_freedom * 256;
        let mut flat = vec![0.0f64; cap];
        let written = unsafe { ffi::input_get_intermediate_positions(&self.inner, flat.as_mut_ptr(), cap) };
        self.intermediate_positions = if self.degrees_of_freedom > 0 && written > 0 {
            flat[..written].chunks(self.degrees_of_freedom).map(|c| c.to_vec()).collect()
        } else {
            Vec::new()
        };
    }

    fn apply(&mut self) {
        ffi::input_set_control_interface(self.inner.pin_mut(), self.control_interface as i32);
        ffi::input_set_synchronization(self.inner.pin_mut(),   self.synchronization as i32);

        unsafe {
            ffi::input_set_current_position(self.inner.pin_mut(),     self.current_position.as_ptr(),     self.current_position.len());
            ffi::input_set_current_velocity(self.inner.pin_mut(),     self.current_velocity.as_ptr(),     self.current_velocity.len());
            ffi::input_set_current_acceleration(self.inner.pin_mut(), self.current_acceleration.as_ptr(), self.current_acceleration.len());
            ffi::input_set_target_position(self.inner.pin_mut(),      self.target_position.as_ptr(),      self.target_position.len());
            ffi::input_set_target_velocity(self.inner.pin_mut(),      self.target_velocity.as_ptr(),      self.target_velocity.len());
            ffi::input_set_target_acceleration(self.inner.pin_mut(),  self.target_acceleration.as_ptr(),  self.target_acceleration.len());
            ffi::input_set_max_velocity(self.inner.pin_mut(),         self.max_velocity.as_ptr(),         self.max_velocity.len());
            ffi::input_set_max_acceleration(self.inner.pin_mut(),     self.max_acceleration.as_ptr(),     self.max_acceleration.len());
            ffi::input_set_max_jerk(self.inner.pin_mut(),             self.max_jerk.as_ptr(),             self.max_jerk.len());
            ffi::input_set_max_position(self.inner.pin_mut(),         self.max_position.as_ptr(),         self.max_position.len());
            ffi::input_set_min_position(self.inner.pin_mut(),         self.min_position.as_ptr(),         self.min_position.len());
        }

        if self.intermediate_positions.is_empty() {
            ffi::input_clear_intermediate_positions(self.inner.pin_mut());
        } else {
            let dofs = self.intermediate_positions[0].len();
            let flat: Vec<f64> = self.intermediate_positions.iter().flat_map(|r| r.iter().cloned()).collect();
            unsafe { ffi::input_set_intermediate_positions(self.inner.pin_mut(), flat.as_ptr(), dofs, self.intermediate_positions.len()); }
        }
    }
}

pub struct OutputParameter {
    pub degrees_of_freedom: usize,
    pub new_position: Vec<f64>,
    pub new_velocity: Vec<f64>,
    pub new_acceleration: Vec<f64>,
    pub time: f64,
    pub new_section: usize,
    pub new_calculation: bool,
    pub calculation_duration: f64,
    inner: UniquePtr<ffi::OutputParameter>,
}

impl OutputParameter {
    pub fn new_direct(dofs: usize) -> Self {
        Self {
            degrees_of_freedom: dofs,
            new_position: vec![0.0; dofs],
            new_velocity: vec![0.0; dofs],
            new_acceleration: vec![0.0; dofs],
            time: 0.0,
            new_section: 0,
            new_calculation: false,
            calculation_duration: 0.0,
            inner: ffi::output_new_direct(dofs),
        }
    }

    pub fn new(dofs: usize, waypoints: usize) -> Self {
        Self {
            degrees_of_freedom: dofs,
            new_position: vec![0.0; dofs],
            new_velocity: vec![0.0; dofs],
            new_acceleration: vec![0.0; dofs],
            time: 0.0,
            new_section: 0,
            new_calculation: false,
            calculation_duration: 0.0,
            inner: ffi::output_new(dofs, waypoints),
        }
    }

    pub fn pass_to_input(&self, input: &mut InputParameter) {
        ffi::output_pass_to_input(&self.inner, input.inner.pin_mut());
        input.sync();
    }

    fn sync(&mut self) {
        self.new_position         = ffi::output_get_new_position(&self.inner).iter().cloned().collect();
        self.new_velocity         = ffi::output_get_new_velocity(&self.inner).iter().cloned().collect();
        self.new_acceleration     = ffi::output_get_new_acceleration(&self.inner).iter().cloned().collect();
        self.time                 = ffi::output_get_time(&self.inner);
        self.new_section          = ffi::output_get_new_section(&self.inner);
        self.new_calculation      = ffi::output_get_new_calculation(&self.inner);
        self.calculation_duration = ffi::output_get_calculation_duration(&self.inner);
    }

    pub fn trajectory(&self) -> TrajectoryRef<'_> {
        TrajectoryRef { inner: ffi::output_get_trajectory(&self.inner) }
    }
}

pub struct Ruckig(UniquePtr<ffi::Ruckig>);

impl Ruckig {
    pub fn new_direct_and_offline(dofs: usize) -> Self {
        Self(ffi::ruckig_new_direct_and_offline(dofs))
    }

    pub fn new_direct(dofs: usize, delta_time: f64) -> Self {
        Self(ffi::ruckig_new_direct(dofs, delta_time))
    }

    pub fn new_offline(dofs: usize, waypoints: usize) -> Self {
        Self(ffi::ruckig_new_offline(dofs, waypoints))
    }

    pub fn new(dofs: usize, delta_time: f64, waypoints: usize) -> Self {
        Self(ffi::ruckig_new(dofs, delta_time, waypoints))
    }

    pub fn reset(&mut self) {
        ffi::ruckig_reset(self.0.pin_mut());
    }

    pub fn update(&mut self, input: &mut InputParameter, output: &mut OutputParameter) -> Result {
        input.apply();
        let raw = ffi::ruckig_update(self.0.pin_mut(), &input.inner, output.inner.pin_mut());
        output.sync();
        Result::from(raw)
    }

    pub fn calculate(&mut self, input: &mut InputParameter, trajectory: &mut Trajectory) -> Result {
        input.apply();
        let raw = ffi::ruckig_calculate(self.0.pin_mut(), &input.inner, trajectory.0.pin_mut());
        Result::from(raw)
    }
}
