^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package ruckig
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Forthcoming
-----------
* [ci] don't check python for c++11
* [ci] move gcc-5 to ubuntu 20.04
* bump version
* add example for per section minimum duration
* Budan's theorem
* use valid profile iterator in step1
* add jerk-minimizing profiles in step2
* fix c++11 patch
* fix c++11 patch
* fix patch c++11
* c++11 fix inplace patching
* improve support for c++11
* update pybind11 in ci
* fix optional in tests with c++11 patch
* update ci checkout version
* fix msvc warning in online calculator
* fix msvc compiler warnings
* expose calculators
* bump version
* performance improvement in step2 vel
* Contributors: pantor

0.8.4 (2022-09-13)
------------------
* robustify step2 vel uddu
* fix readme tracking example
* add tracking interface to readme
* fix brake duration reset
* fix brake when current acceleration is on limit
* improve tracking example
* fix pyproject.toml
* clean tracking example
* clean examples
* add pyproject file
* add DOFs to traj serialization
* nlohmann/json based serialization for trajectory
* point out Eigen 3.4 or later
* add custom vector types to readme
* disable vector types example with online client
* include deque
* add examples for custom vector types
* dont support custom data types with online client
* add tests and fixes for custom vector type
* use template specalization without using
* custom vector type as template template parameter
* clean error_at_max_duration, add StandardVector
* clean licenses
* update header dependencies
* clean comparison benchmarks
* extend phase synchronization to velocity control
* move src to src/ruckig
* clean gitignore
* remove -Werror
* Contributors: pantor

0.7.1 (2022-07-10)
------------------
* bump to 0.7.1
* fix python 3.10 deployment
* Merge branch 'master' of github.com:pantor/ruckig
* bump to 0.7.0
* allow user to force new computation for the same input (`#132 <https://github.com/pantor/ruckig/issues/132>`_)
  * allow user to force computation again
  Otherwise sending the same target with non-zero min_duration
  multiple times in a control loop will not trigger a new trajectory.
  * rename to reset
  Co-authored-by: pantor <lars.berscheid@online.de>
* profile precision as constant
* clean notebooks
* fix c++11 std::clamp
* use steady_clock, minor code cleaning of roots.hpp
* fix numeric for time sync with discrete durations
* fix independent min duration with brake trajectory
* clean header includes
* improve stability of velocity control
* msvc warnings
* fix msvc warnings
* static cast for std::div
* Merge branch 'master' of github.com:pantor/ruckig
* fix warnings, limiting_dof to std::optional
* doc: fixed used-by entry (`#131 <https://github.com/pantor/ruckig/issues/131>`_)
* move roots into ruckig namespace
* refactor tests and benchmarks
* step2 vel catch root at boundary
* update benchmark plot
* add algorithm header, use c++ array
* make delta_time non const
* add smoothing of non-limiting dofs
* const degrees_of_freedom
* print error messages from online api
* check that current_input is initialized
* Contributors: Matthias Seehauser, Michael Görner, pantor

0.6.5 (2022-03-06)
------------------
* revert to manylinux2010_x86_64
* bump to v0.6.5
* fix pypi build from source
* remove 9_trajectory.png
* use .pdf for example trajectories
* numerical stability in velocity control
* fix for zero-to-zero traj with velocity control
* invalidate discrete duration
* fix numerical instability in acc0
* update to nlohmann/json v3.10.5
* bump to 0.6.4
* clarify Online API in examples and Readme
* fix example docs
* fix ci
* build waypoints example only for online client
* add joint example for dynamic dofs and waypoints
* fix capacity / actual waypoints mismatch
* disable ros cd
* retry
* retry
* retry
* debug publish ros
* Contributors: pantor

0.6.3 (2022-01-21)
------------------
* bump to v0.6.3
* activaten open_pr for bloom release
* publish ros on release
* test bloom ci
* add bloom release
* several perf optimizations, step2 stability
* clear up waypoints in readme
* fix time sync for discrete durations, rename step1
* Contributors: pantor

0.6.0 (2021-12-06)
------------------
* fix python upload
* bump version to 0.6.0
* filter_intermediate_positions
* add braces to if function
* fix error in step 2
* remove filter positions
* remote api for intermediate waypoints
* fix trajectories with zero duration
* use integrated instead target values after traj
* use back() instead of explicit number
* ci build matrix
* BUILD_ONLINE_CLIENT in python package
* add brake in online calculator
* fix ci online-calculator
* auto ci name
* add online calculator for intermediate waypoints
* add httplib and build option
* create third_party directory
* update demo notebook
* update notebook demo
* add jupyter demo notebook
* change brief description of calculator
* expose internal
* add note to ruckig pro examples
* clear up section term
* clean brake class
* refactor integrate to utils
* prepare accel phase
* use dynamic dofs const
* improve input validation
* clean up CD
* Contributors: pantor

0.5.0 (2021-11-14)
------------------
* debug publish
* publish pypi package on released release
* bump version
* add hint for Step 1 None
* optimize block class
* improve readme
* per_section_limits in readme
* add per_section_limits
* fix c+11 patch
* fix per dof synchronization
* split CMakeLists examples
* fix per dof synchronization with none
* separate trajectory and calculator class
* fix windows  build
* code cleaning
* intermediate waypoints readme
* fix number of waypoints
* avoid pow completely
* update maintainer
* simplify readme
* use brake profile class
* fix finished trajectory for disabled dofs
* minor code cleaning
* Merge branch 'master' of github.com:pantor/ruckig
* add to_string for output parameters (`#77 <https://github.com/pantor/ruckig/issues/77>`_)
* add ref to vel add_profile
* positional limits
* min limits for intermediate positions
* extend phase synchronization
* performance improvements
* add section to output in readme
* pass_to_input, did_section_change
* fix nullopt with c++11 patch
* fix nullopt in c++11 patch
* fix c++11
* per dof control interface and synchronization `#53 <https://github.com/pantor/ruckig/issues/53>`_
* add section index to output
* Notes about Intermediate Waypoints
* interrupt calculation duration in microseconds
* add ruckig pro examples
* add ruckig toppra comparison
* improve readme and examples
* introduce Ruckig Pro
* remove generate docs ci
* link docs to ruckig.com
* add examples doc pages
* fix example names
* add more examples
* Instantaneous Motion Generation
* add calculation interruption
* add doxyfile again
* step1: numeric stability
* Contributors: Lars Berscheid, lpdon, pantor

0.4.0 (2021-08-23)
------------------
* update version to 0.4
* code cleaning
* add was_calculation_interrupted
* step 1: performance optimizations
* interrupt_calculation_duration
* Add CITATION.cff file (`#63 <https://github.com/pantor/ruckig/issues/63>`_)
  * add CITATION.cff
  * CITATOION wip
  * fix cite
* Update README.md
* update to doctest 2.4.6
* code cleaning
* performance optimizations
* step 2: performance optimization
* step 2: performance optimization
* performance optimization: step 1
* performance optimization: positive set of roots
* performance optimization in step1
* code cleaning
* set find_package reflexxes to quiet
* remove BSD license text, why was it there anyway?
* clean plot trajectory
* add printing intermediate_positions
* code cleaning
* add degrees_of_freedom to python trajectory
* rename interface to control_interface
* set_limits for ACC1, code cleaning
* improve numeric stability
* in vel interface, change acc threshold
* code cleanup
* add DynamicDOFs constant
* numerical stability of velocity interface
* improve numerical stability
* fix variable name redeclaration
* fix jerk violation in step2, some optimizations
* clean up check methods of profile
* fix min_velocity with phas e synchronization
* fix phase synchronization in python
* improve numerical stability for high jerks
* fix newton step in acc0/acc1 for t=0
* clean up plot_trajectory
* validate for isnan
* fix python path in examples
* fix position extrema for dynamic number of dofs
* Added python example for non-realtime context (`#43 <https://github.com/pantor/ruckig/issues/43>`_)
  * added python example for non-realtime context
  * fixed ci workflow
  * need to reset t_start
  * rename
  * change example
  Co-authored-by: Max Dhom <md@adva.store>
  Co-authored-by: Lars Berscheid <lars.berscheid@kit.edu>
* Dynamic Dofs (`#47 <https://github.com/pantor/ruckig/issues/47>`_)
  * vectorize
  * use vector for dof==0
  * add tests and checks
  * include vector
  * default dofs = 0, fix reflexxes
  * redo default dofs template parameter
  * add readme
  * improve readme
  * fix grammar
* add tests for invalid input
* add offline trajectory generation
* add paren to get time at position
* add get_first_time_at_position method
* Contributors: Lars Berscheid, Mathis, Max Dhom, pantor

0.3.3 (2021-06-25)
------------------
* Merge branch 'master' of github.com:pantor/ruckig
* version 0.3.3
* Set CMAKE_OUTPUT_DIRECTORY for windows build (`#41 <https://github.com/pantor/ruckig/issues/41>`_)
  * check windows
  * add library
  * add debug log
  * try to fix
  * try to fix2
  * try to fix 3
  * try to fix 4
  * fix 5
  * rest test number
  * fix setup.py
  * remove log
* hard-code build directory of python library
* fix windows packge, version 0.3.2
* pre-compiled packages for windows on pypi
* Contributors: Lars Berscheid, pantor

0.3.1 (2021-06-24)
------------------
* set version 0.3.1
* add manifest.in
* double newton step in step2 vel udud
* Fix windows python package (`#38 <https://github.com/pantor/ruckig/issues/38>`_)
  * fix windows
* update benchmark figure
* c++11 dont treat warnings as errors
* fix three step
* performance improvements
* vectorize dof
* Fix Patch for C++11 (`#36 <https://github.com/pantor/ruckig/issues/36>`_)
  * add ci for c++11
  * remove maybe_unused
  * patch in-place
  * fix c++11
  * replace make_unique
  * find error in ci
  * try to fix gcc-5
  * dont build python
  * dont patch cpp files
  * deactivate cmake flags in patch script
  * test python example
* add C++11 patch in readme
* add patch script for C++11
* Contributors: Lars Berscheid, pantor

0.3.0 (2021-06-16)
------------------
* update version number
* add python at_time comment for doxygen
* check for v_min, fix directional tests
* python return at_time method
* fix max target acceleration
* fix and test extremal positions
* synchronize function to trajectory class
* fix max target acceleration when min_velocity
* clean up docs
* fixed bug in check_position_extremum (`#33 <https://github.com/pantor/ruckig/issues/33>`_)
  * fixed copy & paste error in Trajectory::is_phase_synchronizable
  * fixed obvious copy & paste error in check_position_extremum
* fixed copy & paste error in Trajectory::is_phase_synchronizable (`#32 <https://github.com/pantor/ruckig/issues/32>`_)
* fix negative, near zero pd cases
* Update README.md
* Update README.md
* Update README.md
* check limits with equal sign for numeric stability
* copy jerk_signs for phase synchronization
* remove alternative otgs
* add citation and link to paper
* remove alternative otgs from python
* fix numerical issues on time-optimal traj
* Merge branch 'master' of github.com:pantor/ruckig
* fix numeric error for long durations
* Scale to 5e9 random trajectories
* fix numerical issues
* add step through tests
* fix for braking
* double newton step
* fix macos and windows build
* recalculate after error
* fix some numerical issues
* use checkout@v2 in ci
* use ARCHIVE_OUTPUT_NAME for python wrapper
* msvc: warning as errors
* fix msvc compiler warnings
* Update README.md
* Update README.md
* Fix BUILD_PYTHON_MODULE option on Windows/MSVC (`#18 <https://github.com/pantor/ruckig/issues/18>`_)
* fix ci args
* add phase synchronization
* set boundary method
* simplify python example
* Add pip install to Readme
* Contributors: Lars Berscheid, Silvio Traversaro, pantor, stefanbesler

0.2.6 (2021-03-29)
------------------
* remove env
* fix python cd
* use static library
* test python executable
* fix python package
* python cd: fix windows
* add setup.py, version 0.2.1
* rename python module to ruckig
* generate python classes for multiple dof
* Add a ROS package manifest (`#10 <https://github.com/pantor/ruckig/issues/10>`_)
  Enables building Ruckig as a plain CMake package in a Catkin/Colcon workspace.
* fix end of brake trajectory integration
* Include GNU install dirs earlier (`#11 <https://github.com/pantor/ruckig/issues/11>`_)
  Otherwise 'CMAKE_INSTALL_INCLUDEDIR' is empty/undefined when it's used to setup ruckig::ruckig's target_include_directories(..).
* readme: some minor typos (`#9 <https://github.com/pantor/ruckig/issues/9>`_)
  I happened to notice them.
* privatize trajectory class members
* privatize some class members
* use cmath
* code cleaning
* show enum in docs
* split parameter files, calculate in trajectory
* code cleaning
* add python example
* Move options to API documentation
* Fix Readme Code
* fix undefined output for zero duration
* indicate default values in Readme
* add discrete durations
* Add Windows and macOS build to CI (`#4 <https://github.com/pantor/ruckig/issues/4>`_)
  * windows and mac ci
  * use cmake action for generator
  * fix ci build directory
  * run tests only on linux
  * test example
* Add support to compile on Windows (`#3 <https://github.com/pantor/ruckig/issues/3>`_)
* Merge pull request `#2 <https://github.com/pantor/ruckig/issues/2>`_ from traversaro/patch-1
  Use BUILD_SHARED_LIBS to select if compile C++ library as static or shared
* document examples/CMakeLists.txt
* add ci to PRs
* Use BUILD_SHARED_LIBS to select if compile as static or shared
* Merge pull request `#1 <https://github.com/pantor/ruckig/issues/1>`_ from traversaro/add-install-support
  Add support for installation of the C++ library
* Add support for installation of the C++ library
* set correct cmake and doxygen version
* fix more edge cases
* fix some edge cases
* document velocity interface
* added velocity interface
* improve readme
* improve benchmark
* fix finding subdirectory
* check ci
* build python module in ci
* use own set class on stack
* fix synchronization enum, better python support
* add time synchronization parameter
* fix motion finished reset
* fix immediate reaction
* fix more edge cases
* refine min acceleration
* add min_acceleration
* fix some edge cases
* Merge branch 'master' of github.com:pantor/ruckig
* decrease required cmake to 3.10
* fix ci
* introduce trajectory class
* position extrema
* scale tests to 1e9
* several optimizations
* compile with warnings
* step2: code cleaning
* fix numeric edge cases
* move test suite to doctest
* fix cmake as submodule
* Merge branch 'master' of github.com:pantor/ruckig
* fix optional minimum time
* code documentation, more tests
* fix benchmark
* build benchmark in ci
* fix gcc
* remove eigen dep in cmake
* code cleaning
* code cleaning
* add comparison with multiple DoF
* rix ci: braking
* fix interval selection
* clean braking
* fix ci, more tests
* add benchmark, numeric stability
* fix block in step1
* clean root finding
* Increase stability
* improve tests, remove eigen
* code style
* increase test coverage
* add min_velocity
* fix ci
* Step1: Use Newton step
* fix ci
* fix time sync
* update tests
* more tests
* add example
* more tests
* increase number of tests
* simplify equations
* simplify equations
* simplify equations
* add tuple header
* further code cleaning
* remove eigen dependency, code cleaning
* clean brake code
* code cleaning
* code cleaning
* add doxygen
* improve acceleration target, readme
* refine max time deviation to 1e-8
* code cleaning
* improve time sync
* block synchronization
* stability for target acceleration in 1 dof
* update readme for ruckig
* add license
* improve tests
* fix eigen ci folder name
* fix eigen git repository
* fix 1dof vf comparison
* remove complex algorithmic
* code cleaning
* initial commit
* Contributors: G.A. vd. Hoorn, Lars Berscheid, Silvio Traversaro, pantor
