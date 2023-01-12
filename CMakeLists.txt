cmake_minimum_required(VERSION 3.10)


project(ruckig VERSION 0.9.2 LANGUAGES CXX)


include(GNUInstallDirs)


option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_PYTHON_MODULE "Build python module with pybind11" OFF)
option(BUILD_ONLINE_CLIENT "Build online client to calculate Ruckig Pro trajectories" OFF)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_BENCHMARK "Build benchmark" OFF)
option(BUILD_SHARED_LIBS "Build as shared library" ON)
option(EXPOSE_INTERNAL "Expose some internal classes and methods" OFF)

if(WIN32 AND BUILD_SHARED_LIBS)
  option(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "On Windows, export all symbols when building a shared library." ON)
endif()


add_library(ruckig
  src/ruckig/brake.cpp
  src/ruckig/position-step1.cpp
  src/ruckig/position-step2.cpp
  src/ruckig/velocity-step1.cpp
  src/ruckig/velocity-step2.cpp
)
add_library(ruckig::ruckig ALIAS ruckig)

target_compile_features(ruckig PUBLIC cxx_std_17)
target_include_directories(ruckig PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ruckig PUBLIC)


if(MSVC)
  target_compile_definitions(ruckig PUBLIC _USE_MATH_DEFINES)
endif()

if(EXPOSE_INTERNAL)
  target_compile_definitions(ruckig PUBLIC WITH_EXPOSE_INTERNAL)
endif()


if(BUILD_ONLINE_CLIENT)
  target_include_directories(ruckig PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  target_compile_definitions(ruckig PUBLIC WITH_ONLINE_CLIENT WITH_SERIALIZATION)
endif()


if(BUILD_EXAMPLES)
  set(EXAMPLES_LIST 1_position 2_position_offline 5_velocity 6_stop 7_minimum_duration 9_dynamic_dofs)

  if(BUILD_ONLINE_CLIENT)
    list(APPEND EXAMPLES_LIST 3_waypoints 4_waypoints_online 8_per_section_minimum_duration 10_dynamic_dofs_waypoints)
  else()
    list(APPEND EXAMPLES_LIST 12_custom_vector_type 13_custom_vector_type_dynamic_dofs)
    if(TARGET Eigen3::Eigen)
      list(APPEND EXAMPLES_LIST 11_eigen_vector_type)
    endif()
  endif()

  foreach(example IN LISTS EXAMPLES_LIST)
    add_executable(example-${example} examples/${example}.cpp)
    target_link_libraries(example-${example} PRIVATE ruckig)
  endforeach()
endif()


if(BUILD_PYTHON_MODULE)
  # Check if pybind11 exists as a subdirectory
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
    add_subdirectory(pybind11)
  else()
    find_package(pybind11 2.6 REQUIRED)
  endif()

  pybind11_add_module(python_ruckig src/ruckig/python.cpp)
  target_compile_features(python_ruckig PUBLIC cxx_std_17)
  target_link_libraries(python_ruckig PUBLIC ruckig)
  if(BUILD_ONLINE_CLIENT)
    target_compile_definitions(python_ruckig PUBLIC WITH_ONLINE_CLIENT)
  endif()

  set_target_properties(python_ruckig PROPERTIES OUTPUT_NAME ruckig)
  set_target_properties(python_ruckig PROPERTIES ARCHIVE_OUTPUT_NAME python_ruckig)
endif()


# Add support for installation
include(CMakePackageConfigHelpers)

# Install headers
install(DIRECTORY include/ruckig DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install library
install(TARGETS ruckig
  EXPORT ${PROJECT_NAME}-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install CMake config files
set(ruckig_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

install(EXPORT ${PROJECT_NAME}-targets
  DESTINATION ${ruckig_INSTALL_CONFIGDIR}
  NAMESPACE ruckig::
)

configure_package_config_file(cmake/ruckig-config.cmake.in ruckig-config.cmake
  INSTALL_DESTINATION ${ruckig_INSTALL_CONFIGDIR}
)

write_basic_package_version_file(ruckig-config-version.cmake
  VERSION ${ruckig_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/ruckig-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/ruckig-config-version.cmake"
  DESTINATION ${ruckig_INSTALL_CONFIGDIR}
)

install(FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/package.xml"
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}
)


if(BUILD_TESTS)
  enable_testing()

  add_executable(test-target test/test-target.cpp test/test-target-known.cpp)

  target_link_libraries(test-target PRIVATE ruckig)
  target_include_directories(test-target PRIVATE third_party)
  if(MSVC)
    target_compile_options(test-target PRIVATE /W4)
  else()
    target_compile_options(test-target PRIVATE -Wall -Wextra)
  endif()

  add_test(NAME test-target COMMAND test-target)
endif()


if(BUILD_BENCHMARK)
  add_executable(benchmark-target test/benchmark-target.cpp)

  target_link_libraries(benchmark-target PRIVATE ruckig)
endif()
