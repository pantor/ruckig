cmake_minimum_required(VERSION 3.11)


project(ruckig VERSION 0.1.0 LANGUAGES CXX)
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)


option(BUILD_EXAMPLES "Build example programs" OFF)
option(BUILD_PYTHON_MODULE "Build python module" ON)
option(BUILD_TESTS "Build tests" ON)


find_package(Eigen3 3.3.9 REQUIRED NO_MODULE)
find_package(Reflexxes)

message("Found Eigen Version: ${Eigen3_VERSION}")

if(Reflexxes)
  set(REFLEXXES_TYPE "ReflexxesTypeII" CACHE STRING "Type of Reflexxes library") # or ReflexxesTypeIV

  message("Found Reflexxes ${REFLEXXES_TYPE}")
else()
  message("Did not found Reflexxes.")
endif()


add_library(ruckig SHARED
  src/brake.cpp
  src/step1.cpp
  src/step2.cpp
)
target_compile_features(ruckig PUBLIC cxx_std_17)
target_include_directories(ruckig PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ruckig PUBLIC Eigen3::Eigen)


if(Reflexxes)
  target_compile_definitions(ruckig PUBLIC WITH_REFLEXXES)
  target_link_libraries(ruckig PUBLIC Reflexxes::Reflexxes)
endif()


if(BUILD_EXAMPLES)
  foreach(example IN ITEMS position)
    add_executable(example-${example} "examples/${example}.cpp")
    target_link_libraries(example-${example} PRIVATE ruckig)
  endforeach()
endif()


if(BUILD_PYTHON_MODULE)
  # Check if pybind11 exists as a subdirectory
  if(EXISTS pybind11)
    add_subdirectory(pybind11)
  else()
    find_package(Python COMPONENTS Interpreter Development)
    find_package(pybind11 2.6 REQUIRED)
  endif()

  pybind11_add_module(_ruckig src/python.cpp)
  target_compile_features(_ruckig PUBLIC cxx_std_17)
  target_link_libraries(_ruckig PUBLIC ruckig)
endif()


if(BUILD_TESTS)
  enable_testing()

  find_package(Catch2 REQUIRED)

  foreach(test IN ITEMS otg-test)
    add_executable(${test} "test/${test}.cpp")
    if(Reflexxes)
      target_compile_definitions(${test} PUBLIC WITH_REFLEXXES)
    endif()
    target_link_libraries(${test} PRIVATE ruckig Catch2::Catch2)
    add_test(NAME ${test} COMMAND ${test})
  endforeach()
endif()
