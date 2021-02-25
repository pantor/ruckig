cmake_minimum_required(VERSION 3.10)

project(ruckig_examples)

find_package(ruckig REQUIRED)

add_executable(ruckig-example-position position.cpp)
target_link_libraries(ruckig-example-position PRIVATE ruckig::ruckig)
