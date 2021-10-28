cmake_minimum_required(VERSION 3.10)

project(ruckig_examples)

# (1) When installed system-wide
find_package(ruckig REQUIRED)

# (2) Or find the header-files and library
include_directories(include)
link_directories(lib)

# Build the position example
add_executable(example-position 1_position.cpp)
target_compile_features(example-position PUBLIC cxx_std_17)

# For (1)
target_link_libraries(example-position PRIVATE ruckig::ruckig)

# For (2)
target_link_libraries(example-position ruckig)