cmake_minimum_required(VERSION 3.10)


project(ruckig VERSION 0.4.0 LANGUAGES CXX)


list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(GNUInstallDirs)


option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_PYTHON_MODULE "Build python module" OFF)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_BENCHMARK "Build benchmark" OFF)
option(BUILD_SHARED_LIBS "Build as shared library" ON)

if(WIN32 AND BUILD_SHARED_LIBS)
  option(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "On Windows, export all symbols when building a shared library." ON)
endif()


find_package(Reflexxes QUIET)


add_library(ruckig
  src/brake.cpp
  src/position-step1.cpp
  src/position-step2.cpp
  src/velocity-step1.cpp
  src/velocity-step2.cpp
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
  target_compile_options(ruckig PRIVATE /W4 /WX)

else()
  target_compile_options(ruckig PRIVATE -Werror -Wall -Wextra)
endif()


if(Reflexxes)
  set(REFLEXXES_TYPE "ReflexxesTypeII" CACHE STRING "Type of Reflexxes library") # or ReflexxesTypeIV
  message("Found Reflexxes ${REFLEXXES_TYPE}")

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
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
    add_subdirectory(pybind11)
  else()
    find_package(pybind11 2.6 REQUIRED)
  endif()

  pybind11_add_module(python_ruckig src/python.cpp)
  target_compile_features(python_ruckig PUBLIC cxx_std_17)
  target_link_libraries(python_ruckig PUBLIC ruckig)

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

  foreach(test IN ITEMS otg-test)
    add_executable(${test} "test/${test}.cpp")
    if(Reflexxes)
      target_compile_definitions(${test} PUBLIC WITH_REFLEXXES)
    endif()
    target_link_libraries(${test} PRIVATE ruckig)
    add_test(NAME ${test} COMMAND ${test})
  endforeach()

  if(BUILD_BENCHMARK)
    add_executable(otg-benchmark "test/otg-benchmark.cpp")
    if(Reflexxes)
      target_compile_definitions(otg-benchmark PUBLIC WITH_REFLEXXES)
    endif()
    target_link_libraries(otg-benchmark PRIVATE ruckig)
  endif()
endif()
