cmake_minimum_required(VERSION 3.24)
project(FluidSim)

set(CMAKE_CXX_STANDARD 17)

add_executable(FluidSim main.cpp fluid.cpp fluid.h projection.cpp advection.cpp transfer_to_particles.cpp transfer_to_grid.cpp particles.h particles.cpp grid.h grid.cpp)


find_package(SFML 2.5 COMPONENTS system window graphics network audio REQUIRED)
include_directories(${SFML_INCLUDE_DIRS})
target_link_libraries(FluidSim sfml-system sfml-window sfml-graphics)
