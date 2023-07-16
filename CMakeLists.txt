cmake_minimum_required(VERSION 3.24)
project(FluidSim)

set(CMAKE_CXX_STANDARD 17)

add_executable(FluidSim main.cpp Fluid.cpp Fluid.h Utils.h Projection.cpp Advection.cpp TransferToParticles.cpp TransferToGrid.cpp)


find_package(SFML 2.5 COMPONENTS system window graphics network audio REQUIRED)
include_directories(${SFML_INCLUDE_DIRS})
target_link_libraries(FluidSim sfml-system sfml-window sfml-graphics)
