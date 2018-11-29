cmake_minimum_required(VERSION 3.9)
project(threebody)

set(CMAKE_CXX_STANDARD 11)

add_executable(threebody main.cpp SolarSystem.cpp SolarSystem.h Coords.cpp Coords.h)