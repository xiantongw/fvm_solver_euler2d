cmake_minimum_required (VERSION 3.0)

project(solver.exe)

include_directories(./include)

file(GLOB SRC_FILES "./src/*.cpp")

add_executable(${PROJECT_NAME} ${SRC_FILES})

add_definitions(-std=c++11)

set(CMAKE_CXX_COMPILER      "g++" )
set(CMAKE_CXX_FLAGS_DEBUG   "-O3" )