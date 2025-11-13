cmake_minimum_required(VERSION 3.16)
project(TrajectorySimplification CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(main
  main.cpp
)

# Add include path for headers in inc/
target_include_directories(main PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/inc
)