cmake_minimum_required(VERSION 3.16)
project(simplification_of_trajectory_stream CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()

# Headless simplify binary — no Qt dependency. Web-trace emission lives in
# web_trace.cpp, kept separate from simplify.cpp per refactor.
find_package(CGAL REQUIRED)
add_executable(simplify
  simplify.cpp
  web_trace.cpp
)
target_link_libraries(simplify
  CGAL::CGAL
)

# Explicitly requested only: a serial, long-lived sampling-profiler workload.
add_executable(profile_sots EXCLUDE_FROM_ALL scripts/profile_sots.cpp)
target_include_directories(profile_sots PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(profile_sots CGAL::CGAL)
target_compile_options(profile_sots PRIVATE -g)

add_executable(test_clip_batch EXCLUDE_FROM_ALL scripts/test_clip_batch.cpp)
target_include_directories(test_clip_batch PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(test_clip_batch CGAL::CGAL)

# Baseline binaries used by the web compare pane and benchmark.py.
# dp / squish need no Qt. dots needs Qt6 Core only (see its CMakeLists).
set(DP_SRC "${CMAKE_CURRENT_SOURCE_DIR}/traj-compression/batch/DP/DP_adapted.cpp")
set(SQUISH_SRC "${CMAKE_CURRENT_SOURCE_DIR}/traj-compression/online/SQUISH/SQUISH_adapted.cpp")
set(DOTS_CMAKE "${CMAKE_CURRENT_SOURCE_DIR}/traj-compression/online/DOTS/CMakeLists.txt")
if(EXISTS "${DP_SRC}")
  add_executable(dp "${DP_SRC}")
endif()
if(EXISTS "${SQUISH_SRC}")
  add_executable(squish "${SQUISH_SRC}")
endif()
if(EXISTS "${DOTS_CMAKE}")
  add_subdirectory(traj-compression/online/DOTS)
endif()
