cmake_minimum_required(VERSION 3.16)
project(simplification_of_trajectory_stream CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC OFF)
set(CMAKE_AUTOUIC OFF)

find_package(CGAL REQUIRED COMPONENTS Qt6)

option(SIGNPOSTS "Emit os_signpost points-of-interest for Instruments profiling" OFF)
if(SIGNPOSTS)
    add_compile_definitions(SIGNPOSTS)
endif()

add_executable(simplify
  simplify.cpp
)

target_link_libraries(simplify
  CGAL::CGAL
)

add_executable(simplify_with_gui
  simplify_with_gui.cpp
  drawing.cpp
)

target_link_libraries(simplify_with_gui
  CGAL::CGAL_Qt6
)

target_compile_definitions(simplify_with_gui PRIVATE CGAL_USE_BASIC_VIEWER)

# Normalize raw taxi data
add_executable(normalize
  scripts/normalize.cpp
)
