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
find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets OpenGL OpenGLWidgets Svg)

add_executable(simplify
  simplify.cpp
)

target_link_libraries(simplify
  CGAL::CGAL
)

# Instrumented twin of `simplify`: same algorithm with TIMER hooks and a
# per-run timing breakdown printed to stderr.  Kept as a separate binary so the
# production `simplify` hot path carries zero timing overhead.
add_executable(simplify_with_time
  simplify_with_time.cpp
)

target_link_libraries(simplify_with_time
  CGAL::CGAL
)

add_executable(simplify_with_gui
  simplify_with_gui.cpp
  drawing.cpp
)

target_link_libraries(simplify_with_gui
  CGAL::CGAL_Qt6
  Qt6::Core
  Qt6::Gui
  Qt6::Widgets
  Qt6::OpenGL
  Qt6::OpenGLWidgets
  Qt6::Svg
)

target_compile_definitions(simplify_with_gui PRIVATE CGAL_USE_BASIC_VIEWER)

add_subdirectory(traj-compression/online/DOTS)
