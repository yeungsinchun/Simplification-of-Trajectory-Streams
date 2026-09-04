cmake_minimum_required(VERSION 3.16)
project(simplification_of_trajectory_stream CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()

# Build the Qt viewer and the DOTS baseline used by the benchmark. Default
# ON for local development and the GitHub Actions runners. Set to OFF to
# produce a headless-only build with no Qt6 dependency (used by the Docker /
# Cloud Run image).
option(BUILD_GUI "Build the Qt viewer and DOTS baseline" ON)

# Headless simplify binary — always built, no Qt dependency.
find_package(CGAL REQUIRED)
add_executable(simplify
  simplify.cpp
)
target_link_libraries(simplify
  CGAL::CGAL
)

# Baseline binaries used by the web compare pane (no Qt).
add_executable(dp traj-compression/batch/DP/DP_adapted.cpp)
add_executable(squish traj-compression/online/SQUISH/SQUISH_adapted.cpp)

if(BUILD_GUI)
  # Qt-only builds need AUTOMOC for the Q_OBJECT types in the GUI target.
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTORCC OFF)
  set(CMAKE_AUTOUIC OFF)

  find_package(CGAL REQUIRED COMPONENTS Qt6)
  find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets OpenGL OpenGLWidgets Svg)

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
endif()
