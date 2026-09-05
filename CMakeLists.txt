cmake_minimum_required(VERSION 3.16)
project(simplification_of_trajectory_stream CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()

# Optional Qt viewer (Widgets / OpenGL). Default ON for local development.
# Set to OFF for headless builds (Docker / Cloud Run). DOTS only needs Qt6
# Core and is built with the other baselines below, not gated on this flag.
option(BUILD_GUI "Build the Qt viewer (simplify_with_gui)" ON)

# Headless simplify binary — always built, no Qt dependency.
find_package(CGAL REQUIRED)
add_executable(simplify
  simplify.cpp
)
target_link_libraries(simplify
  CGAL::CGAL
)

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
endif()
