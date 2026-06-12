cmake_minimum_required(VERSION 3.16)
project(simplification_of_trajectory_stream CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Debug CACHE STRING "BUILD TYPE" FORCE)
# endif()
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer -g)
    add_link_options(-fsanitize=address)
endif()

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC OFF)

find_package(CGAL REQUIRED COMPONENTS Qt6)

find_package(OpenCV COMPONENTS core imgproc PATHS /opt/homebrew NO_DEFAULT_PATH)
find_package(Boost REQUIRED)

add_definitions(-DCGAL_USE_BASIC_VIEWER)

add_executable(simplify
  simplify.cpp
  drawing.cpp
)

add_executable(simplify_old
  simplify_old.cpp
  drawing_old.cpp
)

add_executable(simplify_opencv
  simplify_opencv.cpp
  drawing.cpp
)

add_executable(simplify_boost
  simplify_boost.cpp
  drawing.cpp
)

add_executable(simplify_SH
  simplify_SH.cpp
  drawing_old.cpp
)

target_include_directories(simplify_boost PRIVATE ${OpenCV_INCLUDE_DIRS})
target_include_directories(simplify_opencv PRIVATE ${OpenCV_INCLUDE_DIRS})

target_link_libraries(simplify
  CGAL::CGAL_Qt6
)

target_link_libraries(simplify_old
  CGAL::CGAL_Qt6
)

target_link_libraries(simplify_opencv
  CGAL::CGAL_Qt6
  ${OpenCV_LIBS}
)

target_link_libraries(simplify_boost
  CGAL::CGAL_Qt6
  Boost::boost
)

target_link_libraries(simplify_SH
  CGAL::CGAL_Qt6
)

# Viewer for plotting two curves from files
add_executable(plot_curves
  plot_curves.cpp
  drawing.cpp
)

target_link_libraries(plot_curves
  CGAL::CGAL_Qt6
)

# Algorithms runner (single build tree)
add_executable(main
  algorithms/main.cpp
)

target_include_directories(main PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/algorithms/inc
)

# Data normalization tool (C++)
add_executable(normalize
  tools/normalize.cpp
)