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
find_package(OpenCV COMPONENTS core imgproc PATHS /opt/homebrew NO_DEFAULT_PATH)
# find_package(Boost REQUIRED)

add_executable(simplify
  simplify.cpp
  drawing.cpp
)

target_include_directories(simplify PRIVATE ${OpenCV_INCLUDE_DIRS})

target_link_libraries(simplify
  CGAL::CGAL_Qt6
  ${OpenCV_LIBS}
)

# Normalize raw taxi data
add_executable(normalize
  tools/normalize.cpp
)

# Plot all curves in data/<id>/ in one window, clickable to toggle.
add_executable(plot_curve
  tools/plot_curve.cpp
  drawing.cpp
)

target_include_directories(plot_curve PRIVATE ${CMAKE_SOURCE_DIR})

target_compile_definitions(simplify PRIVATE CGAL_USE_BASIC_VIEWER)
target_compile_definitions(plot_curve PRIVATE CGAL_USE_BASIC_VIEWER)

target_link_libraries(plot_curve
  CGAL::CGAL_Qt6
)

# Compare O'Rourke's convex-convex intersection with the CGAL implementation
# used in simplify.cpp.  Headless; no Qt, no OpenCV.
add_executable(test_convex_intersect test.cpp)
target_link_libraries(test_convex_intersect CGAL::CGAL)