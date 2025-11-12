cmake_minimum_required(VERSION 3.16)
project(simplification_of_trajectory_stream CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "BUILD TYPE" FORCE)
endif()

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC OFF)

find_package(CGAL REQUIRED COMPONENTS Qt6)

add_definitions(-DCGAL_USE_BASIC_VIEWER)

add_executable(simplify
  simplify.cpp
  drawing.cpp
)

target_link_libraries(simplify
  PRIVATE
    CGAL::CGAL_Qt6
    Qt6::Core Qt6::Gui Qt6::Widgets
)

# Viewer for plotting two curves from files
add_executable(plot_curves
  plot_curves.cpp
  drawing.cpp
)

target_link_libraries(plot_curves
  PRIVATE
    CGAL::CGAL_Qt6
    Qt6::Core Qt6::Gui Qt6::Widgets
)

#[===[
create_single_source_cgal_program("generate.cpp")
target_link_libraries(generate PRIVATE CGAL::CGAL_Qt6 Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets)

#[===[
add_executable(test

#[===[ pp
  drawing.cpp
)

target_link_libraries(test
  PRIVATE
    CGAL::CGAL_Qt6           # brings CGAL core + viewer
    Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets
)

add_executable(pikachu
  pi
achu.cpp
  drawing.cpp
)

target_link_libraries(pikachu
  PRIVATE
    CGAL::CGAL_Qt6           # brings CGAL core + viewer
    Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets
)
#]===]