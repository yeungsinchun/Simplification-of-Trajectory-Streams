cmake_minimum_required(VERSION 3.16)
project(simplification_of_trajectory_stream CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
endif()

# Enable Qt auto-moc/resources
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC OFF)

find_package(CGAL REQUIRED COMPONENTS Qt6)
find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets OpenGL OpenGLWidgets)
find_package(Boost REQUIRED)

add_definitions(-DCGAL_USE_BASIC_VIEWER)

# simplify target (multi-source)
add_executable(simplify
  simplify.cpp
  drawing.cpp
)

target_link_libraries(simplify
  PRIVATE
    CGAL::CGAL_Qt6           # brings CGAL core + viewer
    Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets
)

target_compile_options(simplify PRIVATE -Wall -Wextra -Wpedantic)

if (CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT MSVC)
  target_compile_options(simplify PRIVATE -fsanitize=address,undefined -g)
  target_link_options(simplify PRIVATE -fsanitize=address,undefined -g)
endif()

# generate target (single source; you can keep macro or make it explicit)
create_single_source_cgal_program("generate.cpp")
target_link_libraries(generate PRIVATE CGAL::CGAL_Qt6 Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets)

#[===[
add_executable(test
  test.cpp
  drawing.cpp
)

target_link_libraries(test
  PRIVATE
    CGAL::CGAL_Qt6           # brings CGAL core + viewer
    Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets
)

add_executable(pikachu
  pikachu.cpp
  drawing.cpp
)

target_link_libraries(pikachu
  PRIVATE
    CGAL::CGAL_Qt6           # brings CGAL core + viewer
    Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets
)
#]===]