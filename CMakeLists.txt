# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)

project( Mobile )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()

endif()

# QT6
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

add_definitions(-DCGAL_USE_BASIC_VIEWER)


# include for local directory

# include for local package


# Creating entries for all C++ files with "main" routine
# ##########################################################


create_single_source_cgal_program( "test.cpp" )

target_link_libraries(test PUBLIC CGAL::CGAL_Qt6)