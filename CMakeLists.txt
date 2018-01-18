cmake_minimum_required(VERSION 3.0)
project(nonlinear_systems)
set(CMAKE_BUILD_TYPE Release)

# if g++ show all warnings
if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-w")
endif(CMAKE_COMPILER_IS_GNUCXX)

# version string 0.1
set(MY_LIBRARY_VERSION_MAJOR 0 CACHE STRING "major version" FORCE)
set(MY_LIBRARY_VERSION_MINOR 1 CACHE STRING "minor version" FORCE)
set(MY_LIBRARY_VERSION ${MY_LIBRARY_VERSION_MAJOR}.${MY_LIBRARY_VERSION_MINOR} CACHE STRING "version" FORCE)

# use c++11 standard
set(CMAKE_CXX_STANDARD 11)

# include boost
find_package(Boost 1.61.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(test)

