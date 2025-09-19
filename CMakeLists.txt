cmake_minimum_required(VERSION 3.23 FATAL_ERROR)
project(ExaGOOP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(EXAGOOP_DIM "3" CACHE STRING "Number of physical dimensions")
set(EXAGOOP_PRECISION DOUBLE)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

# HPC options
option(EXAGOOP_ENABLE_MPI "Enable MPI" ON)
option(EXAGOOP_ENABLE_OPENMP "Enable OpenMP" OFF)
option(EXAGOOP_ENABLE_CUDA "Enable CUDA" ON)
option(EXAGOOP_ENABLE_HIP "Enable HIP" OFF)
option(EXAGOOP_ENABLE_SYCL "Enable SyCL" OFF)

if(EXAGOOP_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 10.0 or newer.")
  endif()
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 86)
  endif()
endif()

if(EXAGOOP_ENABLE_HIP)
   enable_language(HIP)
endif ()

#MPI Flags
if(EXAGOOP_ENABLE_MPI)
  message(STATUS "MPI Configure Section")
  find_package(MPI REQUIRED CXX)
endif()

message(STATUS "AMReX Configure Section")
set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/amrex")
include(SetAmrexOptions)
list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
add_subdirectory(${AMREX_SUBMOD_LOCATION})
include(SetAmrexCompileFlags)

# Include your MPM source
include_directories(${PROJECT_SOURCE_DIR}/Source)

# Collect source files
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/Source/*.cpp)

# Create executable
add_executable(ExaGOOP ${SOURCES})

# General information about machine, compiler, and build type
message(STATUS "${PROJECT_NAME} Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "PELE_PRECISION = ${PELE_PRECISION}")

#if(DEFINED ENV{CONDA_PREFIX})
#    link_directories($ENV{CONDA_PREFIX}/lib)
#    target_link_libraries(ExaGOOP PRIVATE amrex $<$<BOOL:${EXAGOOP_ENABLE_CUDA}>:CUDAToolkit::cudart)
#endif()
target_link_libraries(ExaGOOP PRIVATE amrex)
