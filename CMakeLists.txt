cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(ExaGOOP LANGUAGES CXX)

#Setting CMAKE basics
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Setting ExaGOOP permanent settings. Do NOT change them
set(EXAGOOP_DIM "3" CACHE STRING "Number of spatial dimensions")
set(EXAGOOP_PRECISION DOUBLE)
set(EXAGOOP_ENABLE_EB ON)
set(EXAGOOP_ENABLE_PARTICLES ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

# HPC options
option(EXAGOOP_ENABLE_MPI "Enable MPI" ON)
option(EXAGOOP_ENABLE_OPENMP "Enable OpenMP" ON)
option(EXAGOOP_ENABLE_CUDA "Enable CUDA" OFF)
option(EXAGOOP_ENABLE_HIP "Enable HIP" OFF)
option(EXAGOOP_ENABLE_SYCL "Enable SyCL" OFF)

# Setting CUDA
if(EXAGOOP_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 10.0 or newer.")
  endif()
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 90)
  endif()
endif()

# Setting HIP
if(EXAGOOP_ENABLE_HIP)
   enable_language(HIP)
endif ()

# Setting MPI Flags
if(EXAGOOP_ENABLE_MPI)
  message(STATUS "MPI Configure Section")
  find_package(MPI REQUIRED CXX)
endif()

# Setting AMReX Configuration
message(STATUS "AMReX Configure Section")
set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/amrex")
include(SetAmrexOptions)
list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
add_subdirectory(${AMREX_SUBMOD_LOCATION})
include(SetAmrexCompileFlags)

# General information about machine, compiler, and build type
message(STATUS "${PROJECT_NAME} Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "EXAGOOP_PRECISION = ${EXAGOOP_PRECISION}")

include(BuildExeAndLib)
