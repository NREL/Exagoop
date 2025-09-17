cmake_minimum_required(VERSION 3.14)
project(ExaGOOP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# HPC options
option(EXAGOOP_ENABLE_MPI "Enable MPI" ON)
option(EXAGOOP_ENABLE_OPENMP "Enable OpenMP" OFF)
option(EXAGOOP_ENABLE_CUDA "Enable CUDA" OFF)
option(EXAGOOP_ENABLE_HIP "Enable HIP" OFF)
option(EXAGOOP_ENABLE_SYCL "Enable SyCL" OFF)

# AMREX flags
set(AMReX_EB ON)
set(AMReX_EB2 ON)
set(AMReX_PARTICLES ON)
set(AMReX_SPACEDIM 3)
set(AMReX_FORTRAN OFF)
set(AMReX_SUNDIALS OFF)
set(AMReX_LINEAR_SOLVERS OFF)
set(AMReX_OMP ${EXAGOOP_ENABLE_OPENMP})

#cmake_policy(SET CMP0104 NEW)

#MPI Flags
if(EXAGOOP_ENABLE_MPI)
  message(STATUS "MPI Configure Section")
  find_package(MPI REQUIRED CXX)
endif()

#CUDA Flags
if(EXAGOOP_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 10.0 or newer.")
  endif()
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 86 75 70)
  endif()
endif()

if(EXAGOOP_ENABLE_HIP)
   enable_language(HIP)
endif ()

# Fetch AMReX
include(FetchContent)
FetchContent_Declare(
  amrex
  GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git
  GIT_TAG        development
)
FetchContent_MakeAvailable(amrex)

# Include your MPM source
include_directories(${PROJECT_SOURCE_DIR}/Source)

# Collect source files
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/Source/*.cpp)

# Create executable
add_executable(ExaGOOP ${SOURCES})
target_link_libraries(ExaGOOP PRIVATE amrex)