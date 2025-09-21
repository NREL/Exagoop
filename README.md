# EXAGOOP - A material point method (MPM) solver based on AMReX Framework
## Overview

![MPM_Github](https://github.com/NREL/Exagoop/assets/98907926/59a06fe9-56c3-4822-9d57-996e7beebd0f)

ExaGOOP is a multi-phase solver based on the Material Point Method (MPM), developed by the Scalable Algorithms, Modeling, and Simulation (SAMS) team at the National Renewable Energy Laboratory in Golden, Colorado. This solver, which can be classified as a particle-based method, is designed to tackle continuum physical problems involving gaseous, liquid, and solid phases, particularly those that experience large deformations.

Key features of ExaGOOP MPM solver include:

- Built on the AMReX framework, utilizing a single-level, block-structured Cartesian grid as the background grid. Material points are represented using the particle class in AMReX.
- Support for MPI and GPU computing through CUDA and HIP.
- Availability of linearly elastic and compressible fluid material models.
- Support for linear, quadratic, and cubic B-Spline grid shape functions.
- Use of explicit time integration for time advancement.
- Simulation of complex geometries using the embedded boundary method.
- Rigid material points are available for simulating either stationary or moving rigid walls.

For more details on the governing equations, numerical methods, validation tests, and tutorials, please visit: [ExaGOOP Documentation](https://nrel.github.io/Exagoop/).

## Installing ExaGOOP

### Requirements:
ExaGOOP is designed to run on macOS and Linux environments and can be executed on various computing devices, including laptops, workstations, and high-performance computing (HPC) systems. Although we do not provide official support for running ExaGOOP on Windows, we recommend that Windows users utilize the Windows Subsystem for Linux (WSL) to build and run ExaGOOP.

ExaGOOP requires a C++17 compatible compiler, specifically GCC version 8 or higher and Clang version 3.6 or higher; Microsoft Visual C++ (MSVC) is not supported. Additionally, it requires CMake version 3.20 or higher, or GNU Make version 3.81 or higher to manage the build process. 

ExaGOOP is developed using the AMReX framework, which provides essential APIs for handling grid and particle functionalities. The built-in capabilities of AMReX enable ExaGOOP to be built on heterogeneous computing architectures that utilize MPI+X, where X refers to GPU architectures such as NVIDIA (CUDA version 11 or higher) and AMD (ROCm version 5.2 or higher) GPUs.

### Downloading ExaGOOP

To clone the ExaGOOP source code to your local environment, you can directly access the GitHub repository at this link: https://github.com/NREL/Exagoop. The AMReX framework is included as a submodule within the ExaGOOP repository. Use the following command to clone the complete source code:

```
git clone --recurse-submodules https://github.com/NREL/Exagoop.git
```

This command will create a folder named `Exagoop` containing the full ExaGOOP source code. The AMReX sources can be found in the `Exagoop/Submodules/amrex` directory.

Next, set up the necessary environment variables for building ExaGOOP. You can do this by typing the following commands in your terminal. If you'd like these settings to persist, consider adding them to your `.bashrc` or `.zshrc` file.

```
export MPM_HOME=<path_to_ExaGOOP>
export AMREX_HOME=${MPM_HOME}/Submodules/amrex
```
### Building ExaGOOP
Once the above steps are completed, you can build ExaGOOP using either CMake or GNU Make.
#### Building using CMake
1. Navigate to the `Build_Cmake` directory within `$MPM_HOME`:

   ```bash
   cd $MPM_HOME/Build_Cmake
   ```
2. Edit the `cmake.sh` file to configure the build environment. Be sure to set the following options as needed:
   - Set `-DEXAGOOP_ENABLE_MPI:BOOL` to `ON` if you are using MPI.
   - Set `-DEXAGOOP_ENABLE_CUDA:BOOL` and `-DEXAGOOP_ENABLE_HIP:BOOL` to `ON` or `OFF` based on your GPU usage.
   - Define `-DAMReX_CUDA_ARCH` and `-DAMReX_AMD_ARCH` according to your GPU architecture. 
   - Always set `-DAMReX_SPACEDIM` to `3`, as ExaGOOP operates in three spatial dimensions.
 
3. Build the project by executing:

   ```bash
   sh cmake.sh
   ```
Upon successful completion, the ExaGOOP executable named `ExaGOOP.exe` will be available in the `$MPM_HOME/Build_Cmake` directory.

#### Building using GNUmake
1. Navigate to the `Build_Gnumake` directory within `MPM_HOME`:

   ```bash
   cd $MPM_HOME/Build_Gnumake
   ```

2. Adjust the `Gnumake` file according to your build environment:
   - Set the C++ compiler variable `COMP` to either `gnu` or `Clang`.
   - If you are using MPI, set `USE_MPI = TRUE`.
   - Similarly, set `USE_OMP`, `USE_CUDA`, and `USE_HIP` to `TRUE` if you plan to use OpenMP, CUDA, or HIP.

3. Compile and link the project by running:

   ```bash
   make -j<num_proc>
   ```

   Replace `<num_proc>` with the number of processors you want to utilize for parallel compilation. After a successful build, the executable will be located in the `$MPM_HOME/Build_Gnumake` directory.



### Visualization Instructions

- The simulation output files are in the form of AMReX plotfiles.
- Paraview can be used to load and view the particle ( __plt__ files) and nodal ( __nplt__ files) solution files.

## Getting help

- To engage with the development team, please use the [GitHub discussions link](https://github.com/NREL/Exagoop/discussions).
- If you encounter a bug that you would like to report, kindly use the [GitHub issues page](https://github.com/NREL/Exagoop/issues). When reporting, please provide as much detail as possible regarding the major compilation and runtime options you are using.

## Acknowledgment

The development of this software was supported by the National Alliance for Water Innovation (NAWI), funded by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy (EERE), Advanced Manufacturing Office, under Funding Opportunity Announcement Number DE-FOA-0001905. All of the research was performed using computational resources sponsored by the Department of Energy’s Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy Laboratory. This work was authored in part by the National Renewable Energy Laboratory, operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes.
