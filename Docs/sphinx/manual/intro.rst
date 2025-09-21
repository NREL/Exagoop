Key Features
=============================

EXAGOOP is a Material Point Method (MPM) based solver developed by the Scalable Algorithms, Modeling and Simulation (SAMS) team at the National Renewable Energy Laboratory, Colorado. The solver, which can be broadly defined as a particle-based method, can be used to solve physical problems involving gaseous, liquid, and solid phases. The Exagoop solver is developed based on the AMReX framework and has demonstrated excellent scalability performance on both CPU and GPU cores.

Key features of ExaGOOP MPM solver include:

- Single-level, Cartesian, background grid generation and manipulation using the AMReX library
- Material point functionalities implemented using AMReX particle library
- CPU, GPU and CPU+GPU implementation on multiple heterogeneous architectures
- Explicit time integration scheme
- Bilinear, Quadratic B-Spline, Cubic B-Spline shape functions
- Particle in Cell (PIC) & Fluid Implicit Particle (FLIP) methods
- Update stress last (USL) and Modified update stress last (MUSL) stress update schemes
- Multiple constitutive models for solids and fluids

