Key Features
=============================

The EXAGOOP solver suite is a Material Point Method (MPM) based multi-phase solver developed by the HPACF team at National Renewable Energy Laboratory, Colorado. The solver, which can be broadly defined as a particle-based method, can be used to solve gaseous, liquid, and solid phase physical problems. The Exagoop solver is developed based on the AMReX framework and has demonstrated excellent scale-up performance both on CPU and GPU cores.

Key features of ExaGOOP MPM solver include:

- Single-level, Cartesian, background grid generation and manipulation using AMReX library
- Material point functionalities implemented using AMReX particles
- CPU+GPU implementation on multiple heterogenous architectures
- Explicit time integration scheme
- Bilinear, Quadratic B-Spline, Cubic B-Spline shape functions
- Particle in Cell (PIC) & Fluid Implicit Particle (FLIP) methods
- Update stress last (USL) and Modified update stress last (MUSL) stress update schemes
- Multiple constitutive models implemented

