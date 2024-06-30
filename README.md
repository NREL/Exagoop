# EXAGOOP - A material point method (MPM) solver based on AMReX Framework
## Introduction

The EXAGOOP solver suite is a Material Point Method (MPM) based multi-phase solver developed by the HPACF team at National Renewable Energy Laboratory, Colorado. The solver, which can be broadly defined as a particle-based method, can be used to solve gaseous, liquid, and solid phase continuum, physical problems.

Most of the fluid and solid dynamics solvers are developed based on the Eulerian framework of the governing equations. These methods solve the governing equations on a collection of structured or unstructured grid elements. However, when materials undergo large deformations, the grid elements also stretch and deform, leading to inaccurate and often unstable computations. On the other hand, MPMs, like many particle-based methods, is based on the Lagrangian framework of the governing equations. The continuum material under study is modeled as a collection of particles (or material points). Unlike the Eulerian solvers, the material's properties are stored on these material points. The Lagrangian form of the governing equations are solved using the material points and a reference background grid. This eliminates the problem posed by grid element deformations. Hence, MPM methods are well suited for all continuum mechanics problems in general and for those involving large material deformations, in particular.

The basic components and the terminologies used in EXAGOOP MPM solver is shown in the figure below. A cartesian Eulerian grid (no adaptive mesh refinement) is used along with material points simulated as spherical particles.


The various steps involved in one time integration stage in EXAGOOP is shown in the following figure and are described in the following 4 steps.

$$
\begin{itemize}
\end
$$

## EXAGOOP features
## Build Instructions
## Run Instructions
## Visualization Instructions
