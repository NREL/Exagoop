# EXAGOOP - A material point method (MPM) solver based on AMReX Framework
## Introduction

The EXAGOOP solver suite is a Material Point Method (MPM) based multi-phase solver developed by the HPACF team at National Renewable Energy Laboratory, Colorado. The solver, which can be broadly defined as a particle-based method, can be used to solve gaseous, liquid, and solid phase continuum, physical problems.

Most of the fluid and solid dynamics solvers are developed based on the Eulerian framework of the governing equations. These methods solve the governing equations on a collection of structured or unstructured grid elements. However, when materials undergo large deformations, the grid elements also stretch and deform, leading to inaccurate and often unstable computations. On the other hand, MPMs, like many particle-based methods, is based on the Lagrangian framework of the governing equations. The continuum material under study is modeled as a collection of particles (or material points). Unlike the Eulerian solvers, the material's properties are stored on these material points. The Lagrangian form of the governing equations are solved using the material points and a reference background grid. This eliminates the problem posed by grid element deformations. Hence, MPM methods are well suited for all continuum mechanics problems in general and for those involving large material deformations, in particular.

The basic components and the terminologies used in EXAGOOP MPM solver is shown in the figure below. A cartesian Eulerian grid (no adaptive mesh refinement) is used along with material points simulated as spherical particles.


The various steps involved in one time integration stage in EXAGOOP is shown in the following figure and are described in the following 4 steps.

Steps:

0. In the initialization stage, material point mass, position, velocity and stresses are initialized.
1. The material point (subscript p) mass and momentum are mapped onto the grid node (subscript I) using grid shape functions $\phi$. Similarly particle forces (external forces such as gravity and internal forces from stresses) are also mapped to grid nodes. Mathematically, 

$$
m_I^t = \sum_p \Phi_I (x_p^t) m_p
$$

$$
(m\mathbf{v})_I^t = \sum_p \Phi_I (x_p^t) (mv)_p
$$

$$
\mathbf{f}_I^{ext,t} = \sum_p \Phi_I (x_p^t) m_p \mathbf{b}(x_p)
$$

$$
\mathbf{f}_I^{int,t} = \sum_p V_p^t \mathbf{\Sigma}_p^t \nabla \Phi_I (x_p^t)
$$

where, $\sum_p$ denotes summation over all material points. $m_p$, $v$, $V_p$, $b$, $\Sigma$ and $\Phi$ denote particle mass, particle velocity vector, particle volume, body force, stress and grid shape functions respectively.

2. The updated grid nodal velocity is calculated using explicit Euler time integration. 

$$
\mathbf{v}_I^{t+\Delta t} = \mathbf{v}_I^{t} + \frac{\mathbf{f}_I^{ext,t}+\mathbf{f}_I^{int,t}}{m_I^t} 
$$

3. Particle velocity is updated with new nodal velocity at time $t+\Delta t$ using a blend of Particle in Cell (PIC) and Fluid Implicit Particle (FLIP) update. The blending factor used is $\alpha$.

$$
\mathbf{v}_p^{t+\Delta t}=\alpha\left(\mathbf{v}_p^t+\sum_I \Phi_I\left[\mathbf{v}_I^{t+\Delta t}-\mathbf{v}_I^t\right]\right)+\left(1-\alpha\right) \sum_I \Phi_I \mathbf{v}_I^{t+\Delta t} \nabla \mathbf{v}_p^{t+\Delta t}=\sum_I^{n g} \nabla \Phi_I \mathbf{v}_I^{t+\Delta t}
$$

4. Particle positions are updated using Euler integration.

$$
x_p^{\{t+\Delta t\}}=x_p^t+\Delta t \sum_I \phi_I\left(x_p^t\right) v_I^{\{t+\Delta t\}}
$$

## EXAGOOP features
## Build Instructions
## Run Instructions
## Visualization Instructions
