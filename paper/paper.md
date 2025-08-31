---
title: "ExaGOOP: an AMReX-based material point method solver"
tags:
  - C++
  - material point method
  - Exascale
  - heterogenous computing
authors:
  - name: Sreejith N. A.
    orcid: 0000-0001-5685-4070
    affiliation: 1
  - name: Nicholas Deak
    orcid: 0000-0003-1543-9334
    affiliation: 1
  - name: Yudong Li
    orcid: 0000-0002-7024-433X
    affiliation: 2
  - name: Hariswaran Sitaraman
    orcid: 0000-0001-5304-1664
    affiliation: 1
  - name: Marc Day
    orcid: 0000-0002-1711-3963
    affiliation: 1
affiliations:
  - name: Scalable Algorithms, Modeling and Simulation (SAMS) Group, National Renewable Energy Laboratory, USA
    index: 1
  - name: Catalytic Carbon Transformation and Scale-Up Center, National Renewable Energy Laboratory, USA
    index: 2  
date: 30 August 2025
bibliography: paper.bib
---
#Summary

ExaGOOP is a versatile continuum mechanics solver based on the material point method (MPM). Traditional numerical solvers typically employ methods such as the finite volume method, finite element method, or finite difference method. These approaches require the entire computational domain to be divided into a large number of computational units called cells, where the governing equations are solved in integral, variational, or differential forms. The necessity of a computational grid makes it challenging to solve problems involving highly deformable and history-dependent materials. In contrast, MPM, like many particle-based methods, is based on the Lagrangian framework of the governing equations. Unlike conventional solvers that rely on grids, MPM stores the material's properties on a collection of particles (also called as material points). While MPM does require a background grid, it is used only as a temporary tool for calculating gradients and can be reset after each time integration step. This approach effectively eliminates issues associated with grid element deformations. As a result, MPM methods are particularly well-suited for a wide range of continuum mechanics problems, especially those that involve significant material deformations.

ExaGOOP heavily depends on the AMReX library, which has been widely utilized in software applications [@PeleLMeX_JOSS, @PeleSoftware, @Sitaraman2021] employing adaptive mesh refinement (AMR) methods. The AMReX library  facilitates the generation of a block-structured, Cartesian background grid within ExaGOOP. Currently, the implementation supports only a single AMR level, but there are plans to extend ExaGOOP to incorporate higher AMR levels in the future. The material points related operations are managed by the particle classes provided by AMReX. The various steps in an MPM time update include particle-to-grid (P2G), nodal update, grid-to-particle (G2P), and position update. ExaGOOP offers users the flexibility to select the spatial discretization scheme, allowing for the use of linear-hat, quadratic B-spline, or cubic B-spline shape functions for both the P2G and G2P operations. Currently, the nodal update is performed using explicit Euler time integration; however, there are plans to incorporate implicit time stepping schemes in the future. In addition to these options, ExaGOOP allows users to select various numerical input parameters, such as the particle-in-cell (PIC)-Fluid Implicit Particle (FIP) blending factor in the G2P step, and whether to use Update Stress Last (USL) or Modified Update Stress Last (MUSL) for stress calculations. The software also supports CFL-based adaptive time-stepping. At present, ExaGOOP supports barotropic fluid and linear elastic solid constitutive models. However, adding new constitutive models is relatively straightforward for users, requiring only the development of the new constitutive model function without necessitating changes to other parts of the code. Complex, static wall boundaries are simulated using the level set method, while moving boundaries can be simulated with fictitious rigid material points.

ExaGOOP has undergone extensive validation and verification using 1D, 2D, and 3D test cases, all of which are available in the GitHub repository. Preprocessing scripts in the repository enable users to generate initial material point distribution with the desired number of material points per cell for either for user-defined simple geometries or based on user-provided images of complex bodies. Users can specify the constitutive model for each material point, facilitating multi-body and multi-phase simulations with ease.

ExaGOOP is developed in C++ and utilizes parallelization subroutines from the AMReX library. It employs an MPI+X approach, where MPI is used to distribute adaptive mesh refinement (AMR) grid patches across different CPU ranks. Each grid can be further divided into logical tiles, which can be distributed among threads using OpenMP on multi-core CPU machines or among GPU threads using CUDA, HIP, or SYCL on GPU-accelerated systems.

#Statement of Need
Other MPM solvers can be found online such as Karamelo [@devaucorbeil2021karamelo], Matter [@blatny2025matter], [@kumar2019geosmpm] and [@hu2018taichimpm]. The unique feature of ExaGOOP is its reliance on the AMReX framework for execution on CPU, GPU and heterogenous architectures. Exagoop’s integration with AMReX gives it access to advanced memory management, parallelism, and embedded boundary support. The demonstrated exascale performance of the AMReX framework on other solvers [cite] suggests the potential of ExaGOOP MPM solver to scale to billions of particles on GPU accelerated and heterogenous computing machines.

ExaGOOP is intended for students, researchers and engineers interested in simulating multi-material dynamics involving severe deformations. Initially developed as a tool to study membrane compaction in high-pressure reverse osmosis application [@nrel2023amrexmpm, @nrel2023exagoop], ExaGOOP is 


#Acknowledgements
The development of this software was supported by the National Alliance for Water Innovation (NAWI), funded by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy (EERE), Advanced Manufacturing Office, under Funding Opportunity Announcement Number DE-FOA-0001905. All of the research was performed using computational resources sponsored by the Department of Energy’s Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy Laboratory. This work was authored in part by the National Renewable Energy Laboratory, operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. The views expressed in the article do not
necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes.
