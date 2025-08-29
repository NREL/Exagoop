.. PeleLMeX documentation master file, created by
   sphinx-quickstart on Fri Mar 18 15:03:51 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ExaGOOP's documentation!
====================================
Numerical simulations have become an indispensable part of engineering analysis today. They serve as an effective alternative to costly experiments and are often used to selectively identify optimal designs for experimental testing. The conventional and oldest numerical methods applied in continuum mechanics include Finite Volume Methods (FVMs), Finite Element Methods (FEMs), and Finite Difference Methods (FDMs). The first two methods rely on the integral form of the governing equations, while the last method is based on the differential form, approximating derivatives using algebraic expressions. A common aspect of these three methods is the necessity of a computational grid to solve the governing equations. This grid consists of a collection of numbered points (or nodes) that are interconnected by edges (in 2D) or faces (in 3D) to form the entire computational domain. In solid mechanics, these methods are typically used in conjunction with the Lagrangian form of the governing equations. In this context, the grid nodes are attached to the material and move as the material deforms. However, when the material experiences severe deformation, the movement of the mesh nodes can lead to grid entanglement, resulting in numerical instabilities and causing the solver to stop abruptly. In contrast, Eulerian-based methods do not require grid movement; instead, the material moves across fixed grid cells, decoupling it from the grid nodes. While this formulation is numerically stable, it is less suitable for problems involving moving interfaces and materials with history-dependent properties.

The Material Point Method (MPM) is a ‘mesh-free’ technique that has gained considerable attention for its ability to simulate problems involving significant deformations. MPM is a Lagrangian, particle-based numerical method that draws inspiration from the Particle-In-Cell (PIC) and Fluid-Implicit Particle (FLIP) methods. Although MPM was initially developed for solid mechanics problems, it has since been expanded for use in fluid simulations as well. In the MPM, the entire material domain is discretized using particles or material points, where all material properties—such as velocity, density, strain rates, and stresses—are stored. This aspect of MPM makes it particularly useful for modeling history-dependent constitutive models. The absence of a grid connecting these material points allows MPM to effectively simulate scenarios with severe material deformations, such as solid fractures, foam deformations, and granular flows. Although a background grid is necessary for MPM, it is merely used as a scratch pad to perform certain numerical operations that are not computationally intensive. These unique features make MPM a promising candidate for studying membrane compaction problems.

.. toctree::
   :maxdepth: 2
   :caption: The ExaGOOP MPM Solver:
   
   intro.rst

.. toctree::
   :maxdepth: 2
   :caption: Theory:
   
   theory.rst
   mpm_equations.rst
   
.. toctree::
   :maxdepth: 2
   :caption: Validation and Verification:

   

.. toctree::
   :maxdepth: 2
   :caption: Tutorials:

   Tutorials.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
