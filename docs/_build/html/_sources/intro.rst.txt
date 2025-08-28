.. role:: raw-latex(raw)
   :format: latex
..

.. _`sec:intro`:

Introduction
============

Disposal of concentrated brines is a significant contributor to
operating costs in industrial desalination and current methods to
achieving zero-liquid/minimum-liquid discharge (ZLD/MLD) is through
energy and carbon intense thermal desalination. An alternative energy
and cost efficient pathway is through Ultra-High Pressure Reverse
Osmosis (UHPRO) that can operate at high brine concentrations on the
order of 200-250 g/L. The high osmotic pressure requirement for UHPRO
leads to operating pressures on the order of 200 bar compared to
traditional salt water RO that requires :math:`\sim` 70 bar. RO
operation at very high pressure (:math:`\sim` 200 bar) is challenging
due to the degradation of membrane performance from 1) mechanical
compaction of Polysulfone (PSF) and polyester support layer thus
reducing porosity and permeability, 2) embossing of the membrane into
permeate spacers that impedes permeate flow and increases pressure
losses and 3) increased concentration polarization in viscous
hyper-saline brines that leads to reduced mass-transfer.

In order to address these challenges, an improved understanding of the
mechanisms of compaction and embossing in spiral wound elements at UHPRO
operating pressures is required. This report documents the development
of a computational model for membrane structural mechanics using the
material-point-method (MPM), using which membrane deformation under high
pressure is simulated. Our simulations agrees well with experimental
measurements on the overall displacement of the PSF layer and
qualitatively shows the reduction in porosity through macrovoid closure.
This report is organized as follows. Section `[] <#>`__ describes the
model equations and the MPM algorithm used to solve them. Section
`[] <#>`__ describes the verification of our solver with canonical solid
mechanics problems, Section `[] <#>`__ presents the validation of our
simulation methodology with experimental measurements of membrane
compaction under pressure.

Numerical simulations have become an indespensible part of engineering
analysis today. They are a good substitute to costly experiments and are
often used to downselect prime designs for performing experimental
tests. They are also desirable due to the availability of numerical
datasets that can shed great insights to spatial and temporal scales
which are hitherto inaccessible in experiments. Numerical simulations of
HPRO membrane compaction and flow in HPRO devices have been attempted in
the past
:raw-latex:`\cite{Gu2017,Pankaj2016,Abdelbaky2019,Liang2019,Lelong2019,Mao2021,Benjamin2022,Aschmoneit2022}`.
Recently, numerical simulations are even used in computer animations to
be used in movies such as *Frozen, Big Hero 6* and *Zootopia*.

The conventional and the oldest numerical methods applied to continuum
mechanics are based on either Finite Volume Methods(FVMs), Finite
Element Methods (FEMs) or Finite Difference Methods (FDMs). While the
first two methods are based on the integral form of the governing
quations, the last one is based on the differential form by
approximating the derivatives using algebraic expressions. An aspect
which is common to the above three methods is the need for a
computational grid to solve the governing equations. The grid refers to
a collection of numbered points (or nodes) which are connected using
edges (in 2D) or faces (in 3D) to form the full computational domain. It
is common in solid mechanics to use the above class of methods in
conjunction with Lagrangian form of the governing equations. In such
applications, the grid nodes are attached to the material and move as
the material deforms. When the material undergoes severe deformation,
the mesh nodes move along with the material and lead to grid
entanglement. This ultimately leads to numerical instabilities and the
solver stopping abrptly. On the other hand, Eulerian description based
methods do not require grid movement as the material advects across grid
cells and are not attached to grid nodes. Hene, this formulation
although numerically stable, is not a good candidate for problems with
moving interfaces and with materials with history dependant properties.

Material point method (MPM)
:raw-latex:`\cite{SULSKY1994179, Bardenhagen2004, Zhang2015TheMP, Vaucorbeil}`
is a ’mesh-free’ method that has received a lot of attention recently
for its ability to simulate problems with severe deformations. MPM is a
Lagrangian, particle-based numerical method developed by drawing
inspiration from the particle-in-cell (PIC)
:raw-latex:`\cite{osti_4769185}` and fluid implicit particle method
(FLIP) :raw-latex:`\cite{BRACKBILL1986314}` methods. Although MPMs were
initially applied to study solid mechanics problems, it has been
considerably extended to be applied to fluid simulations as well. In MPM
method, the full material domain is descretised using particles or
material points. All the material properties such as the velocity,
density, strainrates and stresses are stored at the material points.
This makes this method attractive to model history dependant
constitutive models. The lack of a grid connecting these material points
make MPM suitable to simulate problems with severe material deformations
such as in solid fracture, foam deformations and in granular flows.
Although a background grid is necessary in MPM, it is simply used as a
scratch pad to carry out certain numerical operations that are not
computationally intensive. These characteristics of MPM make it a good
candidate for studying membrane compaction problems.
