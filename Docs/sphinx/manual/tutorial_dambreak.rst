.. highlight:: rst
	
This section presents a straightforward tutorial for using ExaGOOP to simulate the 2D dam break problem. While we have selected the dam break scenario for this example, the steps outlined here can be adapted to simulate any continuum mechanics problem.

To run a simulation with ExaGOOP, you will need an ExaGOOP input file and a particle (or material point) file, in addition to the ExaGOOP executable, which is built as described in the README file found on the `GitHub repository page <https://github.com/NREL/Exagoop>`_. Sample input and particle files can be found in the <ExaGOOP folder>/Tests/Dam_Break directory.

Setting up the simulation folder
--------------------------------

Create a folder (say 'Sim1_DamBreak') in your chosen location and navigate into it.

.. code-block:: bash

   # Create a directory 'Sim1_DamBreak' for dam break simulation
   mkdir Sim1_DamBreak
   cd Sim1_DamBreak
	
Copy the input and particle files from the Tests folder. Ensure the environment variable ``$MPM_HOME`` is set as detailed in the `README <https://github.com/NREL/Exagoop>`_  file.

.. code-block:: bash

   # Copy input and particle file from Tests/Dam_Break to current folder
   cp $MPM_HOME/Tests/Dam_Break/inputs_dambreak.in .
   cp $MPM_HOME/Tests/Dam_Break/mpm_particles.dat .
	
	
Copy the ExaGOOP executable from the build folder (depending on how it was built):

.. code-block:: bash

   # Copy ExaGOOP executable from Build_Cmake folder
   cp $MPM_HOME/Build_Cmake/ExaGOOP.exe .
	
or,

.. code-block:: bash

   # Copy ExaGOOP executable from Build_Gnumake folder
   cp $MPM_HOME/Build_Gnumake/ExaGOOP3d<suffix>.ex .
	
	
The <suffix> string is automatically decided based on the build environments.

Now that the simulation folder is prepared, we will discuss how to set up the input file in the next section.


Setting up the input file
-------------------------

As mentioned in the :ref:`mpmsect` section, MPM requires a background grid. ExaGOOP utilizes a single, Adaptive Mesh Refinement (AMR) level Cartesian background grid. The following section provides a detailed description of the necessary inputs required for setting up the background grid. Please note that all input parameter definitions in the input file should begin with the prefix ``mpm.`` to be read correctly by the solver.

.. code-block:: bash

   #Geometry parameters
   mpm.prob_lo = 0.0 0.0 0.0                       #Lower corner of physical domain
   mpm.prob_hi = 0.4 0.4 0.012                     #Upper corner of physical domain
   mpm.ncells  = 100 100 3                         #number of cells in each direction
   mpm.max_grid_size = 101                         #Max grid size
   mpm.is_it_periodic = 0  0  1                    #Periodicity
   

In the input file section above, ``prob_lo`` specifies the 3-dimensional coordinates of the lower left corner of the background grid. In this case, the x, y, and z coordinates are all set to 0.0. Similarly, ``prob_hi`` indicates the upper right corner of the computational domain (0.4, 0.4 and 0.012). It is important to note that, although the problem is inherently 2-dimensional, ExaGOOP models it as a 3-dimensional problem by including cells in the z-direction.

With the background grid domain fully defined, the Cartesian grid parameters are set using ``mpm.ncells``, which specifies the number of cells in each of the three directions. Since we will be using a linear hat shape function, it is recommended to set the number of cells in the z-direction to 3. The ``max_grid_size`` parameter defines the maximum number of cells in an AMReX box and can be used to adjust the number of boxes that are passed to each MPI rank in a parallel simulation. Additionally, because this is a 2-dimensional problem, the z-direction is simulated using a translational periodic boundary by setting ``mpm.is_it_periodic`` to ``0 0 1``, which indicates that the z-boundary is periodic (indicated by 1) while the other two boundaries are not (indicated by 0).

The next step is to define how the material points, also referred to as particles, will be specified. ExaGOOP offers two options for specifying material points. The first option is 'autogen', where ExaGOOP generates the material points internally. The second option allows users to specify their own material point file generated externally. If the user chooses to have ExaGOOP generate the material points during runtime, set the autogen flag to 1 by including ``mpm.use_autogen=1`` in the input file. After this, the user will need to specify the locations of the material points and their constitutive properties as outlined below.


.. code-block:: bash

   mpm.use_autogen=1                               #Use particle autogeneration tool
   mpm.mincoords_autogen=0.0 0.0 0.0
   mpm.maxcoords_autogen=1.0 1.0 1.0
   mpm.vel_autogen=0.0 0.0 0.0                     #Velocity components of particle
   mpm.constmodel_autogen=0                        #0->Elastic solid,1->Compressible fluid
   mpm.dens_autogen=1.0                            #Density
   mpm.E_autogen=1e6                               #Youngs modulous
   mpm.nu_autogen=0.3                              #Poisson's ratio
   mpm.bulkmod_autogen=2e6                         #Bulk modulous
   mpm.Gama_pres_autogen=7                         #Gamma
   mpm.visc_autogen=0.001                          #Viscosity
   mpm.multi_part_per_cell_autogen=1               #Number of particles per cell

In this context, ``mincoords_autogen`` and ``maxcoords_autogen`` define the lower and upper corners of a rectangular sub-domain within the background grid where the material points will be embedded. The variable ``vel_autogen`` represents the initial velocity components of the material points, while ``constmodel_autogen`` indicates the constitutive model used for these material points.

Alternatively, if the user opts to use an already existing particle file, the input statement ``mpm.particle_file=<particle filename>`` can be utilized to indicate the material point file.

Next, the numerical parameters for the simulation are defined. ExaGOOP offers three options for specifying the nodal shape function: 1 for linear hat functions, 2 for quadratic B-splines, and 3 for cubic B-splines. These options are set using the command ``mpm.order_scheme=<shape function flag>``. The PIC-FLIP blending factor, defined previously as :math:`\alpha_{P-F}`, can be specified with ``mpm.alpha_pic_flip``, while the stress update scheme (either MUSL or USL) can be indicated using the variable ``mpm.stress_update_scheme=1``  (for MUSL) or ``mpm.stress_update_scheme=1`` for USL. ExaGOOP employs an Explicit Euler time integration scheme. The time step size can either be fixed, indicated by ``mpm.fixed_timestep = 1``, or adaptive, with ``mpm.fixed_timestep = 0`` followed by the specification of the CFL number using ``mpm.CFL=<CFL number>``.

To prevent simulations from running into numerical instabilities, the parameters ``mpm.dt_min_limit`` and ``mpm.dt_max_limit`` can be set to constrain the time step values. The complete numerical setup block in the input file will appear as shown below:

.. code-block:: bash

   #Timestepping parameters
   mpm.fixed_timestep = 0                          #1=> fixed time step, 0=> CFL based adaptive time step
   mpm.timestep = 1.0e-5                           #Timestep value in case mpm.fixed_timestep=1
   mpm.CFL=0.1                                     #CFL number
   mpm.dt_min_limit=1e-12                          #The minimum value of timestep. Throw a warning if lower. Continues iteration with dt_min_limit
   mpm.dt_max_limit=1e-0                           #The maximum value of timestep. if(dt>dt_max_limit) dt=dt_max_limit
   
   #Numerical schemes
   mpm.order_scheme=1                              #1->Linear shape func. 3->Spline shape function
   mpm.alpha_pic_flip = 0.99                       #1.0->Pure FLIP scheme. 0.0->Pure PIC scheme. Recommended:(0.95-0.99)
   mpm.stress_update_scheme= 1                     #1->MUSL,0->USL
   mpm.mass_tolerance = 1e-18                      #Stability parameters to avoid zero mass on nodes. Recommended: 1e-18
   mpm.mpm.calculate_strain_based_on_delta=0       #Stress calculation based on delta epsilon if this flag =1
   
   
Problem-specific parameters are set next. ``mpm.final_time`` and ``mpm.max_steps`` specify the maximum duration of the simulation and the maximum number of iterations, with the simulation stopping when either condition is met first. ``screen_output_time`` indicates how often (in flow-time) the simulation log and diagnostics will be output, where as, ``write_output_time`` sets the frequency for writing output files.
 
.. code-block:: bash
	
	mpm.final_time=0.3                      #Maximum simulation time
	mpm.max_steps=10                                #Maximum number of iterations for the simulation
	mpm.screen_output_time = 0.001                  #How frequently to output iteration msgs
	mpm.write_output_time=0.01                      #How frequently to write output files
	mpm.num_redist = 1                              #How frequently to redistribute
	
Finally, the problems specific boundary conditions are specified using the Boundary conditions block as shown below,

.. code-block:: bash
	
	#Boundary conditions
	mpm.bc_lower= 2 2 0                      #0->Periodic 1->Noslipwall 2->Slipwall 3->Outflow
	mpm.bc_upper= 2 2 0                      #0->Periodic 1->Noslipwall 2->Slipwall 3->Outflow
	
``bc_lower`` and ``bc_upper`` refer to the lower and higher three faces of the computational domain. The boundary flags are defined as follows: 0 for periodic, 1 for no-slip wall, 2 for slip wall, and 3 for outflow.
	
Generating initial material point file
--------------------------------------
When ``mpm.use_autogen`` is set to 0, the user must specify the initial material point file. ExaGOOP expects this file to be in ASCII format. For a detailed understanding of the required format, the user may refer to the example script found at ``./Tests/Dam_Break/generate.py``.
 
Running ExaGOOP and viewing output files
----------------------------------------
Once the input file and material point files are correctly set up, the ExaGOOP executable can be run as:

.. code-block:: bash
	
	./ExaGOOP input.in 

when run in serial mode or,

.. code-block:: bash
	
	mpirun -n <nproc> ./ExaGOOP input.in 
	
or when run using MPI with ``<nproc>`` being the number of MPI ranks. ExaGOOP generates two sets of output files. The first is a particle file, which can be named by setting the option ``mpm.prefix_particlefilename=<output-particle-filename>``. The second is a grid file, with its name specified using ``mpm.prefix_gridfilename=<output-grid-filename>`` in the input file. Both files can be opened and visualized using ParaView, utilizing the AMReX grid and particle file readers.