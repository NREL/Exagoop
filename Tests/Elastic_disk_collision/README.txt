To run this test:
1. Open generate_elastic_collision_2d.py
  1.a) fill in appropriate values for:
    1.a.1) no_of_cells_in_z--> How many cells should be present in the z-direction
    1.a.2) ncells_x,ncells_y--> No of cells in x,y
    1.a.3) L-->Side of the square domain
    1.a.4) nparticle_per_cells_eachdir--> how many particles should be present in each cell in each direction.nparticle_per_cells_eachdir=1=> 1 particle in cell, nparticle_per_cells_eachdir=2=> 8 particles in cell
  1.b) run script and note down the number of particles generated (this would be printed on the output console)
  1.c) enter this number on the first line of mpm_particles.dat file generated in step 1.b
2. Open inputs_elasticdisk.in and fill in appropriate input data
3. Run exagoop executible with inputs_elasticdisk.in file
4. Plot of the energies with time can be obtained by running plot_energy.py script with input as the plot image filename.
5. Enjoy!
