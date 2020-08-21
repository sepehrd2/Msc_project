# Msc_project

During my Master's project I was working on a C++ code, which I wrote from scratch, capable of simulating a quasi 2D model of biological cells containing rigid and flexible polymers.
Rigid plymers (or rods) were added to the code in order to model the rigidity of cytoskleton and actin filments. While, the flexible polymer was used to model the lipid bilayer. The code contains many functions and features capable of perfoming many tasks, including: 

1- non-bonded interactions (Lenard-Johnes and electostatics)

2- bond and angle forces

3- Area and lenght contraint forces to keep the area and volume of the cell fixed

4- Brownian dynamics and langevan dynamics simulaitons

5- Peridodic boundary conditoins (PBC)

6- Any additional biasing force 

The code was parallelized with openMP and ran on multiple CPU cores to get a better speed up.
