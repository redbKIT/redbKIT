This folder contains an implementation of the acoustic horn problem presented in 

[NMA] F. Negri, A. Manzoni, D. Amsallem - Efficient model reduction of 
parametrized systems by matrix discrete empirical interpolation. Journal 
of Computational Physics 303, 2015.

In particular, we consider here an affine representation of the test case
of Sect. 5.1 of [NMA] involving only one parameter: the frequency.

Instructions:

0) Generate the meshes using gmsh and the .geo file given in the subfolder 
   "mesh"

1) To build the FOM and the ROM using POD (by default) or Greedy, launch
   >> Reduced_Horn 

2) To assess the quality of the ROM vs the FOM, launch
   >> Online_Analysis
   which compute the output of interest over the parameter space.
