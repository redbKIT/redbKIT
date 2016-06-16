This folder contains an implementation of the acoustic horn problem presented in 

[NMA] F. Negri, A. Manzoni, D. Amsallem - Efficient model reduction of 
parametrized systems by matrix discrete empirical interpolation. Journal 
of Computational Physics 303, 2015.

Instructions:

0) Generate the meshes using gmsh and the .geo file given in the subfolder 
   "mesh" (the mesh used in [NMA] corresponds to the setting "fine2" in AcousticHorn.geo)

1) Launch 
   >> add_this_path

2) Navigate into the subfolder "FrequencyParam" which implements the consistency 
   test presented in Sect. 5.1 of [NMA]

   To build the FOM and the HROM using MDEIM+DEIM+POD(or greedy), launch
   >> main
   
   To assess the quality of the HROM vs the FOM, launch
   >> Online_Analysis
   which compute the output of interest over the parameter space.


3) Navigate into the subfolder "AllParam" which implements the 
   test presented in Sect. 5.3 of [NMA]
    
   Same as before.