##redbKIT : a MATLAB(R) library for reduced-order modeling of parametrized PDEs

redbKIT is a [MATLAB](http://www.mathworks.com/products/matlab/) library (developed at [EPFL](https://www.epfl.ch/) - [Chair of Modeling and Scientific Computing](http://cmcs.epfl.ch/)) which implements some Reduced Basis (RB) methods for parametrized Partial Differential Equations (PDEs). In particular, it includes straightforward implementations of many of the algorithms presented in the book

>[**[QMN16] A. Quarteroni, A. Manzoni, F. Negri. Reduced Basis Methods for Partial Differential Equations. An Introduction**, Springer, 2016.](http://www.springer.com/us/book/9783319154305#aboutBook)

`redbKIT` consists of three main packages: [`RB_library`](https://github.com/redbKIT/redbKIT/tree/master/RB_library), [`FEM_library`](https://github.com/redbKIT/redbKIT/tree/master/FEM_library) and [`Problems`](https://github.com/redbKIT/redbKIT/tree/master/Problems), which are briefly described below.

#### RB_library

Implements in the MATLAB language many of the algorithms and methods presented in Chapters 3, 6, 7 and 10 of the book **[QMN16]**, such as
- proper orthogonal decomposition
- greedy algorithm
- radial basis function interpolation of stability factors
- Galerkin and least-squares reduced basis methods
- the empirical interpolation method and its discrete variant

The implementation is almost independent of the underlying high-fidelity approximation, provided that the high-fidelity model is described in a specified format. Here, we rely on the finite element method as high-fidelity approximation. See below.

As of Release 2.0, an implementation of the Matrix Discrete Empirical Interpolation Method (see [[NMA15]](http://www.sciencedirect.com/science/article/pii/S0021999115006543)) as well as an example of its application to the Helmholtz equation are provided.


#### FEM_library
Provides a flexible implementation for the following families of problems:
- 2D/3D steady and unsteady diffusion-transport-reaction equations, with P1/P2 finite elements.

- 2D/3D steady and unsteady  Navier-Stokes equations approximated by 
    - P2-P1 or P1Bubble-P1 finite elements for velocity and pressure spaces, respectively;
    - P1-P1 finite elements stabilized with the SUPG stabilization (implemented as in the framework of the Variational MultiScale Method).
For the steady case, Newton iterations are provided. For the unsteady case, time advancing is performed via BDF integrator, while the convection term can be treated either implicitly (with Newton subiterations) or with a semi-implicit scheme with extrapolation of the convective term.

- 2D/3D steady and unsteady structural problem. Both linear elasticity and nonlinear hyperelastic St. Venant Kirchhoff, nearly incompressible Neo-Hookean  and Raghavan-Vorp material models are implemented. For the steady case, a basic Newton algorithm with backtracking is provided. In the unsteady case, time advancing is performed via the generalized alpha-scheme suitably combined with Newton subiterations in the nonlinear case.


The package can load linear triangular meshes either in the .msh format or in the [MATLAB Partial Differential Equation Toolbox(R)](http://www.mathworks.com/products/pde/index.html?s_tid=gn_loc_drop) format. Results can be either visualized in MATLAB or exported in the (binary) [VTK format](http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf) for post-processing with [Paraview](http://www.paraview.org/).
                  

The assembly routines for the FE vectors and matrices are written in C code using suitable mex interfaces.  For this reason, before using the library you need to first compile the FE assemblers, as explained in the INSTALL file.  Loops over the elements are parallelized via OpenMP, while the global assembly of  sparse matrices from local contributes can be speeded-up by installing the [FAST package](http://user.it.uu.se/~stefane/freeware).
                  
The default linear solver is Matlab backslash (sparse direct solver). If available, also [MUMPS](http://mumps.enseeiht.fr/index.php?page=home) can be used in a straightforward fashion.  For moderately large size problems, a One-Level (geometric) Additive Schwarz preconditioner can be employed  in combination with a suitable iterative solver (usually gmres). In this case, mesh partitioning is done through the [Meshpart toolbox](http://www.cerfacs.fr/algor/Softs/MESHPART/) and [Metis library](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).

#### Problems
Contains a gallery of examples and applications. Many of them are described in Chapters 8, 9 and 10 of the book **[QMN16]**.


Download and Installation
-------

You can directly [**`download the current release`**](https://github.com/redbKIT/redbKIT/archive/master.zip) or clone the git repository

	$ git clone https://github.com/redbKIT/redbKIT.git /folder_path


`redbKIT` contains the following files and folders

[`FEM_library/`](https://github.com/redbKIT/redbKIT/tree/master/FEM_library)  [`LICENSE`](https://github.com/redbKIT/redbKIT/blob/master/LICENSE)  [`Problems/`](https://github.com/redbKIT/redbKIT/tree/master/Problems)  [`RB_library/`](https://github.com/redbKIT/redbKIT/tree/master/RB_library)  [`README.md`](https://github.com/redbKIT/redbKIT/blob/master/README.md)  [`setPath.m`](https://github.com/redbKIT/redbKIT/blob/master/setPath.m)  [`INSTALL`](https://github.com/redbKIT/redbKIT/blob/master/INSTALL)

To check that your operating system and MATLAB installation are supported, and that the code is correctly working, please follow these instructions. Start MATLAB and navigate to the `redbKIT` folder. Then, type in the MATLAB prompt

	>> setPath
	>> cd Problems
	>> test_all

The `test_all.m` functions launches a series of tests and generates a log file `test_log.txt`. If all the tests are marked as *passed*, you can start enjoying `redbKIT`. Otherwise, do not hesitate to contact us by email at <redbkit@gmail.com>.

`redbKIT` has been successfully tested on the following operating systems/MATLAB versions:

| MATLAB version / OS  |  Linux 64 bit  |  MAC 64 Bit | Windows 64 Bit | Windows 32 Bit |
|:--------------------:|:--------------:|:-----------:|:--------------:|:--------------:|
|           **R2015a** |                |      x      |                |                |
|           **R2014b** |        x       |      x      |       x        |                |
|           **R2014a** |        x       |             |                |                |
|           **R2013a** |        x       |             |                |                |
|           **R2012a** |        x       |             |                |                |
|           **R2011a** |        x       |             |                |                |
|           **R2010a** |        x       |             |                |                |
|           **R2009b** |                |             |                |        x       |

Usage
-------

Start MATLAB(R) and navigate to the `redbKIT` folder. Then, type in the MATLAB prompt

	>> setPath

to add the [`FEM_library/`](https://github.com/redbKIT/redbKIT/tree/master/FEM_library) and [`RB_library/`](https://github.com/redbKIT/redbKIT/tree/master/RB_library) folders to the current path.
A welcome message should appear: you can now start using redbKIT!

A gallery of applications is provided in the folder [`Problems/`](https://github.com/redbKIT/redbKIT/tree/master/Problems).

To check that the code is correctly working, please follow these instructions, please type in the MATLAB prompt

	>> setPath
	>> cd Problems
	>> test_all

The `test_all.m` functions launches a series of tests and generates a log file `test_log.txt`. If all the tests are marked as *passed*, you can start enjoying `redbKIT`. Otherwise, do not hesitate to contact us by email at <redbkit@gmail.com>.


License
-------

**redbKIT is distributed under BSD 2-clause license**

Copyright (c) 2015-2016, Ecole Polytechnique Fédérale de Lausanne (EPFL)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


How to cite
-------
Please use the following citation to reference `redbKIT`

	A. Quarteroni, A. Manzoni and F. Negri.
    Reduced Basis Methods for Partial Differential Equations. An Introduction.
    Unitext, vol. 92. Springer, 2016.

BibTex entry

    @book{QMN_RBspringer,
  	author  = {Quarteroni, A. and  Manzoni, A. and Negri, F.},
  	title   = {Reduced Basis Methods for Partial Differential Equations. An Introduction},
  	year    = {2016},
  	publisher = {Springer},
      series = {Unitext},
      volume = {92},
	}

Development
-------

redbKIT is developed and maintained by [`Federico Negri`](http://cmcs.epfl.ch/people/negri) (EPFL).

Paola Gervasio (Università degli Studi di Brescia ) is gratefully acknowledged for granting the use of parts of the finite element code MLife.


Contact
-------
Should you have any questions regarding `redbKIT`, do not hesitate to contact us by email at <redbkit@gmail.com>.
We also encourage contributions that can help to improve the code or the documentation, and to make `redbKIT` more useful.