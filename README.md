##redbKIT : a MATLAB(R) library for reduced-order modeling of parametrized PDEs

redbKIT is a [MATLAB](http://www.mathworks.com/products/matlab/) library (developed at [EPFL](https://www.epfl.ch/) - [Chair of Modeling and Scientific Computing](http://cmcs.epfl.ch/)) which implements some Reduced Basis (RB) methods for parametrized Partial Differential Equations (PDEs). In particular, it includes straightforward implementations of many of the algorithms presented in the book

>[**[QMN15] A. Quarteroni, A. Manzoni, F. Negri. Reduced Basis Methods for Partial Differential Equations**, Springer, 2015.](http://www.springer.com/us/book/9783319154305#aboutBook)

`redbKIT` consists of three main packages: [`RB_library`](https://github.com/redbKIT/redbKIT/tree/master/RB_library), [`FEM_library`](https://github.com/redbKIT/redbKIT/tree/master/FEM_library) and [`Problems`](https://github.com/redbKIT/redbKIT/tree/master/Problems), which are briefly described below.

#### RB_library

Implements in the MATLAB language many of the algorithms and methods presented in Chapters 3, 6, 7 and 10 of the book **[QMN15]**, such as
- proper orthogonal decomposition
- the greedy algorithm
- radial basis function interpolation of stability factors
- Galerkin and least-squares reduced basis methods
- the empirical interpolation method

The implementation is almost independent of the underlying high-fidelity approximation, provided that the high-fidelity model is described in a specified format. Here, we rely on the finite element method as high-fidelity approximation. See below.


#### FEM_library
Provides a flexible implementation of the finite element method for two-dimensional, stationary  advection-diffusion-reaction PDEs. The package can load linear triangular meshes either in the .msh format or in the [MATLAB Partial Differential Equation Toolbox(R)](http://www.mathworks.com/products/pde/index.html?s_tid=gn_loc_drop) format. Results can be either visualized in MATLAB or exported in the [VTK format](http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf) for post-processing with [Paraview](http://www.paraview.org/).

#### Problems
Contains a gallery of examples and applications. Many of them are described in Chapters 8, 9 and 10 of the book **[QMN15]**.


Download and Installation
-------

You can directly [**`download the current release`**](https://github.com/redbKIT/redbKIT/archive/master.zip) or clone the git repository

	$ git clone https://github.com/redbKIT/redbKIT.git /folder_path


`redbKIT` contains the following files and folders

[`FEM_library/`](https://github.com/redbKIT/redbKIT/tree/master/FEM_library)  [`LICENSE`](https://github.com/redbKIT/redbKIT/blob/master/LICENSE)  [`Problems/`](https://github.com/redbKIT/redbKIT/tree/master/Problems)  [`RB_library/`](https://github.com/redbKIT/redbKIT/tree/master/RB_library)  [`README.md`](https://github.com/redbKIT/redbKIT/blob/master/README.md)  [`setPath.m`](https://github.com/redbKIT/redbKIT/blob/master/setPath.m)

To check that your operating system and MATLAB installation are supported, and that the code is correctly working, please follow the following instructions. Start MATLAB and navigate to the `redbKIT` folder. Then, type in the MATLAB prompt

	>> setPath
	>> cd Problems
	>> test_all

The `test_all.m` functions launches a series of tests and generates a log file `test_log.txt`. If all the tests are marked as *passed*, you can start enjoying `redbKIT`. Otherwise, do not hesitate to contact us by email at <redbkit@gmail.com>.

Usage
-------

Start MATLAB(R) and navigate to the `redbKIT` folder. Then, type in the MATLAB prompt

	>> setPath

to add the [`FEM_library/`](https://github.com/redbKIT/redbKIT/tree/master/FEM_library) and [`RB_library/`](https://github.com/redbKIT/redbKIT/tree/master/RB_library) folders to the current path.
A welcome message should appear: you can now start using redbKIT!

A gallery of applications is provided in the folder [`Problems/`](https://github.com/redbKIT/redbKIT/tree/master/Problems).

If you want to speedup the finite element assembling operations, you can "mexify" the function
`ADR_mex_assembler.m` by running the following command in the MATLAB prompt ([`MATLAB Coder`](http://www.mathworks.com/products/matlab-coder/?refresh=true) needed):

	>> coder -build FEM_library/ADR_assembly.prj

Then replace `ADR_mex_assembler` with `ADR_mex_assembler_mex` at line 109
of `Assembler_2D.m`

An HTML documentation of redbKIT can be automatically generated using [`M2HTML`](http://www.artefact.tk/software/matlab/m2html/). To this end, you should download, extract and add to the path [`M2HTML`](http://www.artefact.tk/software/matlab/m2html/) using the following commands:

	>> url_M2HTML = 'http://www.artefact.tk/software/matlab/m2html/m2html.zip';
	>> unzip(url_M2HTML);
	>> addpath(genpath(strcat(pwd,'/m2html')));

Then generate the documentation

	>> m2html('mfiles',{'RB_library', 'RB_library/RBF_interpolation' 'FEM_library'}, ...
	   'htmldir','Documentation', 'recursive','off', 'global','on',...
	   'template','frame', 'index','menu', 'graph','on');


License
-------

**redbKIT is distributed under BSD 2-clause license**

Copyright (c) 2015, Ecole Polytechnique Fédérale de Lausanne (EPFL)
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
    Springer, 2015.

BibTex entry

    @book{QMN_RBspringer,
  	author  = {Quarteroni, A. and  Manzoni, A. and Negri, F.},
  	title   = {Reduced Basis Methods for Partial Differential Equations. An Introduction},
  	year    = {2015},
  	publisher = {Springer},
    address = {...}
	}

Contact
-------
Should you have any questions regarding `redbKIT`, do not hesitate to contact us by email at <redbkit@gmail.com>.
We also encourage contributions that can help to improve the code or the documentation, and to make `redbKIT` more useful.
