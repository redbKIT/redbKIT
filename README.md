##redbKIT : a MATLAB(R) library for reduced-order modeling of parametrized PDEs

redbKIT is a MATLAB library (developed at [EPFL](https://www.epfl.ch/) -- [Chair of Modeling and Scientific Computing](http://cmcs.epfl.ch/)) which implements some Reduced Basis (RB) methods for parametrized Partial Differential Equations (PDEs). In particular, it includes straightforward implementations of many of the algorithms presented in the book [**Reduced Basis Methods for Partial Differential Equations (Springer, 2015)** by **A. Quarteroni, A. Manzoni and F. Negri.**](http://www.springer.com/us/book/9783319154305#aboutBook), such as
- the proper orthogonal decomposition
- the greedy algorithm
- radial basis function interpolation of stability factors
- Galerkin and least-squares reduced basis methods
- the empirical interpolation method

In all cases, the construction of the RB approximation is based

Download
-------

You can directly download the latest release at

	https://github.com/redbKIT/redbKIT/archive/master.zip

or clone the git repository

	git clone https://github.com/redbKIT/redbKIT.git /folder_path


`redbKIT` contains the following files and folders

[`FEM_library/`](https://github.com/redbKIT/redbKIT/tree/master/FEM_library)  [`LICENSE`](https://github.com/redbKIT/redbKIT/blob/master/LICENSE)  [`Problems/`](https://github.com/redbKIT/redbKIT/tree/master/Problems)  [`RB_library/`](https://github.com/redbKIT/redbKIT/tree/master/RB_library)  [`README.md`](https://github.com/redbKIT/redbKIT/blob/master/README.md)  [`setPath.m`](https://github.com/redbKIT/redbKIT/blob/master/setPath.m)

Usage
-------

Start MATLAB (R) and navigate to the `redbKIT` folder. Then, type in the MATLAB prompt

	setPath

to add the [`FEM_library/`](https://github.com/redbKIT/redbKIT/tree/master/FEM_library) and [`RB_library/`](https://github.com/redbKIT/redbKIT/tree/master/RB_library) to the current path.
A welcome message should appear: you can now start using redbKIT!

A gallery of applications is provided in the folder [`Problems/`](https://github.com/redbKIT/redbKIT/tree/master/Problems).

If you want to speedup the finite element assembling operations, you can "mexify" the function
`ADR_mex_assembler.m` by running the following command in the MATLAB prompt ([`MATLAB Coder`](http://it.mathworks.com/products/matlab-coder/?refresh=true) needed):

	coder -build FEM_library/ADR_assembly.prj

Then replace `ADR_mex_assembler` with `ADR_mex_assembler_mex` at line 109
of `Assembler_2D.m`

An HTML documentation of redbKIT can be automatically generated using [`M2HTML`](http://www.artefact.tk/software/matlab/m2html/). To this end, you should download, extract and add to the path [`M2HTML`](http://www.artefact.tk/software/matlab/m2html/) using the following commands:

	url_M2HTML = 'http://www.artefact.tk/software/matlab/m2html/m2html.zip';
	unzip(url_M2HTML);
	addpath(genpath(strcat(pwd,sslash,'m2html')));

Then generate the documentation

	m2html('mfiles',{'RB_library', 'RB_library/RBF_interpolation' 'FEM_library'}, ...
	   'htmldir','Documentation', 'recursive','off', 'global','on',...
	   'template','frame', 'index','menu', 'graph','on');


License
-------
redbKIT is distributed under BSD 2-clause license by Ecole Polytechnique Fédérale de Lausanne (EPFL)

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


Contact
-------
Should you have any questions regarding `redbKIT`, do not hestitate to contact
us by email at <redbkit@gmail.com>
