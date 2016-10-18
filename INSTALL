===========================================
redbKIT INSTALLATION
===========================================

Pre-requisites
--------------

For the basic installation, you need:
-> a reasonably recent version of Matlab
-> a supported C compiler, see e.g. http://it.mathworks.com/support/compilers/R2016a/
   (gcc on Unix and gcc/clang on Mac OS X work smoothly)

If you want to speed-up matrix assembly via multithreading, you also need:
-> a C compiler with OpenMP support, see e.g. http://it.mathworks.com/support/compilers/R2016a/
   (gcc is tested for both Unix and Mac OS X)

If you want to speed-up Matlab built-in "sparse" command by installing the FAST package,
you also need:
-> a supported C compiler for the serial version
-> a C compiler with OpenMP support for the multithreaded version

If you want to use parallel Additive Schwarz preconditioners, you also need
-> MATLAB Parallel Computing Toolbox
-> CMake


Basic Installation
------------------

Open Matlab and navigate to redbKIT root directory. Then type 

>> make

to compile the C-assembly routines and "mexify" some other files.


Installation with OpenMP enabled
--------------------------------

Open Matlab and navigate to redbKIT root directory. Then type 

>> make(1)

to compile the C-assembly routines with OpenMP enabled.


Install Metis for Additive Schwarz preconditioners
--------------------------------------------------

Open Matlab and navigate to the redbKIT directory FEM_library/Mesh
Then type 

>> install_metismex


Install FAST to speed-up Matlab "sparse" command
--------------------------------------------------

Open Matlab and navigate to the redbKIT directory FEM_library/Tools.

To install the serial version, type
>> install_FAST

To install the multithreaded version, type
>> install_FAST(1)

In both cases you will be required to accept FAST's license.
See http://user.it.uu.se/~stefane/freeware.


Install MUMPS (parallel sparse linear solver)
---------------------------------------------

The LinearSolver class provides the possibility to use MUMPS (a MUltifrontal Massively 
Parallel sparse direct Solver) as sparse linear solver. 
To this end, please follow these steps (only for Linux machines):

1 - Open a terminal and navigate to the folder redbKIT/FEM_library/LinearSolver/Mumps

2 - Launch the install script

    $ sh INSTALL.sh

This script automatically download and install into a local folder Metis, Scotch, OpenBlas 
and finally MUMPS libraries. To check that the MUMPS Matlab interface was correcly compiled, 
open Matlab, navigate to the folder 
   redbKIT/FEM_library/LinearSolver/Mumps/Libraries/MUMPS/MUMPS_5.0.1/MATLAB 
and run the simple_example.
