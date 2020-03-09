============================================================
 RTA-CUDA
 Version 1.0
 Sep 4 2018
 Author(s): Michael Strickland
============================================================

------------------------------------------------------------
 COMPILING
------------------------------------------------------------

To compile, simply type "make" from the command line.

This code requires a graphics card that is compatible with 
Nvidia's CUDA package and for the CUDA developer libraries
to be installed on your machine for compilation.

------------------------------------------------------------
 USAGE
------------------------------------------------------------

There is a file "input/params.txt" that contains all 
parameters that can be adjusted at runtime.  It includes 
comments describing the various options.  To run with the 
params.txt parameters simply type

./rta-cuda

If you would like to override some parameters in the 
params.txt file from the commandline the syntax is

./rta-cuda -<PARAMNAME1> <value1> ... -<PARAMNAMEn> <valuen>

------------------------------------------------------------
 OUTPUT
------------------------------------------------------------

All output is put in the "output" directory.  Right now the
only things that are dumped are "snapshots" of the a and b
functions at certain points in the iteration process.

------------------------------------------------------------
 CONTRIBUTORS
------------------------------------------------------------

Michael Strickland

------------------------------------------------------------
 LICENSE
------------------------------------------------------------

GNU General Public License (GPLv3)
See detailed text in license directory 

------------------------------------------------------------
 ATTRIBUTION
------------------------------------------------------------

If you use this code, you should cite (at minimum):

[1] W. Florkowski, R. Ryblewski, and M. Strickland, Testing 
viscous and anisotropic hydrodynamics in an exactly solvable 
case, Phys. Rev. C 88, 024903 (2013).

[2] W. Florkowski, R. Ryblewski, and M. Strickland, 
Anisotropic Hydrodynamics for Rapidly Expanding Systems, 
Nuclear Physics A 916, 249 (2013).

[3] M. Strickland, The non-equilibrium attractor for kinetic 
theory in relaxation time approximation, Forthcoming (2018).
