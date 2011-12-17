This is the README file for the Particle-in-Cell code BOPS.

BOPS - Boosted Oblique-incidence Particle Simulation - 
is a one and three halves (1 spatial, 3 velocity coordinates: 1D3V)
particle-in-cell code originally created by myself together with
Tony Bell at in the Plasma Physics Group of Imperial College, London.

It employs a Lorentz transformation, or 'boost' along the target surface
to mimic the standard 2D, periodic-in-y geometry common to much of the 
early PIC work on resonance absorption in high-power laser-plasma interactions.
The technique was presented at the ECLIM conference in 1990, and later 
applied to absorption of femtosecond laser pulses on solid targets in 
PRL 68, 1535 (1992).  A longer description of the method including the 
transformation subtleties can be found in Phys. Plasmas 6, 947 (1999).
While restricting the simulations to a special class of problems - in which
the light and all its harmonics are reflected in the specular direction only -
the reduction from 2D->1D brings huge savings in computational effort and/or
increased spatial and temporal resolution.



Compiling BOPS
--------------

The directory structure resulting from unpacking the tar files should look
like this:

	src		...	containing the fortran77 source code
	doc		...	some documentation in html and ps
	run_scripts	...	sample scripts for running the code
	example_plots	...	sample output graphics
	

To compile, change directory to 'src' and select the Makefile corresponding 
closest to your machine/OS.  Eg. if you are running Linux, type

	cd src
	cp Makefile.linux Makefile
	
Then compile the code with:

	make bops
	

On machines other than Sun or Linux-PC, you will probably get complaints about the timing routine
'etime' in the file cputime.f.  If this happens, edit cputime.f and either 
replace etime with something which the compiler knows, or comment it out
altogether - this is not essential to run the code, but handy to know how
long it's going to run for.



Running BOPS
------------

Once compiled, go back to the base (or top) directory and edit one of the examples in
run_scripts (e.g. resabs).  Change the $BOPS variable to the directory where the bops.tar
file was unpacked (e.g. $HOME/bops) and the $RUN variable to where you want the data to be 
placed (e.g.: resabs1).  To run from the top directory, just type

	run_scripts/resabs
	
	
This will create a run directory 'resabs1' and start running the code.
All graphical output etc will be generated as a series of ASCII files
in the run directory.  Actual graphics are NOT supplied - this is up to 
the user!  Some sample plots roughly corresponding to the sample input 
files can be found in the examples directory.

A more detailed description of the input parameters can be found in

	doc/bops_man.html


------------------------------------------------------------------------
DISCLAIMER

This code is a scientific tool and should not be used for commercial purposes.
The author makes no guarantee of the correctness of results obtained with 
this code: it is up to the user to ensure that they make physical sense and
are not compromised by numerical instabilities etc.
This package may be freely distributed for academic purposes: the author 
requests, however, that due acknowledgement be given in any published work
which includes results obtained with the BOPS code.  (For example, a reference
to the Phys. Plasmas paper above).

------------------------------------------------------------------------
Last updated 4 September 2002  P. Gibbon
