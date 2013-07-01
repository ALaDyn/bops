

This is the README file for the Particle-in-Cell code BOPS.
----------------------------------------------------------

BOPS - Boosted Oblique-incidence Particle Simulation - 
is a one and three halves (1 spatial, 3 velocity coordinates: 1D3V)
particle-in-cell code originally created by Paul Gibbon together with
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
increased spatial and temporal resolution.  This type of code has
become a 'workhorse' for high-intensity laser-matter interaction
studies, giving relatively easy access to some extremely nonlinear,
kinetic plasma phenomena, such as hot electron generation, ion acceleration and
high-harmonic generation from solid surfaces.


1. Prerequisites
----------------
Versions 2.X and below are written in Fortran 77 and can be compiled with the GNU g77 or Intel compilers.
The run scripts are designed for Unix systems and will work on a Linux PC. 

Windows users:

The code can be run 'under' Windows using an appropriate Fortran environment, but I have not 
attempted this myself yet.  (Instructions for successful builds v. welcome).

Here is a tried-and-tested method which may require a bit of space disc space.

i) Install  CYGWIN (www.cygwin.com).  This creates a fully-fledge Unix environment under Windows.  
In addition to the default tools/packages you will also need:
	- Gnu compilers gcc, g77 etc (devel package)
	- make (devel)
	- vi, emacs (editors - optional, but handy!)

ii) Download and unpack the bopsXX.tar.gz file with an archiving tool (eg PowerArchiver: www.powerarchiver.com).
You can put this anywhere, but a convenient location is your 'home' directory under CYGWIN, e.g.:

	C:\cygwin\home\foo\bops


iii) Open a Cygwin terminal/shell and 'cd' to the bops directory.



2. Compiling BOPS
-----------------

The directory structure resulting from unpacking the tar files should look
like this:

	src/		...	containing the fortran77 source code
	doc/		...	some documentation in html and ps
	run_scripts/	...	sample scripts for running the code
	tools/          ...     postprocessing tools
	example_plots/	...	sample output graphics
	

To compile, change directory to 'src' and select the Makefile corresponding 
closest to your machine/OS (in src/makes/).  Linux is selected by default.

	cd src
	(cp makes/Makefile.machine Makefile)
	
Then compile the code with:

	make bops
	

On machines other than Sun or Linux-PC, you will probably get complaints about the timing routine
'etime' in the file cputime.f.  If this happens, edit cputime.f and either 
replace etime with something which the compiler knows, or comment it out
altogether - this is not essential to run the code, but handy to know how
long it's going to run for.



3. Running BOPS
---------------

Once compiled, go back to the base (or top) directory and edit one of the examples in
run_scripts (e.g. resabs).  Change the $BOPS variable to the directory where the bops.tar
file was unpacked (e.g. $HOME/bops) and the $RUN variable to where you want the data to be 
placed (e.g.: resabs1).  To run from the base directory, just type

	run_scripts/resabs
	
	
This will create a new run directory 'resabs1' and start executing the code.
All graphical output etc will be generated as a series of ASCII files
in the run directory.  Actual graphics are NOT supplied at present,
though there is a postprocessor in the tools/ directory (od2gle; copy the script odpp
to the run directory) which will create GLE-readable output.
You can get this program (Graphics Layout Engine, Windows and Linux) from:

http://glx.sourceforge.net

Otherwise, this is up to the user - gnuplot or xmgrace will usually 
suffice to get started.  
Some sample plots roughly corresponding to the sample input 
files can be found in the examples directory.

To do a series or parameter study, you might prefer to modify the
script to sit inside a 'project' directory and create subdirectories for 
each run.

A more detailed description of the input parameters can be found in

	doc/bops_man.pdf


4. Example scripts
------------------

resabs		Long scale-length, classical resonance absorption demo
snells_law	Refractive index transition (underdense plasma)
gb_prl92	Vacuum heating demo: steep density gradient, fixed ions
foil		Thin foil
foil+ramp	Foil + exponential leading ramp
hhg		High-harmonic generation from plasma surface
	


------------------------------------------------------------------------
DISCLAIMER

This code is a scientific tool and should not be used for commercial
or military purposes. The author makes no guarantee of the correctness
of results obtained with this code: it is up to the user to ensure
that they make physical sense and are not compromised by numerical
instabilities etc.

This package may be freely distributed for academic purposes: the author 
requests, however, that due acknowledgement be given in any published work
which includes results obtained with the BOPS code.  (For example, a reference
to the Phys. Plasmas paper above).

------------------------------------------------------------------------
Last updated Feb 2007  P. Gibbon











