#
#    Makefile for BOPS3 and derivatives
#
#  Last update: October 2011
#
#------------
#  Compiler optimisation flags 
#  - tune these to your own machine!


# IBM Power 4
#------------
#FC = /usr/bin/mpxlf90 
#TUNE = -qarch=pwr4 -qtune=pwr4 -qhot -qstrict
#DEBUG= -qcheck -g -qfullpath
#TUNE = -g -qcheck
#MPITRACE=  -g -lmpihpm -lpmapi

#
#  Intel f90 compiler 
#--------------------
#FC = ifort
#TUNE = -O3 -FI -w90 -w95 
#TUNE = -FI -w90 -w95 -72
#TUNE = -O3 -traceback
#F90FLAGS = -w
#DEBUG = -g -check bounds -traceback 

# GNU fortran
#--------
#FC= ifort
#FC=gfortran
#TUNE = -O5  
DEBUG = -g  #-check bounds -traceback
#TUNE = -O3  

# cygwin or MINGW
#-------------------------------
FC = g95 
TUNE = -O3

# JUROPA@ FZJ
#--------------------------------
#INT_HOST = juropa
#FC = mpif90
#TUNE = -axSSE4.2 -O5
#DEBUG = -g -check bounds -traceback 
 


#LDFLAGS = -bnoquiet -lxlf90 
#LDFLAGS = $(MPITRACE)
#FFLAGS = $(DEBUG) $(F90FLAGS) 
FFLAGS = $(DEBUG) $(TUNE) 

#
MODULES.f90 = bopsvars.f90

SOURCES.f90 = \
addcurr.f90  emstore.f90  fourier.f90  init.f90     parload.f90  reinji.f90   trkini.f90\
addneut.f90  energy.f90   pbcs.f90     restrt.f90   udist.f90  ioutils.f90\
arinit.f90   escout.f90   gcon.f90     ixq.f90      phatrk.f90   rhoiav.f90   uescd.f90\
avlas.f90    esf1.f90     ghist.f90    ixq2.f90     bops.f90   savear.f90   uescdi.f90\
boost.f90    esfield.f90  grav.f90     ixq4.f90     plotrk.f90   scatrk.f90   uidist.f90\
bunif.f90    esft.f90     grps.f90     ixq8.f90     pout.f90     scramb.f90   urelmax.f90\
denprof.f90  f2v.f90      grxy.f90     jeden.f90    pouti.f90    sim2lab.f90  vdist.f90\
derdat.f90   field_ave.f90      gsnap.f90    jiden.f90    pshif.f90    splot.f90    wguide.f90\
diagnostics.f90   fftrc.f90    gtrak.f90    lab2sim.f90  pushc.f90    stamp.f90    xion.f90\
fill3v.f90   header.f90   laser.f90    pushp.f90    stotrk.f90   xq.f90\
eden.f90     filter1.f90  hmove.f90    load2u.f90   pushs.f90    surfie.f90  files.f90\
emfield.f90  filter2.f90  hop.f90      move.f90     rand.f90     teft.f90     gr3d.f90\
empert.f90   finabs.f90   iden.f90     movie.f90    reinj.f90    tmft.f90	currft.f90   rhoe_ave.f90\
cputime.f90  findxc.f90   empond.f90  push_pond.f90 fpond.f90 push_pondes.f90 \
em_helmholtz.f90 trisolve.f90 setup_arrays.f90 \
add_layer.f90  slab.f90 multislab.f90

# Derived parameters.
.SUFFIXES: .c .f90 

MODS = $(MODULES.f90:.f90=.o)
OBJS = $(SOURCES.f90:.f90=.o) $(MODS)

default: bops
all: mods bops



.f90.o:
	$(FC) -c $(FFLAGS) $<

#  explicit rules for f90 modules to avoid conflict with built-in ones
bopsvars.o:
	$(FC) -c $(FFLAGS) bopsvars.f90

mods: $(MODS)

bops:   $(OBJS) $(MODS)
	@echo "Loading object files using flags: $(FFLAGS)"
	$(FC) $(FFLAGS) -o bops $(OBJS)

debug:	$(OBJS)
	@echo "Loading object files\n\tusing flags\n$(FFLAGS)\n---------"
	@$(FC) $(FFLAGS) -o bops_debug $(OBJS)

ft:	ftana.o fftrc.o ioutils.o
	$(FC) $(FFLAGS) -o ftan ftana.o fftrc.o ioutils.o

clean:
	rm -f *.o bops *.mod

cleanall:   
	rm -f bops bops_debug ftan 

## to specify the dependencies, the corresponding object file is used 
## instead of the module file
## (can be automatically generated with "makedepf90 *.f90 -I../../modules")
addcurr.o: addcurr.f90 bopsvars.o
setup_arrays.o: setup_arrays.f90 bopsvars.o
init.o: names.h init.f90
