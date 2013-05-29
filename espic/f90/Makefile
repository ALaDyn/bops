FFLAGS = -w -g 
DBFLAGS = -g

#  compiler

FC = gfortran
#FC = f90

#  sources

SRCS =  espic.f90  init.f90  loadx.f90  loadv.f90 density.f90 field.f90 \
	diagnostics.f90 push.f90 boundaries.f90 plots.f90 histories.f90

#  objects

#OBJS = espic.o  init.o  loadx.o   loadv.o density.o field.o \
#	diagnostics.o push.o boundaries.o plots.o histories.o

OBJS = $(SRCS:%.f90=%.o)

#  rules for compiling objects


%.o:  %.f90
#	@ echo "Compiling $<" $(VERBOSE)
	$(FC) -c $(FFLAGS) -o $@  $*.f90


espic:   $(OBJS)
	$(FC) $(FFLAGS) -o espic $(OBJS)


clean:
	rm -f pices *.o *~ *.data

