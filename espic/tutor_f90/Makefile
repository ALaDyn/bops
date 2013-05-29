FFLAGS = -w -g #-free
DBFLAGS = -g

#  compiler

FC = gfortran

#  sources

SRCS =  espic.f  init.f  loadx.f  loadv.f density.f field.f \
	diagnostics.f push.f particle_bc.f plots.f histories.f

#  objects

OBJS = espic.o  init.o  loadx.o   loadv.o density.o field.o \
	diagnostics.o push.o particle_bc.o plots.o histories.o

.f.o:
	$(FC) -c $(FFLAGS) $<

#$(OBJS):
#	$(FC) $(FFLAGS) -c  $(SRCS)

master:   $(OBJS)
	@echo "\nLoading object files using flags:\t$(FFLAGS)\n"
	$(FC) $(FFLAGS) -o picf $(OBJS)


clean:
	rm -f picf *.o *~ *.dat

