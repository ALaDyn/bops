#FFLAGS = -w -O2
FFLAGS = -w -g

od2gle: od2gle.o 
	gfortran $(OBJ) od2gle.o -o od2gle.exe

#	ln -sf `pwd`/od2gle ~/bin/od2gle

od2gle.o: od2gle.f lench.f
	gfortran $(FFLAGS) -c od2gle.f
clean:
	rm *.o od2gle.exe *% 

