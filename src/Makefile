f90=gfortran
# f90=/opt/mpich2/ch3_nemesis_tcp-pgi103/bin/mpif90
# f90=mpif90
objects=xsmcnp.o 
copt=-O3 -ffree-line-length-0
# copt=-fast
xsmcnp:    $(objects) 
	$(f90) $(copt) -o xsmcnp  $(objects)
xsmcnp.o:  xsmcnp.f90  
	$(f90) $(copt) -c  xsmcnp.f90
clean:
	rm -f *.o
	rm -f *.mod
	rm -f xsmcnp   
clear:
	rm -f *.o
	rm -f *.mod
