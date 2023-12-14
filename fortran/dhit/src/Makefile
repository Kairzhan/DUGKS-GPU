FC=nvfortran
#FFLAGS=-O0 -g -C -v
#FFLAGS=-O0 -g -C
FFLAGS=-O3 -lm -lfftw3  #-Mcuda=ptxinfo -Mcuda=keepptx -Minfo=accel #-r8
INCLUDE=
LDFLAGS=-L /opt/nvidia/hpc_sdk/Linux_x86_64/22.3/math_libs/lib64/ 

all: clean main

main: main.cuf dugks.cuf stat.f90 dft.f90
	$(FC) $(FFLAGS) $(LDFLAGS) dft.f90 stat.f90 dugks.cuf main.cuf

clean:
	rm -f *.o *.mod main *.tec *.dat
