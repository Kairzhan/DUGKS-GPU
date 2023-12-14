FC=nvfortran
#FFLAGS=-O0 -g -C -v
#FFLAGS=-O0 -g -C
FFLAGS=-O3 -r8 #-Mcuda=ptxinfo -Mcuda=keepptx -Minfo=accel #-r8
INCLUDE=
LDFLAGS=-L /opt/nvidia/hpc_sdk/Linux_x86_64/22.3/math_libs/lib64/ 

all: clean main

main: main.cuf dugks.cuf
	$(FC) $(FFLAGS) $(LDFLAGS) dugks.cuf main.cuf

clean:
	rm -f *.o *.mod main *.tec *.dat
