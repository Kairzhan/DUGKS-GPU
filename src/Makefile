# by default use MPI compiler from NVIDIA HPC 
CC=mpic++

# host cpu architecture
#CPU_ARCH=zen3
CPU_ARCH=native

# uncomment PREC4 for single precision
CFLAGS=-O3 --c++20 -tp=${CPU_ARCH} -gpu=ptxinfo -g -Wall -Wno-unused-parameter -Minfo #-DPREC4=1

LFLAGS=-lm 

main: main.cu parallel.o cuda.o lbm_d3q19.o numerics.o checks.o dugks.o output.o config.o
	$(CC) $(CFLAGS) $(LFLAGS) -o main main.cu parallel.o cuda.o numerics.o checks.o dugks.o lbm_d3q19.o output.o

cuda.o: cuda.cu cuda.h
	$(CC) $(CFLAGS) -c -o cuda.o cuda.cu

numerics.o: numerics.cu numerics.h lbm_d3q19.o
	$(CC) $(CFLAGS) -c -o numerics.o numerics.cu lbm_d3q19.o

parallel.o: parallel.cu parallel.h
	$(CC) $(CFLAGS) -c -o parallel.o parallel.cu

constants.o: constants.cu precision.h
	$(CC) $(CFLAGS) -c -o constants.o constants.cu

checks.o: checks.cu
	$(CC) $(CFLAGS) -c -o checks.o checks.cu

lbm_d3q19.o: lbm_d3q19.cu
	$(CC) $(CFLAGS) -c -o lbm_d3q19.o lbm_d3q19.cu

dugks.o: dugks.cu lbm_d3q19.o 
	$(CC) $(CFLAGS) -c -o dugks.o dugks.cu lbm_d3q19.o

output.o: output.cu
	$(CC) $(CFLAGS) -c -o output.o output.cu

config.o: config.cu config.h
	$(CC) $(CFLAGS) -c -o config.o config.cu

clean:
	rm *.o main *~
