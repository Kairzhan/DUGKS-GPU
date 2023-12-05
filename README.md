# DUGKS-CUDA
Multi-GPU implementation of the DUGKS method[^1][^2].

[^1]: Guo Z. et al, Discrete unified gas kinetic scheme for all Knudsen number flows: Low-speed isothermal case, Phys. Rev. E 88, 033305, 2013, https://doi.org/10.1103/PhysRevE.88.033305
[^2]: Bo Y. et al, DUGKS simulations of three-dimensional Taylor–Green vortex flow and turbulent channel flow, Computers & Fluids, Volume 155, 2017, https://doi.org/10.1016/j.compfluid.2017.03.007

## Installation

To compile the code installation of the [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) is required. The code was developed on a workstation with NVIDIA HPC SDK version 22.3. For visualization purposes one may use [Paraview](https://www.paraview.org/) or [Gnuplot](http://www.gnuplot.info/). Some postprocessing scripts require [Python3](https://python.org).

Repository contains two versions of the code: single-GPU CUDA-Fortran version and Multi-GPU CUDA-C version. The CUDA-C version also can be used in single GPU run.

To compile the codes tweak appropriate Makefiles and run commands `make clean` followed by `make`. By default executable is named `main`. 

## Usage

After compiling the codes use the following command to run simulation using CUDA-C version:

```
mpirun -np 2 ./main
```
where it is assumed that computer has two GPUs and code will automatically use both GPUs through MPI parallelization.

To run the code on GPU clusters please refer to local GPU job submissino policy. Example jobfile used on [Qiming](https://newshub.sustech.edu.cn/en/html/202007/26934.html) supercomputer looks like following:

```
#!/bin/bash
#BSUB -q hgx
#BSUB -J channel
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -gpu "num=4/host"
#BSUB -W 48:00
#BSUB -o stdout_%J.out
#BSUB -e stderr_%J.err

module load nvhpc/22.11

/share/nvhpc/Linux_x86_64/22.11/comm_libs/mpi/bin/mpirun -np 4 ./main
```
where a node with 4 GPUs was requested.

## Examples

## Code details

### main.cu
