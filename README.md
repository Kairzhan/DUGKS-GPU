# DUGKS-GPU
Multi-GPU implementation of the DUGKS method[^1][^2].

[^1]: Guo Z. et al, Discrete unified gas kinetic scheme for all Knudsen number flows: Low-speed isothermal case, Phys. Rev. E 88, 033305, 2013, https://doi.org/10.1103/PhysRevE.88.033305
[^2]: Bo Y. et al, DUGKS simulations of three-dimensional Taylor–Green vortex flow and turbulent channel flow, Computers & Fluids, Volume 155, 2017, https://doi.org/10.1016/j.compfluid.2017.03.007

## Installation

To compile the code installation of the [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) is required. The code was developed on a workstation with NVIDIA HPC SDK version 22.3. For visualization purposes one may use [Paraview](https://www.paraview.org/) or [Gnuplot](http://www.gnuplot.info/). Some postprocessing scripts require [Python3](https://python.org).

Repository contains two versions of the code: single-GPU CUDA-Fortran version and Multi-GPU CUDA-C version. The CUDA-C version also can be used in a single GPU run.

To compile the codes tweak appropriate Makefiles and run commands `make clean` followed by `make`. By default executable is named `main`. 

## Usage

After compiling the codes use the following command to run simulation using CUDA-C version:

```
mpirun -np 2 ./main
```
where it is assumed that computer has two GPUs and code will automatically use both GPUs through MPI parallelization.

To run the code on GPU clusters please refer to the local GPU job submission policy. Shown below is an example jobfile used on the [Qiming](https://newshub.sustech.edu.cn/en/html/202007/26934.html) supercomputer:

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
This jobfile requests a node with 4 GPUs.

## Examples

### Laminar channel flow

Laminar channel flow is considered in this example. Channel center max. velocity is set to 0.1. To run this case use following commands:

```bash
# Cd into example folder
cd examples/laminar

# check run configuration
# vim config.txt

# Run simuation on 1 GPU ...
mpirun -np 1 ./main

# edit the file params.py, which contains output frequency,
# this script creates the file VELMEAN.dat
# which contains various profiles
python3 ProcessProfiles.py # Averages profiles

# then create plots
gnuplot velprof.plt
```
Output is saved in file U.pdf.

### Channel flow
In this example turbulent flow in the channel is considered. Shear velocity based Reynolds number is **Re<sub>τ</sub>=180**. For this case, perturbed initial values for distribution function are provided and can be downloaded from the [Zenodo dataset](https://zenodo.org/doi/10.5281/zenodo.10377131). After downloading the file [ddf00000000.dat.gz](https://zenodo.org/records/10377132/files/ddf00000000.dat.gz?download=1) place it to the directory from which simulation will be started later. Note that this directory should also contain executable **main** and **config.txt** files. Main parameters which control this simulation are (see config.txt):

```
restart=1
nrestart=0
nsteps=6000000
ndump=6000000
```

Where *restart=1* indicates that initial condition will be loaded from dump file, and the dump file contains data for timestep 0 (*nrestart=0*). Total of 6 000 000 steps will be performed, and every 6 000 000 steps code will create dumpfiles for restarting purposes.

After compiling the code, place the executable main into the directory examples/channel, and run the simulation using mpirun or using batch system of your computer cluster. Output will be saved into the following files:

```
diag.dat
profiles-2.dat
s2X0000001.tec
s2Y0000001.tec
s2Z0000001.tec
s2X0005000.tec
s2Y0005000.tec
s2Z0005000.tec
...
```

In a summary, perform the following steps to run channel flow case:

```bash
# Cd into example folder
cd examples/channel

# Download input dataset
curl -o ddf00000000.dat.gz https://zenodo.org/records/10377132/files/ddf00000000.dat.gz?download=1
gunzip ddf00000000.dat.gz

# Run simuation on 4 GPUs ...
mpirun -np 4 ./main

# edit the file params.py, which contains output frequency,
python3 ProcessProfiles.py # Averages profiles

# Download benchmark data & split it into separate files
curl -o ch180.dat https://jaxa-dns-database.jaxa.jp/channelflow/ch180.dat
cat ch180.dat | awk 'NR >  80 && NR <145' > JAXA1.dat 
cat ch180.dat | awk 'NR > 147 && NR <212' > JAXA2.dat 

# then create plots
gnuplot contour.plt 
gnuplot log_plot.plt
gnuplot stress_plot.plt
gnuplot rms_plot.plt
```
<!--
# ... then copy gnuplot and postprocessing files
cp ../../plots/*.plt .
cp ../../plots/*.py .
-->

Following files will then contain plots:

```
Contour.pdf
Ulog.pdf
Ustress.pdf
Urms.pdf
```

<!-- ## Code details

### main.cu
-->
