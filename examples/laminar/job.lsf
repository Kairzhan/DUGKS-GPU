#!/bin/bash
#BSUB -q hgx
#BSUB -J laminar
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -gpu "num=1/host"
#BSUB -W 12:00
#BSUB -o stdout_%J.out
#BSUB -e stderr_%J.err

if [ -d old ]; then
  mv std* old/
fi

module load nvhpc/22.11
nvidia-smi -L

/share/nvhpc/Linux_x86_64/22.11/comm_libs/mpi/bin/mpirun -np 1 ./main
