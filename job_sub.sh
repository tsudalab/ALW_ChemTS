#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe impi 2049 #Number of cores. minimum 2
#$ -jc pcc-skl.168h #pcc-normal.72h # pcc-skl.168h
## Load the Intel MPI environment variables
. /fefs/opt/x86_64/Gaussian/envset.sh

ulimit -s unlimited

#g16 H2O
#. /fefs/opt/x86_64/intel/parallel_studio_xe_2017/impi/2017.2.174/bin64/mpivars.sh

#Please set your environment
source /home/terayama/.bashrc
source activate py2
export KERAS_BACKEND=tensorflow

#Please set the above number of cores to the -n option.
mpiexec -bootstrap sge -n 2049 python mpi_thread_chemts_tree_vl.py
source deactivate
