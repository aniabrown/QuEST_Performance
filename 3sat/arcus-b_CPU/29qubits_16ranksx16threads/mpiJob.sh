#!/bin/bash

# set the number of nodes and processes per node. We are running one process on a single node
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1

#SBATCH --mem=50Gb
# uncomment if NUM_QUBITS - log2(NUM_NODES) > 30
####SBATCH --mem=100Gb

# set max wallclock time
#SBATCH --time=01:00:00

# set name of job
#SBATCH --job-name QuEST

# set queue
#SBATCH --partition=compute
#SBATCH --reservation=nqit

EXE=demo
export OMP_NUM_THREADS=16

module purge
module load mvapich2

. enable_arcus-b_mpi.sh

make clean
make

mpirun $MPI_HOSTS ./$EXE ../../3sat_input/q29/equ29qb.txt out.txt 1 0.78539816339 0
