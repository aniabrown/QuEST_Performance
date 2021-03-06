#!/bin/bash

# set the number of nodes and processes per node. We are running one process on a single node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2

#SBATCH --mem=50Gb
# uncomment if NUM_QUBITS - log2(NUM_NODES) > 30
####SBATCH --mem=100Gb

# set max wallclock time
#SBATCH --time=00:20:00

# set name of job
#SBATCH --job-name QuEST

# set queue
#SBATCH --partition=compute

NUM_QUBITS=30
EXE=demo
export KMP_AFFINITY=verbose,compact
export OMP_NUM_THREADS=8

module purge
module load mvapich2

. enable_arcus-b_mpi.sh

make clean
make

mpirun $MPI_HOSTS -bind-to socket -map-by socket ./$EXE $NUM_QUBITS
