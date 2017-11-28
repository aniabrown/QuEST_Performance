#!/bin/bash

# set the number of nodes and processes per node. We are running one process on a single node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# uncomment if NUM_QUBITS > 31
####SBATCH --mem=100Gb
# uncomment if NUM_QUBITS > 32
####SBATCH --mem=200Gb

# set max wallclock time
#SBATCH --time=04:00:00

# set name of job
#SBATCH --job-name QuEST

# set queue
#SBATCH --partition=mem6T
#SBATCH --exclusive

EXE=demo
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=128

module purge
module load intel-compilers

make clean
make

top -b -n 1 > top.txt

./$EXE ../../3sat_input/q29/equ29qb.txt out.txt 1 0.78539816339 0
