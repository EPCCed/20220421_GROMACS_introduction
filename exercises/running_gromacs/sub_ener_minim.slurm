#!/bin/bash

#SBATCH --account=ta059
#SBATCH --partition=standard
#SBATCH --qos=short

#SBATCH --time=0:10:0
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

module load gromacs

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

gmx mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -s ener_minim.tpr
