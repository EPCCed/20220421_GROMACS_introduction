#!/bin/bash

#SBATCH --job-name=GMX_test
#SBATCH --account=ta059
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --time=0:10:0

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --distribution=block:block
#SBATCH --hint=nomultithread

module load gromacs

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun -n 1 \
  gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -nsteps 10000 -s npt.tpr
