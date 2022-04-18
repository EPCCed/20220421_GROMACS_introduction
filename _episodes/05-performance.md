---
title: "Performance: getting the best bang for your (computational) buck"
teaching: 30
exercises: 30
questions:
- "How does our system perform as we scale up the number of cores on which we 
  run?"
- "How do we run hybrid MPI and OpenMP jobs on ARCHER2?"
- "Does adding OpenMP to MPI GROMACS affect performance?"
- "Does simultaneous multithreading (SMT) improve GROMACS performance on 
  ARCHER2?"
- "How does load balancing affect GROMACS performance?"
objectives:
- "Gain a basic understanding of key aspects of running molecular dynamics 
  simulations in parallel."
- "Run a hybrid MPI with OpenMP simuation on ARCHER2."
- "See how GROMACS performance is changed by including OpenMP."
- "Understand how to use simultaneous multithreading (SMT) on ARCHER2."
- "Learn how to disable GROMACS dynamic load-balancing and appreciate the 
  effect load balancing can have on performance."
keypoints:
- "Hybrid MPI with OpenMP does affect performance."
- "When running hybrid jobs, placement across NUMA regions is important."
---

## Aims

In this series of exercises, you will be running molecular dynamics 
simulations of the 5PEP protein. You can either use the GROMACS topology 
and portable binary files that you generated earlier in this session or use 
some pre-generated files (that should be equivalent to what you've already 
prepared. 

In the first of these exercises, you will be benchmarking how efficiently a 
small system will run using pure MPI as you increase the number of processors 
on which it runs. 

In the second exercise, you will explore how using hybrid MPI+OpenMPI 
methods can improve the runtime for this system.

The third exercise will have you studying the effects of multithreading.

Finally, in the 4th exercise, you will look at how dynamic load balancing can 
further reduce your simulation runtime.

## Exercise 1: MPI-only runs on ARCHER2

Before starting, you will need to get a copy of the exercises by running:

```bash
  svn checkout https://github.com/EPCCed/20220421_GROMACS_introduction/trunk/exercises
```

Once the file is downloaded, go into the `/exercises/performance` directory. 
This directory contains:

  - A Slurm submission file called `sub.slurm`
  - A GROMACS portable binary file called `npt.tpr` -- this file should be the 
    same as the one generated earlier in the course.

Below is a copy of the Slurm submission script that will run on a single node, 
using a single processor.

```
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
```
{: .language-bash}

Run this script to make sure that it works -- how quickly does the job
complete? You can see the walltime and performance of your code by running
`tail md.log`. The GROMACS `-nsteps 10000` flag should ensue that this 
simulation performs 10,000 steps (instead of the 500,000 steps defined in the 
`.mdp` file). Each step will simulate 1 ps of "real-life" time. The "Walltime" 
tells you how quickly the job ran on ARCHER2, and the "Performance" data will 
let you know how many nanoseconds you can run in a day, and how many hours it 
will take to run a nanosecond of simulation.

How do the "Walltime" and "Performance" data change as you increase the number
of cores being used? You can vary this by changing
`#SBATCH --tasks-per-node=1` to a higher number. Try filling out the table
 below:

 |Number of cores| Walltime | Performance (ns/day) | Performance (hours/ns) |
 |---------------|----------|----------------------|------------------------|
 |   1  | | | |
 |   2  | | | |
 |   4  | | | |
 |   8  | | | |
 |  16  | | | |
 |  32  | | | |
 |  64  | | | |
 | 128  | | | |
 | 256* | | | |
 | 512* | | | |

 ---
 **NOTE**

 Jobs run on more than one node will need to be run with constant
 `#SBATCH --tasks-per-node=128` but varying `#SBATCH --nodes=1`

 ---

## Exercise 2: Hybrid MPI and OpenMP jobs on ARCHER2

MPI and OpenMP are two different methods of programming parallel codes. MPI 
will allow you to parallelise your code across a distributed memory machine 
where each core is assumed to have its own memory that is not visible to any 
other core -- each core will have its own copy of every variable and passing 
information that should be shared from one core to another will require 
inter-core communication. OpenMP will allow you to parallelise your cores on 
shared-memory machines where each core has access to all of the data visible 
to any other core. Each has its advantages and disadvantages.

You can run GROMACS using a hybrid of MPI+OpenMP. Practically, this is done by 
splitting your job into a number of MPI ranks, each of which has a number of 
OpenMP threads assigned to it. If this is done correctly, you should get the 
benefits from having used both methods which will result in a noticeable 
speed-up.

When running hybrid MPI + OpenMP (with multiple threads) jobs you need to 
leave free cores between the parallel tasks launched using `srun` for the 
multiple OpenMP threads that will be associated with each MPI task.

You can use the options to `sbatch` to control how many parallel tasks are 
placed on each compute node and can use the `--cpus-per-task` option to set 
the stride between parallel tasks to the right value to accommodate the 
OpenMP threads. The value of `--cpus-per-task` should usually be the same as 
that for `OMP_NUM_THREADS`.

As an example, consider the job script below that runs across 2 nodes with
8 MPI tasks per node and 16 OpenMP threads per MPI task (so all 256 cores
are used).

```bash
#!/bin/bash

#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --reservation=shortqos
#SBATCH --account=ta059
#SBATCH --time=00:10:00

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16

#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

module load gromacs

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -nsteps 10000 -s npt.tpr
```

Each ARCHER2 compute node is made up of 8 NUMA (*Non Uniform Memory Access*) 
regions (4 per socket) with 16 cores in each region. Programs where the 
threads span multiple NUMA regions are likely to be *less* efficient so we 
recommend using thread counts that fit well into the ARCHER2 compute node 
layout. Effectively, this means one of the following options for nodes where 
all cores are used:

* 8 MPI tasks per node and 16 OpenMP threads per task: equivalent to 1 MPI 
  task per NUMA region
* 16 MPI tasks per node and 8 OpenMP threads per task: equivalent to 2 MPI 
  tasks per NUMA region
* 32 MPI tasks per node and 4 OpenMP threads per task: equivalent to 4 MPI 
  tasks per NUMA region
* 64 MPI tasks per node and 2 OpenMP threads per task: equivalent to 8 MPI 
tasks per NUMA region

### Instructions

For this exercise, you will start by comparing the performance of a simulation
that uses all of the cores on a node. Using the Slurm submission script from 
the first exercise as a template, try running simulations that use varying 
levels of MPI and OpenMP (making sure that the number of MPI ranks and OpenMP 
threads always multiply to 128 or less). How do the simulation times change as 
you increase the numbers change? How do these times change if you do not 
spread the threads over the NUMA regions as suggested?

You may find it helpful to fill out this table

| MPI Ranks | OpenMP threads | Walltime (s) | performance (ns/day) |
|-----------|----------------|--------------|----------------------|
|       128 |              1 | | |
|        64 |              2 | | |
|        42 |              3 | | |
|        32 |              4 | | |
|        25 |              5 | | |
|        16 |              8 | | |
|         8 |             16 | | |

## Exercise 3: Two hardware threads per core

The `--hint=nomultithread` asks SLURM to ignore the possibility of running
two threads per core. If we remove this option, this makes available 256
"cpus" per node (2 threads per core in hardware). To run 8 MPI tasks with
1 task per NUMA region running 32 OpenMP threads, the script would look like:

```
#!/usr/bin/env bash

#SBATCH --partition=standard
#SBATCH --time=00:10:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

#SBATCH --hint=multithread
#SBATCH --distribution=block:block

#SBATCH --cpus-per-task=32

module load gromacs

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

srun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -nsteps 10000 -s npt.tpr
```

> ## Multithreading and GROMACS?
>
> Staring with the MPI-only case first, how does enabling multithreading
> affect GROMACS performance?
>
> What about the performance of hybrid MPI+OpenMP jobs?
{: .challenge}

## Exercise 4: Load balancing

GROMACS performs dynamic load balancing when it deems necessary. Can you tell
from your md.log files so far whether it has been doing so, and what it
calculated the load imbalance was before deciding to do so?

To demonstrate the effect of the load imbalance counteracted by GROMACS’s
dynamic load balancing scheme, investigate what happens when this is turned
off by including the `-dlb no` option to `gmx_mpi mdrun`.

{% include links.md %}
