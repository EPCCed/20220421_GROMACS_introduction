---
title: "ARCHER2 scheduler: Slurm"
teaching: 25
exercises: 20
questions:
- "How do I write job submission scripts?"
- "How do I control jobs?"
- "How do I find out what resources are available?"
objectives:
- "Understand the use of the basic Slurm commands."
- "Know what components make up and ARCHER2 scheduler."
- "Know where to look for further help on the scheduler."
keypoints:
- "ARCHER2 uses the Slurm scheduler."
- "`srun` is used to launch parallel executables in batch job submission scripts."
- "There are a number of different partitions (queues) available."
---

ARCHER2 uses the Slurm job submission system, or *scheduler*, to manage resources and how they are made
available to users. The main commands you will use with Slurm on ARCHER2 are:

* `sinfo`: Query the current state of nodes
* `sbatch`: Submit non-interactive (batch) jobs to the scheduler
* `squeue`: List jobs in the queue
* `scancel`: Cancel a job
* `salloc`: Submit interactive jobs to the scheduler
* `srun`: Used within a batch job script or interactive job session to start a parallel program

Full documentation on Slurm on ARCHER2 can be found in [the *Running Jobs on ARCHER2* section of the User
and Best Practice Guide](https://docs.archer2.ac.uk/user-guide/scheduler.html).

## Finding out what resources are available: `sinfo`

The `sinfo` command shows the current state of the compute nodes known to the scheduler:

```
auser@ln01:~> sinfo
```
{: .language-bash}
```
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
standard     up 1-00:00:00      2 drain$ nid[003548,005023]
standard     up 1-00:00:00      6  down$ nid[002546,003549,003583,004381,005123,005194]
standard     up 1-00:00:00      6  maint nid[001693,002462,003430,003835,004084,006112]
standard     up 1-00:00:00      6 drain* nid[001567,003550-003551,003995,004753,006080]
standard     up 1-00:00:00     13  down* nid[001200,001251,002326,002914,003185,003190,003395,003526,003598,003764,004984,005799,006210]
standard     up 1-00:00:00      1   comp nid002252
standard     up 1-00:00:00      6   drng nid[003992-003994,004752,004754-004755]
standard     up 1-00:00:00     12  drain nid[001117,001638,002085,002123,002315,003177,004778,005508-005511,006759]
standard     up 1-00:00:00    159   resv nid[001256-001383,001608-001637,001639]
standard     up 1-00:00:00   5634  alloc nid[001000-001116,001118-001199,001201-001213,001215-001250,001252-001255,001384-001566,001568-001607,001640-001692,001694-002057,002059-002084,002086-002122,002124-002193,002195-002251,002253,002262-002275,002277-002314,002316-002325,002327-002461,002463-002545,002547-002913,002915-003176,003178-003184,003186-003189,003191-003394,003396-003429,003431-003484,003486-003525,003527-003547,003552-003582,003584-003597,003599-003763,003765-003834,003836-003991,003996-004083,004085-004380,004382-004751,004756-004777,004779-004983,004985-005022,005024-005122,005124-005193,005195-005506,005512-005730,005732-005798,005800-006079,006081-006111,006113-006209,006211-006758,006760-006859]
standard     up 1-00:00:00      9   idle nid[002254-002261,002276]
standard     up 1-00:00:00      6   down nid[001214,002058,002194,003485,005507,005731]
highmem      up 1-00:00:00      1  down* nid002914
highmem      up 1-00:00:00    583  alloc nid[002756-002913,002915-003047,006376-006667]
serial       up 1-00:00:00      2   idle dvn[01-02]
```
{: .output}

There is a row for each node state and partition combination. The default output shows the following columns:

* `PARTITION` - The system partition
* `AVAIL` - The status of the partition - `up` in normal operation
* `TIMELIMIT` - Maximum runtime as `days-hours:minutes:seconds`: on ARCHER2, these are set using *QoS*
  (Quality of Service) rather than on partitions
* `NODES` - The number of nodes in the partition/state combination
* `STATE` - The state of the listed nodes (more information below)
* `NODELIST` - A list of the nodes in the partition/state combination

The nodes can be in many different states, the most common you will see are:

* `idle` - Nodes that are not currently allocated to jobs
* `alloc` - Nodes currently allocated to jobs
* `draining` - Nodes draining and will not run further jobs until released by the systems team
* `down` - Node unavailable
* `fail` - Node is in fail state and not available for jobs
* `reserved` - Node is in an advanced reservation and is not generally available
* `maint` - Node is in a maintenance reservation and is not generally available

> ## Easy viewing of output
> A lot of the Slurm commands can output very wide tables which can be hard to view in smaller terminals.
> You can try piping the output to `less -S` to get a scrollable view without line wrapping. Try it out by running
> for example `sinfo | less -S`. You can scroll with the arrow keys and exit the view by pressing `q`.
{: .callout}

If you prefer to see the state of individual nodes, you can use the `sinfo -N -l` command.

> ## Lots to look at!
> Warning! The `sinfo -N -l` command will produce a lot of output as there are nearly 6000 individual
> compute nodes on the full ARCHER2 system!
{: .callout}

```
auser@ln01:~> sinfo -N -l
```
{: .language-bash}
```
Mon Nov 22 15:54:46 2021
NODELIST   NODES PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT AVAIL_FE REASON
dvn01          1    serial        idle 256    2:64:2 515450        0      1 DVN,AMD_ none
dvn02          1    serial        idle 256    2:64:2 515450        0      1 DVN,AMD_ none
nid001000      1  standard   allocated 256    2:64:2 227328        0      1 COMPUTE, none
nid001001      1  standard   allocated 256    2:64:2 227328        0      1 COMPUTE, none
nid001002      1  standard   allocated 256    2:64:2 227328        0      1 COMPUTE, none
nid001003      1  standard   allocated 256    2:64:2 227328        0      1 COMPUTE, none
nid001004      1  standard   allocated 256    2:64:2 227328        0      1 COMPUTE, none
nid001005      1  standard   allocated 256    2:64:2 227328        0      1 COMPUTE, none
nid001006      1  standard   allocated 256    2:64:2 227328        0      1 COMPUTE, none
nid001007      1  standard   allocated 256    2:64:2 227328        0      1 COMPUTE, none

...lots of output trimmed...

```
{: .output}

> ## Explore a compute node
> Let’s look at the resources available on the compute nodes where your jobs will actually run. Try running this
> command to see the name, CPUs and memory available on the worker nodes (the instructors will give you the ID of
> the compute node to use):
> ```
> [auser@ln01:~> sinfo -n nid001005 -o "%n %c %m"
> ```
> {: .language-bash}
> This should display the resources available for a standard node. 
> 
> It is also possible to search nodes by state. Can you find all the free nodes in the system?
> > ## Solution
> > `sinfo` lets you specify the state of a node to search for, so to get all the free nodes in the system you can use:
> > ```
> > sinfo -N -l --state=idle
> > ```
> > More information on what `sinfo` can display can be found in the `sinfo` manual page, i.e. `man sinfo`
> {: .solution}

{: .challenge}

## Using batch job submission scripts

### Header section: `#SBATCH`

As for most other scheduler systems, job submission scripts in Slurm consist of a header section with the
shell specification and options to the submission command (`sbatch` in this case) followed by the body of
the script that actually runs the commands you want. In the header section, options to `sbatch` should 
be prepended with `#SBATCH`.

Here is a simple example script that runs the `xthi` program, which shows process and thread placement, across
two nodes.

```
#!/bin/bash
#SBATCH --job-name=my_mpi_job
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:0
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account={{site.gid}}
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

# Make the xthi software available
module load xthi

export OMP_NUM_THREADS=1

# srun to launch the executable
srun xthi
```
{: .language-bash}

The options shown here are:

* `--job-name=my_mpi_job` - Set the name for the job that will be displayed in Slurm output.
* `--nodes=2` - Select two nodes for this job.
* `--ntasks-per-node=128` - Set 128 parallel processes per node (usually corresponds to MPI ranks).
* `--cpus-per-task=1` - Number of cores to allocate per parallel process.
* `--time=0:10:0` - Set 10 minutes maximum walltime for this job.
* `--partition=standard` - Submit to the standard set of nodes.
* `--qos=standard` - Submit with the standard quality of service settings.
* `--account={{site.gid}}` - Charge the job to the `{{site.gid}}` budget
* `--hint=nomultithread` - Ensures that work is distributed amongst physical cores.
* `--distribution=block:block` - Defines how work is loaded onto the nodes and processors.

We will discuss the `srun` command further below.

> ## Using a reservation
> For this course we have a reservation in place, a set of nodes that have been set aside specifically
> for our use so we can avoid waiting in the queue. If you want to use it, you will need to use
> the login account you created under the {{site.gid}} project, and you should add the following line
> to your job scripts:
> ```
> #SBATCH --reservation=<reservation-name>
> ```
> {: .language-bash}
> The course instructor and helpers will be able to tell you the name of the reservation.
{: .callout}

### Submitting jobs using `sbatch`

You use the `sbatch` command to submit job submission scripts to the scheduler. For example, if the
above script was saved in a file called `test_job.slurm`, you would submit it with:

```
auser@ln01:~> sbatch test_job.slurm
```
{: .language-bash}
```
Submitted batch job 23996
```
{: .output}

Slurm reports back with the job ID for the job you have submitted

> ## What are the default for `sbatch` options?
> If you do not specify job options, what are the defaults for Slurm on ARCHER2? Submit jobs to find out
> what the defaults are for:
> 1. Partition and QoS?
> 2. Budget (or Account) the job is charged to?
> 3. Tasks per node?
> 4. Number of nodes?
> 5. Walltime? (This one is hard!)
> 
> > ## Solution
> > 
> > (1) Partition and QoS: None! We have to specify these or a job will be immediately rejected by `sbatch`.
> >
> > (2) Budget: This depends -- if you only have one budget associated with an account, or if you have set 
> > up a default budget, Slurm will use that as a default. Otherwise, Slurm will not let your job run.
> >
> > You can get the answers to 3. and 4. with a script like the following:
> > 
> > ```
> > #!/bin/bash
> > #SBATCH --job-name=my_mpi_job
> > #SBATCH --account={{site.gid}}
> > #SBATCH --partition=standard
> > #SBATCH --qos=standard
> >
> > echo "Nodes: $SLURM_JOB_NUM_NODES"
> > module load xthi
> > 
> > export OMP_NUM_THREADS=1
> > 
> > srun xthi
> > ```
> > {: .language-bash}
> >
> > The useful environment variable `$SLURM_NTASKS_PER_NODE` would normally tell us how many tasks we would
> > like to run on each node, but this is only set if the  `--ntasks-per-node` option has been given.
> > However, we can still see the default from xthi's output.
> > 
> > (3) Tasks per node: 1
> >
> > (4) Number of nodes: 1
> >
> > Getting the default time limit is more difficult - we need to use the `sacct` command to query the time limit set for
> > the job. For example, if the job ID was "12345", then we could query the time limit with:
> > 
> > ```
> > auser@ln01:~> sacct -o "TimeLimit" -j 12345
> > ```
> > {: .language-bash}
> > ```
> >  Timelimit 
> > ---------- 
> >   01:00:00
> > ```
> > {: .output}
> >
> > (5) Walltime: One hour.
> {: .solution}
{: .challenge}

### Checking progress of your job with `squeue`

You use the `squeue` command to show the current state of the queues on ARCHER2. Without any options, it
will show all jobs in the queue:

```
auser@ln01:~> squeue
```
{: .language-bash}
```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
...probably lots of jobs!...
```
{: .output}

You can use the `-u` option to look at the jobs being run by a particular user. So, to see the status of
your own jobs, if your user name is still `auser`, you would run:

```
auser@ln01:~> squeue -u auser
```
{: .language-bash}
```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
             12345 standard  my_mpi_j    auser PD       0:00      1 (Priority)
```
{: .output}

The `ST` column gives the state of each job. The three most common statuses you are likely to come across
are:

* `PD` - Pending. The job is still waiting to start.
* `R` - Running. The job is underway.
* `CG` - Completing. The job is in the process of finishing. Some processes may still be active and Slurm will be performing cleanup tasks.

Once the job has started, `NODELIST(REASON)` will give the list of nodes being used. Until then, it will
instead show the reason the job has not yet started. `Priority` means that higher priority jobs
are ahead of it in the queue and it will run once those have been scheduled. `Resources` means
that the resources requested are not yet available.

### Cancelling jobs with `scancel`

You can use the `scancel` command to cancel jobs that are queued or running. When used on running jobs
it stops them immediately. To cancel job `12345` you would run:

```
auser@ln01:~> scancel 12345
```
{: .language-bash}

<!-- Great content, not currently available
> ## Getting notified
> Slurm on ARCHER2 can also send e-mails to notify you when your job starts, ends, fails, etc. Can
> you find out how you would setup your job script to send you an e-mail when your job finishes and
> when it fails? Test your answer, does it work?
> > ## Solution
> > The option `--mail-type=END,FAIL` will send mailings to you when the job ends or fails. You can
> > also use the event `TIME_LIMIT` to notify you if a job reaches its walltime without finishing and
> > the events `TIME_LIMIT_50`, `TIME_LIMIT_80` and `TIME_LIMIT_90` to notify you when your job is
> > 50%, 80% and 90% of the way through the specified walltime.
> {: .solution}
{: .challenge}
-->
### Running parallel applications using `srun`

Once past the header section your script consists of standard shell commands required to run your
job. These can be simple or complex depending on how you run your jobs but even the simplest job
script usually contains commands to:

* Load the required software modules
* Set appropriate environment variables (you should always set `OMP_NUM_THREADS`, even if you are
  not using OpenMP when you should set this to `1`)

After this you will usually launch your parallel program using the `srun` command. At its simplest,
`srun` only needs 1 argument to specify the correct binding of processes to cores (it will use the
values supplied to `sbatch` to work out how many parallel processes to launch). In the example above,
our `srun` command simply looks like:

```
srun xthi
```
{: .language-bash}

> ## Underpopulation of nodes
> You may often want to *underpopulate* nodes on ARCHER2 to access more memory or more memory 
> bandwidth per task. Can you state the `sbatch` options you would use to run `xthi`:
> 
> 1. On 4 nodes with 64 tasks per node?
> 2. On 8 nodes with 2 tasks per node, 1 task per socket?
> 3. On 4 nodes with 32 tasks per node, ensuring an even distribution across the 8 NUMA regions
> on the node?
> 
> Once you have your answers run them in job scripts and check that the binding of tasks to 
> nodes and cores output by `xthi` is what you expect.
> 
> > ## Solution
> > 1. `--nodes=4 --ntasks-per-node=64`
> > 2. `--nodes=8 --ntasks-per-node=2 --ntasks-per-socket=1`
> > 3. `--nodes=4 --ntasks-per-node=32 --ntasks-per-socket=16 --cpus-per-task=4`
> {: .solution}
{: .challenge}

### Hybrid MPI and OpenMP jobs

When running hybrid MPI (with the individual tasks also known as ranks or processes) and OpenMP
(with multiple threads) jobs you need to leave free cores between the parallel tasks launched
using `srun` for the multiple OpenMP threads that will be associated with each MPI task.

As we saw above, you can use the options to `sbatch` to control how many parallel tasks are
placed on each compute node and can use the `--cpus-per-task` option to set the stride 
between parallel tasks to the right value to accommodate the OpenMP threads - the value
for `--cpus-per-task` should usually be the same as that for `OMP_NUM_THREADS`. Finally,
you need to specify `--threads-per-core=1` to ensure that the threads use physical
cores rather than SMT (hardware threading).

As an example, consider the job script below that runs across 2 nodes with 8 MPI tasks
per node and 16 OpenMP threads per MPI task (so all 256 cores across both nodes are used).

```
#!/bin/bash
#SBATCH --job-name=my_hybrid_job
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --threads-per-core=1
#SBATCH --time=0:10:0
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account={{site.gid}}
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

module load xthi

export OMP_NUM_THREADS=16

# Load modules, etc.
# srun to launch the executable

srun xthi
```
{: .language-bash}

Each ARCHER2 compute node is made up of 8 NUMA (*Non Uniform Memory Access*) regions (4 per socket) 
with 16 cores in each region. Programs where the threads of a task span multiple NUMA regions
are likely to be *less* efficient so we recommend using thread counts that fit well into the
ARCHER2 compute node layout. Effectively, this means one of the following options for nodes
where all cores are used:

* 8 MPI tasks per node and 16 OpenMP threads per task: equivalent to 1 MPI task per NUMA region
* 16 MPI tasks per node and 8 OpenMP threads per task: equivalent to 2 MPI tasks per NUMA region
* 32 MPI tasks per node and 4 OpenMP threads per task: equivalent to 4 MPI tasks per NUMA region
* 64 MPI tasks per node and 2 OpenMP threads per task: equivalent to 8 MPI tasks per NUMA region 

## STDOUT/STDERR from jobs

STDOUT and STDERR from jobs are, by default, written to a file called `slurm-<jobid>.out` in the
working directory for the job (unless the job script changes this, this will be the directory
where you submitted the job). So for a job with ID `12345` STDOUT and STDERR would be in
`slurm-12345.out`.

If you run into issues with your jobs, the Service Desk will often ask you to send your job
submission script and the contents of this file to help debug the issue.

If you need to change the location of STDOUT and STDERR you can use the `--output=<filename>`
and the `--error=<filename>` options to `sbatch` to split the streams and output to the named
locations.

## Other useful information

In this section we briefly introduce other scheduler topics that may be useful to users. We
provide links to more information on these areas for people who may want to explore these 
areas more. 

### Interactive jobs: `salloc` 

Similar to the batch jobs covered above, users can also run interactive jobs using the Slurm
command `salloc`. `salloc` takes the same arguments as `sbatch` but, obviously, these are
specified on the command line rather than in a job submission script.

When the job requested with `salloc` starts, the resources you requested are allocated for
your use. You will be returned to the command line 
and can now start parallel jobs on the compute nodes interactively with the `srun` command
in the same way as you would within a job submission script.

For example, to execute `xthi` across all cores on two nodes (1 MPI task per core and no
OpenMP threading) within an interactive job you would issue the following commands:

```
auser@ln01:~> salloc --partition=standard --qos=standard --nodes=2 --ntasks-per-node=128 --cpus-per-task=1 --distribution=block:block --hint=nomultithread --time=0:10:0 --account={{site.gid}}
salloc: Pending job allocation 12345
salloc: job 12345 queued and waiting for resources
salloc: job 12345 has been allocated resources
salloc: Granted job allocation 12345
salloc: Waiting for resource configuration
salloc: Nodes nid[004186-004187] are ready for job
auser@ln01:~> module load xthi
auser@ln01:~> srun xthi
```
{: .language-bash}
```
Node summary for    2 nodes:
Node    0, hostname nid004186, mpi 128, omp   1, executable xthi
Node    1, hostname nid004187, mpi 128, omp   1, executable xthi
MPI summary: 256 ranks
Node    0, rank    0, thread   0, (affinity =    0)
Node    0, rank    1, thread   0, (affinity =    1)
Node    0, rank    2, thread   0, (affinity =    2)
Node    0, rank    3, thread   0, (affinity =    3)
Node    0, rank    4, thread   0, (affinity =    4)
Node    0, rank    5, thread   0, (affinity =    5)
Node    0, rank    6, thread   0, (affinity =    6)
Node    0, rank    7, thread   0, (affinity =    7)
Node    0, rank    8, thread   0, (affinity =    8)
Node    0, rank    9, thread   0, (affinity =    9)
Node    0, rank   10, thread   0, (affinity =   10)
Node    0, rank   11, thread   0, (affinity =   11)
Node    0, rank   12, thread   0, (affinity =   12)
Node    0, rank   13, thread   0, (affinity =   13)
Node    0, rank   14, thread   0, (affinity =   14)
Node    0, rank   15, thread   0, (affinity =   15)
Node    0, rank   16, thread   0, (affinity =   16)
Node    0, rank   17, thread   0, (affinity =   17)
Node    0, rank   18, thread   0, (affinity =   18)

...long output trimmed...
```
{: .output}

Once you have finished your interactive commands, you exit the interactive job with `exit`:

```
auser@ln01:~> exit
exit
salloc: Relinquishing job allocation 12345
salloc: Job allocation 12345 has been revoked.
auser@ln01:~>
```
{: .language-bash}

### Interactive jobs: `srun`

An alternative way to run interactive jobs is with `srun` itself. This works a bit differently
from `salloc`. This method launches you into a new shell running on the lead node of the job.

```
auser@ln01:~> srun --partition=standard --qos=standard --nodes=2 --ntasks-per-node=128 --cpus-per-task=1 --time=0:10:0 --account={{site.gid}} --pty /bin/bash
srun: job 12345 queued and waiting for resources
srun: job 12345 has been allocated resources
auser@nid002158:~> module load xthi
auser@nid002158:~> srun --oversubscribe --ntasks=256 --distribution=block:block --hint=nomultithread xthi
```
{: .language-bash}

You'll notice that once the job has begun and you are in the shell running on the job's lead node,
the prompt will change to show the node's name. You'll also see in the example above that the `srun` used to launch `xthi`
now includes an extra an extra `--oversubscribe` option that is needed to allow Slurm to launch `xthi` on the
allocated resources. It also has to specify the total number of tasks to run and the `distribution` and `hint` options.

The shell will inherit the environment from the shell running on the login node, so you will also likely want to include
an `--export=none` option to the first `srun` to prevent this.

### The `short` QoS

Interactive jobs are very helpful ways to modify and run code on the fly while developing, but if you want to
submit small non-interactive test jobs, the `short` QoS is available. Compared to the `standard` QoS, jobs running
in the `short` QoS can user fewer nodes (32 vs 1024) and have a shorter max walltime (20 minutes vs 24 hours). It
can be easy to forget about these limits when using a different QoS, so if you find your job is immediately rejected,
you should try checking these first.

The `short` QoS is still accessed via the `standard` partition, but the QoS is set to `short`. For example, the
original xthi job from above could be run with:

```
#!/bin/bash
#SBATCH --job-name=my_short_mpi_job
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:0
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account={{site.gid}}
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

# Make the xthi software available
module load xthi

export OMP_NUM_THREADS=1

# srun to launch the executable
srun xthi
```
{: .language-bash}

Though the resource limits are lower than the other QoSes, using `short` can help you to run small test jobs
quickly.

### More information on running jobs on ARCHER2

There is a great deal of flexibility when running jobs on ARCHER2. For example, you may be interested in

* Using the high-memory nodes
* Running very long jobs
* Running job arrays
* Chaining jobs with dependencies
* Running subjobs using several nodes, or each with a fraction of one node

This information can be found in the documentation in the [ARCHER2 scheduler documentation](https://docs.archer2.ac.uk/user-guide/scheduler/).

<!-- Need to add information on the solid state storage and Slurm once it is in place

### Using the ARCHER2 solid state storage

-->


{% include links.md %}


