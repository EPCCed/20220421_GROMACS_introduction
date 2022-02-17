---
title: "Debugging and profiling on ARCHER2"
teaching: 30
exercises: 30
start: True
questions:
- "What debugging tools are available on ARCHER2 and how can I access them?"
- "What profiling tools are available on ARCHER2 and how can I access them?"
- "Where can I find more documentation on and get help with these tools?"
objectives:
- "Know what tools are available to help you debug and profile parallel applications on ARCHER2."
- "Know where to get further help."
keypoints:
- "The main debugging tool on ARCHER2 is *gdb4hpc*"
- "The main profiling tool on ARCHER2 is *CrayPat*"
---

ARCHER2 has a range of debugging and profiling software available. In this section we provide a brief
overview of the tools available, their applicability and links to more information. A detailed tutorial
on the use of these tools is beyond the scope of this course but if you are interested in this, then
you may be interested in the following courses offered by the ARCHER2 service:

* Performance Analysis Workshop

The [Cray Performance Measurement and Analysis Tools User Guide](https://pubs.cray.com/bundle/Cray_Performance_Measurement_and_Analysis_Tools_User_Guide_644_S-2376/page/About_the_Cray_Performance_Measurement_and_Analysis_Tools_User_Guide.html)
and the ARCHER2 [debugging](https://docs.archer2.ac.uk/user-guide/debug/) and
[profiling](https://docs.archer2.ac.uk/user-guide/profile/) documentation will also be useful.

## Debugging tools overview

The following debugging tools are available on ARCHER2:

* **gdb4hpc** is a command-line tool working similarly to [gdb](https://www.gnu.org/software/gdb/)
  that allows users to debug parallel programs. It can launch parallel programs or attach to ones
  already running and allows the user to step through the execution to identify the causes of any
  unexpected behaviour. Available via ``module load gdb4hpc``.
* **valgrind4hpc** is a parallel memory debugging tool that aids in detection of memory leaks and
  errors in parallel applications. It aggregates like errors across processes and threads to simply
  debugging of parallel appliciations. Available via ``module load valgrind4hpc``.
* **STAT** generate merged stack traces for parallel applications. Also has visualisation tools.
  Available via ``module load cray-stat``.
* **ATP** scalable core file and backtrace analysis when parallel programs crash. Available via
  ``module load atp``.
* **CCDB** Cray Comparative Debugger. Compare two versions of code side-by-side to analyse differences.
  Available via ``module load cray-ccdb``.

## Using gdb4hpc to debug an application

For this exercise, we'll be debugging a short program using gdb4hpc. To start, we'll grab a copy of a buggy code from ARCHER2:

```bash
wget {{site.url}}{{site.baseurl}}/files/gdb4hpc_exercise.c
```

You can look at the code if you want -- you might even be able to debug it by inspection (but that defeats the purpose of this exercise). When you're ready, compile the code using the C compiler wrappers and the debugging flag `-g`:

```bash
 cc -g gdb4hpc_exercise.c -o gdb_exercise
 ```

You can choose a different name for your executable, but I'll be using `gdb_exercise` through this exercise for consistency -- if you use a different name, make the appropriate change wherever you see `gdb_exercise`.

We'll be using ``gdb4hpc`` to go through this program and see where errors might arise.

Setup your environment, load and launch ``gdb4hpc``:

```bash
 module load gdb4hpc
 gdb4hpc
```

You will get some information about this version of the program and, eventually, you will get a command prompt:

```
gdb4hpc 4.12 - Cray Line Mode Parallel Debugger
With Cray Comparative Debugging Technology.
Copyright 2007-2021 Hewlett Packard Enterprise Development LP.
Copyright 1996-2016 University of Queensland. All Rights Reserved.

Type "help" for a list of commands.
Type "help <cmd>" for detailed help about a command.
dbg all>
```

We will use ``launch`` to start an application within gdb4hpc. For now, we want to run our simulation on a single process, so we will type:

```bash
 dbg all> launch --launcher-args="--account={{site.gid}} --partition=standard --qos=short --time=0:10:0 --tasks-per-node=1 --cpus-per-task=1 --exclusive --export=ALL" $my_prog{1} ./gdb_exercise
```

This will launch an ``srun`` job for ``gdb_exercise`` on one of the compute nodes. The name ``my_prog`` will be used by ``gdb4hpc`` as a reference to this particular run of the program -- you will not be able to launch another program using this name, but you can use any name you want instead. Once the run is started, you can reference it by prepending it with a ``$`` sign, so ``$my_prog`` in this case. The number in the curly brackets ``{1}`` indicates the number of processes this job will be using (it's  1 here). You could use a larger number if you wanted. If you call for more processes than available on a single compute node, ``gdb4hpc`` will launch the program on an appropriate number of nodes. Note though that the more cores you ask for, the slower ``gdb4hpc`` will be to launch the tasks once the job has begun. We use ``--launcher-args`` to pass all the ``SBATCH`` options we would normally provide in a job script through to the job launcher.

Once the program is launched, gdb4hpc will load up the program and begin to run it. You will get output to screen something that looks like:

```
Starting application, please wait...
Creating MRNet communication network...
Waiting for debug servers to attach to MRNet communications network...
Timeout in 400 seconds. Please wait for the attach to complete.
Number of dbgsrvs connected: [1];  Timeout Counter: [0]
Finalizing setup...
Launch complete.
my_prog{0}: Initial breakpoint, main at /PATH/TO/gdb4hpc_exercise.c:9
```

The line number at which the initial breakpoint is made (in the above example,
line 9) corresponds to the first line within the `main` function.

Once the code is loaded, you can use various commands to move through your code. The following lists and describes some of the most useful ones:

* ``help`` -- Lists all gdb4hpc commands. You can run ``help COMMAND_NAME`` to learn more about a specific command (*e.g.* ``help launch`` will tell you about the launch command
* ``list`` -- Will show the current line of code and the 9 lines following. Repeated use of ``list`` will move you down the code in ten-line chunks.
* ``next`` -- Will jump to the next step in the program for each process and output which line of code each process is on. It will not enter subroutines. Note that there is no reverse-step in gdb4hpc.
* ``step`` -- Like ``next``, but this will step into subroutines.
* ``up`` -- Go up one level in the program (*e.g.* from a subroutine back to main).
* ``print var`` -- Prints the value of variable ``var`` at this point in the code.
* ``watch var`` -- Like print, but will print whenever a variable changes value.
* ``backtrace`` -- Prints the stack trace for each process.
* ``quit`` -- Exits gdb4hpc.

For now, we will look at `list`, `next`, `print`, and `watch`. Running:

```bash
 dbg all> list
```

should output the first 10 lines of `main`:

```
 my_prog{0}: 9
 my_prog{0}: 10	  // Initiallise MPI environment
 my_prog{0}: 11	  MPI_Init(NULL,NULL);
 my_prog{0}: 12
 my_prog{0}: 13	  // Get processor rank
 my_prog{0}: 14	  int rank;
 my_prog{0}: 15	  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 my_prog{0}: 16
 my_prog{0}: 17	  int count = rank + 1;
 my_prog{0}: 18
```

Repeating `list` will bring show you the next 10 lines, *etc.*.

At the moment, we are at the start of the program. By running `next`, we will move to the next executable part of the program@

```bash
 dbg all> next
```
```
 my_prog{0}: main at /PATH/TO/gdb4hpc_exercise.c:13
```

Running `list` again will output the ten lines from 11-20. We can jump forward multiple lines by running `next N` -- by replacing *N* with a number, we will jump down *N* executable lines within our code. The `next` command will not allow us to move from one subroutine or function to another.

We can see on line 15 that there is a variable `count` about to be set. If we type:

```bash
 dbg all> print count
```

The current value of variable `count` is printed to screen. If we progress the code past line 15 and print this variable value again, it has changed to 1. If we wanted, we could have used the `watch` command to get a notification whenever the value of the variable changes.

> ## Exercise
> What happens if you keep using `next` and `list`?
> > ## Solution
> > The program will move from executable line to executable line until it reaches line 18, at which point the program is exited due to an MPI error
> {: .solution}
{: .challenge}

Let's now ``quit`` gdb4hpc, start it again, and try launching across multiple processes:

```bash
 dbg all> launch --launcher-args="--account={{site.gid}} --partition=standard --qos=short --time=0:10:0 --tasks-per-node=2 --cpus-per-task=1 --exclusive --export=ALL" $my_prog{2} ./gdb_exercise
```

> ## Exercise
> The code seems to be trying to send the variable `count` from one process to another. Follow `count` (using `watch`) and see how it changes throughout the code. What happens?
> > ## Solution
> > Eventually, both processes will hang: process 0 hangs at an `MPI_Barrier` on line 19 and is stuck waiting for process 1 to reach its barrier. Process 1 is stuck at an `MPI_Recv` on line 21. Further investigation shows that it is waiting for an `MPI_Send` that does not exist -- the source is process 1 (which has not sent anything) and the tag is `1` (there is no MPI_Send with this tag).
> {: .solution}
{: .challenge}

Let's `quit` our program, fix that bug, and go back into `gdb4hpc`. Again, we'll launch our program on 2 processes, and again, we'll watch the variable `count`. This time, both processes are able to get the same value for the variable `count`. There is one final part of the code to look at -- process 0 will try to get the sum of all even numbers between 0 and 20. However, the program begins to hang when process 0 reached line 28 (process 1 also hangs, but it's already at the `MPI_Finalize` part of the routine so we don't need to worry about it). Once we reach this hang, we can't easily keep going. Let's stop this debugging session and restart it by using `release`:

```bash
 dbg all> release $my_prog
 dbg all> launch --launcher-args="--account={{site.gid}} --partition=standard --qos=short --tasks-per-node=2 --cpus-per-task=1 --exclusive --export=ALL" $my_new_prog{2} ./gdb_exercise
```

This time, instead of `next`, we will use `step` -- this does the same as `next` with the added feature that we can go into functions and subroutines where applicable. As the new bug appears to come from the `sum_even` function, let's see where exactly the program hangs.

> ## Exercise
> Having `step`ed into the `sum_even` function, can you find where the code hangs?
> > ## Solution
> > The `i++` should be brought outside of the `if` part of the `while` loop. Changing this will make the code work fully.
> {: .solution}
{: .challenge}

There are other ways to run gdb4hpc. For example. you can attach it to a job which is already running if you would like to examine
its execution. You can also start an interactive session with `salloc` as we saw in the previous session and then `launch` from
within gdb4hpc directly to the job without having to specify the `--launcher-args` options. This may be more efficient if you need
to repeatedly restart the executable, as you won't have to wait for nodes to be allocated every time.

## Profiling tools overview

Profiling on ARCHER2 is provided through the Cray Performance Measurement and Analysis Tools (CrayPat). These have
a number of different components:

* **CrayPat** the full-featured program analysis tool set. CrayPat in turn consists of the following major components.
    * pat_build, the utility used to instrument programs
    * the CrayPat run time environment, which collects the specified performance data during program execution
    * pat_report, the first-level data analysis tool, used to produce text reports or export data for more sophisticated analysis
* **CrayPat-lite** a simplified and easy-to-use version of CrayPat that provides basic performance analysis information automatically, with a minimum of user interaction.
* **Reveal** the next-generation integrated performance analysis and code optimization tool, which enables the user to correlate performance data captured during program execution directly to the original source, and identify opportunities for further optimization.
* **Cray PAPI** components, which are support packages for those who want to access performance counters
* **Cray Apprentice2** the second-level data analysis tool, used to visualize, manipulate, explore, and compare sets of program performance data in a GUI environment.


See the [Cray Performance Measurement and Analysis Tools User Guide](https://pubs.cray.com/bundle/Cray_Performance_Measurement_and_Analysis_Tools_User_Guide_644_S-2376/page/About_the_Cray_Performance_Measurement_and_Analysis_Tools_User_Guide.html).

## Using CrayPat Lite to profile an application

Let's grab and unpack a toy code for training purposes. To do this, we'll use `wget`.

```
auser@ln01:~>  wget {{site.url}}{{site.baseurl}}/files/nbody-par.tar.gz
```
{: .language-bash}

To extract the files from a `.tar.gz` file, we run the command `tar -xvf filename.tar.gz`:
```
auser@ln01:~> tar -xzf nbody-par.tar.gz
```
{: .bash}

Load CrayPat-lite module (`perftools-lite`)
```
auser@ln01:~> module load perftools-lite
```
{: .bash}

and move into the new ``nbody-par`` directory and compile the application normally
```
auser@ln01:~> cd nbody-par
auser@ln01:~/nbody-par> make
```
{: .bash}
```
cc -DMPI -c main.c -o main.o
cc -DMPI -c utils.c -o utils.o
cc -DMPI -c serial.c -o serial.o
cc -DMPI -c parallel.c -o parallel.o
cc -DMPI main.o utils.o serial.o parallel.o -o nbody-parallel.exe
INFO: creating the PerfTools-instrumented executable 'nbody-parallel.exe' (lite-samples) ...OK
```
{: .output}

As the output of the compilation says, the executable `nbody-parallel.exe` binary has been instrumented by CrayPat-lite. A batch script called ``run.slurm`` is in the directory. Change the account to ``{{site.gid}}`` so it can be run. You can also use the course node reservation (check with the instructor or helpers) or the ``short`` QoS.

Once our job has finished, we can get the performance data summarized at the end of the job STDOUT.
```
auser@ln01:~/nbody-par> less slurm-out.txt
```
 {: .language-bash}
```
CrayPat/X:  Version 21.02.0 Revision ee5549f05  01/13/21 04:13:58
.
N Bodies =                     10240
Timestep dt =                  2.000e-01
Number of Timesteps =          10
Number of MPI Ranks =          64
BEGINNING N-BODY SIMULATION
SIMULATION COMPLETE
Runtime [s]:              3.295e-01
Runtime per Timestep [s]: 3.295e-02
interactions:             10
Interactions per sec:     3.182e+09

#################################################################
#                                                               #
#            CrayPat-lite Performance Statistics                #
#                                                               #
#################################################################

CrayPat/X:  Version 21.02.0 Revision ee5549f05  01/13/21 04:13:58
Experiment:                  lite  lite-samples
Number of PEs (MPI ranks):     64
Numbers of PEs per Node:       64
Numbers of Threads per PE:      1
Number of Cores per Socket:    64
Execution start time:  Tue Nov 23 15:18:45 2021
System name and speed:  nid005209  2.250 GHz (nominal)
AMD   Rome                 CPU  Family: 23  Model: 49  Stepping:  0
Core Performance Boost:  All 64 PEs have CPB capability

Avg Process Time:       0.47 secs
High Memory:         4,326.8 MiBytes     67.6 MiBytes per PE
I/O Write Rate:   701.173609 MiBytes/sec


...lots of output trimmed...

```
{: .output}

Instructions are given at the end of the output on how to see more detailed results and how to
visualise them graphically with the Cray Apprentice2 tool.

## Using CrayPat to profile an application

We are now going to use the full CrayPat tools. To do so, we first need to load the required modules

```
auser@ln01:~/nbody-par>  module unload perftools-lite
auser@ln01:~/nbody-par>  module load perftools
```
{: .language-bash}

After loading the modules, we need to recompile the application
```
auser@ln01:~/nbody-par> make clean; make
```
{: .bash}

Once the application has been built, we need to instrument the binary. We do this with `pat_build`
```
auser@ln01:~/nbody-par> pat_build nbody-parallel.exe
```
{: .bash}

and a new binary called `nbody-parallel.exe+pat` it will be generated. In fact, this is the binary that we need
to run in order to obtain the performance data. Replace the ``srun`` line in the Slurm script with the following
```
srun --cpu-bind=cores ./nbody-parallel.exe+pat -n 10240 -i 10 -t 1
```
 {: .language-bash}

After the job has finished, files are stored in an experiment data directory with the following format: `exe+pat+PID-node[s|t]` where:

* `exe`: The name of the instrumented executable
* `PID`: The process ID assigned to the instrumented executable at runtime
* `node`: The physical node ID upon which the rank zero process was executed
* `[s|t]`: The type of experiment performed, either `s` for sampling or `t` for tracing

for example, in our case a new directory called `nbody-parallel.exe+pat+192028-1341s` could be created. It is now time to obtain the performance report. We do this with the `pat_report` command and the new created directory
```
auser@ln01:~/nbody-par> pat_report nbody-parallel.exe+pat+192028-1341s
```
{: .bash}

This will command generate a full performance report and can generate a large amount of data, so you may wish to capture the data in an output file, either using a shell redirect like `>`,  or we could choose to see only some reports. If we want to see only a profile report by function we can do

```
auser@ln01:~/nbody-par> pat_report -v -O samp_profile nbody-parallel.exe+pat+192028-1341s
```
{: .bash}

```
...run details...

Table 1:  Profile by Function

  Samp% | Samp | Imb. |  Imb. | Group
        |      | Samp | Samp% |  Function
        |      |      |       |   PE=HIDE

 100.0% | 33.4 |   -- |    -- | Total
|-------------------------------------------------------
|  86.4% | 28.8 |  1.2 |  3.9% | USER
||------------------------------------------------------
||  86.4% | 28.8 |  1.2 |  3.9% | compute_forces_multi_set
||======================================================
|  12.0% |  4.0 |   -- |    -- | MPI
||------------------------------------------------------
||   4.3% |  1.4 |  1.6 | 52.9% | MPI_File_write_all
||   2.9% |  1.0 |  0.0 |  1.6% | MPI_Recv
||   2.9% |  1.0 |  0.0 |  3.2% | MPI_File_open
||   1.8% |  0.6 |  2.4 | 81.0% | MPI_Sendrecv
||======================================================
|   1.6% |  0.5 |  1.5 | 74.6% | MATH
||------------------------------------------------------
||   1.6% |  0.5 |  1.5 | 74.6% | sqrt
|=======================================================

...more run details...
```
{: .output}

The table above shows the results from sampling the application. Program functions are separated out into different types, `USER` functions are those defined by the application, `MPI` functions contains the time spent in MPI library functions, `ETC` functions are generally library or miscellaneous functions included. `ETC` functions can include a variety of external functions, from mathematical functions called in by the library to system calls.

The raw number of samples for each code section is shown in the second column and the number as an absolute percentage of the total samples in the first. The third column is a measure of the imbalance between individual processors being sampled in this routine and is calculated as the difference between the average number of samples over all processors and the maximum samples an individual processor gave in this routine.

Another useful table can be obtained profiling by Group, Function, and Line
```
auser@ln01:~/nbody-par> pat_report -v -O samp_profile+src nbody-parallel.exe+pat+192028-1341s
```
{: .bash}

```
...run details...

Table 1:  Profile by Group, Function, and Line

  Samp% | Samp | Imb. |  Imb. | Group
        |      | Samp | Samp% |  Function
        |      |      |       |   Source
        |      |      |       |    Line
        |      |      |       |     PE=HIDE

 100.0% | 33.4 |   -- |    -- | Total
|--------------------------------------------------------------
|  86.4% | 28.8 |   -- |    -- | USER
||-------------------------------------------------------------
||  86.4% | 28.8 |   -- |    -- | compute_forces_multi_set
3|        |      |      |       |  work/ta043/nbody-par/parallel.c
||||-----------------------------------------------------------
4|||   3.6% |  1.2 |  2.8 | 70.6% | line.138
4|||   2.1% |  0.7 |  2.3 | 77.8% | line.140
4|||   1.6% |  0.5 |  2.5 | 83.1% | line.141
4|||   1.1% |  0.4 |  1.6 | 83.3% | line.142
4|||   2.6% |  0.9 |  2.1 | 72.5% | line.144
4|||   3.5% |  1.2 |  2.8 | 72.2% | line.145
4|||   1.5% |  0.5 |  1.5 | 77.0% | line.147
4|||  29.7% |  9.9 |  6.1 | 38.7% | line.148
4|||  18.2% |  6.1 |  4.9 | 45.3% | line.149
4|||  11.8% |  3.9 |  5.1 | 57.1% | line.150
4|||  10.5% |  3.5 |  4.5 | 57.1% | line.151
||||===========================================================
||=============================================================
|  12.0% |  4.0 |   -- |    -- | MPI
||-------------------------------------------------------------
||   4.3% |  1.4 |  1.6 | 52.9% | MPI_File_write_all
||   2.9% |  1.0 |  0.0 |  1.6% | MPI_Recv
||   2.9% |  1.0 |  0.0 |  3.2% | MPI_File_open
||   1.8% |  0.6 |  2.4 | 81.0% | MPI_Sendrecv
||=============================================================
|   1.6% |  0.5 |  1.5 | 74.6% | MATH
||-------------------------------------------------------------
||   1.6% |  0.5 |  1.5 | 74.6% | sqrt
|==============================================================

...more run details...
```
{: .output}

If we want to profile by Function and Callers, with Line Numbers then
```
auser@ln01:~/nbody-par> pat_report -O ca+src nbody-parallel.exe+pat+192028-1341s
```
{: .bash}

```
...run details...

Table 1:  Profile by Function and Callers, with Line Numbers

  Samp% | Samp | Group
        |      |  Function
        |      |   Caller
        |      |    PE=HIDE

 100.0% | 33.4 | Total
|--------------------------------------------------------------
|  86.4% | 28.8 | USER
||-------------------------------------------------------------
||  86.4% | 28.8 | compute_forces_multi_set
3|        |      |  run_parallel_problem:parallel.c:line.81
4|        |      |   main:main.c:line.81
||=============================================================
|  12.0% |  4.0 | MPI
||-------------------------------------------------------------
||   4.3% |  1.4 | MPI_File_write_all
3|        |      |  distributed_write_timestep:parallel.c:line.218
4|        |      |   run_parallel_problem:parallel.c:line.73
5|        |      |    main:main.c:line.81
||   2.9% |  1.0 | MPI_Recv
3|        |      |  run_parallel_problem:parallel.c:line.59
4|        |      |   main:main.c:line.81
||   2.9% |  1.0 | MPI_File_open
3|        |      |  run_parallel_problem:parallel.c:line.59
4|        |      |   main:main.c:line.81
||   1.8% |  0.6 | MPI_Sendrecv
3|        |      |  run_parallel_problem:parallel.c:line.75
4|        |      |   main:main.c:line.81
||=============================================================
|   1.6% |  0.5 | MATH
||-------------------------------------------------------------
||   1.6% |  0.5 | sqrt
3|        |      |  run_parallel_problem:parallel.c:line.81
4|        |      |   main:main.c:line.81
|==============================================================

...more run details...
```
{: .output}

The run will generate two more files in the output directory, one with the extension `.ap2` which holds the same data as the report data (`.xf`) but in the post processed form. The other file is called `build-options.apa` and is a text file with a configuration for generating a traced (as opposed to sampled) experiment. The APA (Automatic Program Analysis) configuration is a targeted trace, based on the results from our previous sampled experiment. You are welcome and encouraged to review this file and modify its contents in subsequent iterations, however in this first case we will continue with the defaults.

This `build-options.apa` file acts as the input to the `pat_build` command and is supplied as the argument to the `-O` flag.

```
pat_build -O nbody-parallel.exe+pat+192028-1341s/build-options.apa
```
{: .bash}

This will produce a third binary with extension `+apa`. This binary should once again be run on the back end, so the submission script should be modified and the name of the executable changed to `nbody-parallel.exe+apa`.
Similarly to the sampling process, a new directory called `exe+apa+PID-node[s|t]` will be generated by the application, which should be processed by the `pat_report` tool. The output format of this new directory is similar to the one obtained with sampling, but now this includes `apa` and `t` to indicate that this is a tracing experiment. For instance, with our code we could get a new directory called `nbody-parallel.exe+apa+114297-5209t`.
```
auser@ln01:~/nbody-par> pat_report nbody-parallel.exe+apa+114297-5209t
```
{: .bash}

```
...run details...

Table 1:  Profile by Function Group and Function

  Time% |     Time |     Imb. |  Imb. | Calls | Group
        |          |     Time | Time% |       |  Function
        |          |          |       |       |   PE=HIDE

 100.0% | 0.402275 |       -- |    -- | 690.0 | Total
|-----------------------------------------------------------------
|  71.3% | 0.286987 | 0.003933 |  1.4% |   1.0 | USER
||----------------------------------------------------------------
||  71.3% | 0.286987 | 0.003933 |  1.4% |   1.0 | main
||================================================================
|  27.6% | 0.111100 |       -- |    -- | 686.0 | MPI
||----------------------------------------------------------------
||  14.0% | 0.056473 | 0.003279 |  5.6% |   1.0 | MPI_Recv
||   9.2% | 0.037024 | 0.000008 |  0.0% |  21.0 | MPI_File_write_all
||   2.3% | 0.009178 | 0.000001 |  0.0% |   1.0 | MPI_File_open
||   1.6% | 0.006504 | 0.000890 | 12.2% | 640.0 | MPI_Sendrecv
||================================================================
|   1.0% | 0.004188 |       -- |    -- |   3.0 | MPI_SYNC
|=================================================================

...more tables and details...
```
{: .output}

The new table above is the version generated from tracing data instead of the previous sampling data table. This version makes available true timing information (average per processor) and the number of times each function is called.

We can also get very important performance data from the hardware (HW) counters
```
auser@ln01:~/nbody-par> pat_report -v -O profile+hwpc nbody-parallel.exe+apa+114297-5209t
```
{: .bash}
```
...run details...

Table 1:  Profile by Function Group and Function

Group / Function / PE=HIDE


==============================================================================
  Total
------------------------------------------------------------------------------
  Time%                                        100.0%
  Time                                       0.402275 secs
  Imb. Time                                        -- secs
  Imb. Time%                                       --
  Calls                           0.002M/sec    690.0 calls
  CORE_TO_L2_CACHEABLE_REQUEST_ACCESS_STATUS:
    LS_RD_BLK_C                   0.402M/sec  161,687 req
  L2_PREFETCH_HIT_L2              0.150M/sec   60,526 hits
  L2_PREFETCH_HIT_L3              0.034M/sec   13,861 hits
  REQUESTS_TO_L2_GROUP1:L2_HW_PF  0.599M/sec  240,888 ops
  REQUESTS_TO_L2_GROUP1:RD_BLK_X  0.386M/sec  155,439 ops
  Cache Lines PF from OffCore     0.448M/sec  180,362 lines
  Cache Lines PF from Memory      0.414M/sec  166,501 lines
  Cache Lines Requested from
    Memory                        0.371M/sec  149,261 lines
  Write Memory Traffic GBytes     0.017G/sec     0.01 GB
  Read Memory Traffic GBytes      0.050G/sec     0.02 GB
  Memory traffic GBytes           0.067G/sec     0.03 GB
  Memory Traffic / Nominal Peak                  0.0%
  Average Time per Call                      0.000583 secs
  CrayPat Overhead : Time          0.1%
==============================================================================

...lots more information...
```
{: .output}

Finally, you can examine the results of any CrayPat run (including CrayPat-lite) graphically using Apprentice2.
For this you will need to log in with X11-forwarding enabled -- on Linux and macOS this means logging in with
the `-X` option passed to `ssh`. Then, with the `perftools-base` module loaded, we are able to examine the
results of our traced run as follows:
```
auser@ln01:~/nbody-par> app2 nbody-parallel.exe+apa+114297-5209t
```
{: .bash}


## Getting help with debugging and profiling tools

You can find more information on the debugging and profiling tools available on ARCHER2 in the ARCHER2 Documentation
and the Cray documentation:

* [ARCHER2 Documentation](https://docs.archer2.ac.uk)
* [Cray Technical Documentation](https://pubs.cray.com/)

If the documentation does not answer your questions then please contact
[the ARCHER2 Service Desk](https://www.archer2.ac.uk/support-access/servicedesk.html) and they
will be able to assist you.

{% include links.md %}

