---
title: "ARCHER2 development environment "
teaching: 30
exercises: 45
questions:
- "What does the ARCHER2 development environment look like and how do I access different components?"
- "How can I find out what compilers, tools, and libraries are available?"
- "How can I capture my current environment for reuse or to share with others?"
- "How can I get help with compiling and developing software on ARCHER2?"
objectives:
- "Know how to access different parts of the development environment on ARCHER2 using Lmod modules."
- "Know how to find out what is installed and where to get help."
keypoints:
- "The development environment is controlled through Lmod modules."
- "ARCHER2 supports the GCC and Cray compilers."
- "Compilers are accessed through the `ftn`, `cc` and `CC` wrapper commands."
- "The CSE service can help with software development issues."
---

## Using software modules on ARCHER2

ARCHER2 software modules use the [Lmod](https://lmod.readthedocs.io/en/latest/) system to provide
access to different software and versions on the system. The modules and versions available will
change across the lifetime of the service.

Software modules are provided by both HPE Cray and the ARCHER2 CSE team at
[EPCC](https://www.epcc.ed.ac.uk).

## What modules are loaded when you log into ARCHER2?

All users start with a default set of modules loaded into their environment. These include:

   - Cray Compiler Environment (CCE)
   - Cray MPICH MPI library
   - Cray LibSci scientific and numerical libraries
   - Cray lightweight performance analysis toolkit
   - System modules to enable use of the ARCHER2 hardware

You can see what modules you currently have loaded with the `module list` command:

```
auser@ln01:~> module list
```
{: .language-bash}
```
Currently Loaded Modules:
  1) cce/11.0.4                                  8) cray-mpich/8.1.4
  2) craype/2.7.6                                9) cray-libsci/21.04.1.1
  3) craype-x86-rome                            10) PrgEnv-cray/8.0.0
  4) libfabric/1.11.0.4.71                      11) bolt/0.7
  5) craype-network-ofi                         12) epcc-setup-env
  6) perftools-base/21.02.0                     13) load-epcc-module
  7) xpmem/2.2.40-7.0.1.0_2.7__g1d7a24d.shasta 
```
{: .output}

This means that your environment is already configured to build software using the Cray Compiling
Environment (CCE), the current default version of HPE Cray MPICH, and with numerical libraries provided
by the current default version of HPE Cray LibSci.

> ## Getting back if you purge or make a mistake
> 
> Unlike many other HPC systems you may have used, you should not generally use the
> `module purge` command before starting to use the system. Some of the modules
> loaded by default are required for you to be able to use the system correctly and
> so many things will not work if you use `module purge`. If you need to change the
> setup, you will generally use the `module load` or `module swap` commands instead.
>
> You should be able to return to the original module state by running `module restore`,
> but if you do find yourself with a broken environment you can usually fix things by
> logging out and logging back in again.
{: .callout}

## Finding out what software is available

You can query which modules are currently available to load with the `module avail` command:

```
auser@ln01:~> module avail
```
{: .language-bash}
```
------------------ /work/y07/shared/archer2-lmod/python/core -------------------
   matplotlib/3.4.3    netcdf4/1.5.7    seaborn/0.11.2

------------------- /work/y07/shared/archer2-lmod/libs/core --------------------
   adios/1.13.1                                 metis/5.1.0
   arpack-ng/3.8.0                              mkl/19.5-281
   boost/1.72.0                                 mkl/21.2-2883      (D)
   eigen/3.4.0                                  mumps/5.3.5
   epcc-cray-hdf5-parallel/1.12.0.3      (D)    parmetis/4.0.3
   epcc-cray-hdf5-parallel/1.12.0.7             petsc/3.14.2
   epcc-cray-netcdf-hdf5parallel/4.7.4.3 (D)    scotch/6.1.0
   epcc-cray-netcdf-hdf5parallel/4.7.4.7        slepc/3.14.1
   glm/0.9.9.6                                  superlu-dist/6.4.0
   hypre/2.18.0                                 superlu/5.2.2
   matio/1.5.18                                 trilinos/12.18.1

------------------- /work/y07/shared/archer2-lmod/apps/core --------------------
   castep/20.11                    nektar/5.0.3
   code_saturne/7.0.1-cce12        nwchem/7.0.2
   code_saturne/7.0.1-gcc11 (D)    onetep/6.1.3.7
   cp2k/cp2k-8.1                   openfoam/com/v2106
   elk/elk-7.2.42                  openfoam/org/v8.20200901
   gromacs/2021.3+plumed           openfoam/org/v9.20210903 (D)
   gromacs/2021.3           (D)    quantum-espresso/6.8
   lammps/29_Sep_2021              vasp/5/5.4.4.pl2-vtst
   namd/2.14-nosmp                 vasp/5/5.4.4.pl2         (D)
   namd/2.14                (D)    vasp/6/6.2.1

------------------- /work/y07/shared/archer2-lmod/utils/core -------------------
   bolt/0.7         (L)    genmaskcpu/1.0        visidata/2.1
   cdo/1.9.9rc1            ncl/6.6.2             vmd/1.9.3-gcc10
   cmake/3.21.3            other-software/1.1    xthi/1.2
   epcc-reframe/0.2        reframe/3.8.2
   epcc-setup-env   (L)    usage-analysis/1.1

--- /opt/cray/pe/lmod/modulefiles/mpi/crayclang/10.0/ofi/1.0/cray-mpich/8.0 ----
   cray-hdf5-parallel/1.12.0.3 (D)    cray-parallel-netcdf/1.12.1.3 (D)
   cray-hdf5-parallel/1.12.0.7        cray-parallel-netcdf/1.12.1.7

--------------- /opt/cray/pe/lmod/modulefiles/perftools/21.02.0 ----------------
   perftools                perftools-lite-gpu      perftools-preload
   perftools-lite           perftools-lite-hbm
   perftools-lite-events    perftools-lite-loops

--------- /opt/cray/pe/lmod/modulefiles/comnet/crayclang/10.0/ofi/1.0 ----------
   cray-mpich-abi/8.1.4 (D)    cray-mpich/8.1.4 (L,D)
   cray-mpich-abi/8.1.9        cray-mpich/8.1.9

------------------ /opt/cray/pe/lmod/modulefiles/net/ofi/1.0 -------------------
   cray-openshmemx/11.2.0 (D)    cray-openshmemx/11.3.3

---------------- /opt/cray/pe/lmod/modulefiles/cpu/x86-rome/1.0 ----------------
   cray-fftw/3.3.8.9 (D)    cray-fftw/3.3.8.11

------------ /opt/cray/pe/lmod/modulefiles/compiler/crayclang/10.0 -------------
   cray-hdf5/1.12.0.3 (D)    cray-hdf5/1.12.0.7

-------------------- /usr/share/lmod/lmod/modulefiles/Core ---------------------
   lmod    settarg

---------------------- /opt/cray/pe/lmod/modulefiles/core ----------------------
   PrgEnv-aocc/8.0.0     (D)      cray-libsci/21.08.1.2
   PrgEnv-aocc/8.1.0              cray-pals/1.0.17
   PrgEnv-cray/8.0.0     (L,D)    cray-pmi-lib/6.0.10    (D)
   PrgEnv-cray/8.1.0              cray-pmi-lib/6.0.13
   PrgEnv-gnu/8.0.0      (D)      cray-pmi/6.0.10        (D)
   PrgEnv-gnu/8.1.0               cray-pmi/6.0.13
   aocc/2.2.0                     cray-python/3.8.5.0    (D)
   aocc/2.2.0.1          (D)      cray-python/3.9.4.1
   aocc/3.0.0                     cray-stat/4.10.1       (D)
   atp/3.13.1            (D)      cray-stat/4.11.5
   atp/3.14.5                     craype/2.7.6           (L,D)
   cce/11.0.4            (L,D)    craype/2.7.10
   cce/12.0.3                     craypkg-gen/1.3.14     (D)
   cpe-cuda/21.09                 craypkg-gen/1.3.18
   cpe/21.04             (D)      gcc/9.3.0
   cpe/21.09                      gcc/10.2.0             (D)
   cray-R/4.0.3.0        (D)      gcc/10.3.0
   cray-R/4.1.1.0                 gcc/11.2.0
   cray-ccdb/4.11.1      (D)      gdb4hpc/4.12.5         (D)
   cray-ccdb/4.12.4               gdb4hpc/4.13.5
   cray-cti/2.13.6       (D)      iobuf/2.0.10
   cray-cti/2.15.5                papi/6.0.0.6           (D)
   cray-dsmml/0.1.4      (D)      papi/6.0.0.9
   cray-dsmml/0.2.1               perftools-base/21.02.0 (L,D)
   cray-jemalloc/5.1.0.4          perftools-base/21.09.0
   cray-libpals/1.0.17            valgrind4hpc/2.11.1    (D)
   cray-libsci/21.04.1.1 (L,D)    valgrind4hpc/2.12.4

------------- /opt/cray/pe/lmod/modulefiles/craype-targets/default -------------
   craype-accel-amd-gfx908    craype-hugepages256M    craype-network-none
   craype-accel-amd-gfx90a    craype-hugepages2G      craype-network-ofi  (L)
   craype-accel-host          craype-hugepages2M      craype-network-ucx
   craype-accel-nvidia70      craype-hugepages32M     craype-x86-milan
   craype-accel-nvidia80      craype-hugepages4M      craype-x86-rome     (L)
   craype-hugepages128M       craype-hugepages512M    craype-x86-trento
   craype-hugepages16M        craype-hugepages64M
   craype-hugepages1G         craype-hugepages8M

------------------------- /usr/local/share/modulefiles -------------------------
   load-epcc-module (L)

---------------------------- /opt/cray/modulefiles -----------------------------
   cray-lustre-client/2.12.4.2_cray_63_g79cd827-7.0.1.0_8.1__g79cd827237.shasta
   cray-shasta-mlnx-firmware/1.0.8
   dvs/2.12_4.0.112-7.0.1.0_15.1__ga97f35d9
   libfabric/1.11.0.4.71                                                        (L)
   xpmem/2.2.40-7.0.1.0_2.7__g1d7a24d.shasta                                    (L)

------------------------------- /opt/modulefiles -------------------------------
   aocc/2.2.0      aocc/3.0.0        gcc/8.1.0    gcc/10.2.0
   aocc/2.2.0.1    cray-R/4.0.3.0    gcc/9.3.0

Use "module spider" to find all possible modules and extensions.
Use "module keyword key1 key2 ..." to search for all possible modules matching
any of the "keys".

```
{: .output}

The output lists the available modules and their versions. It also shows you which modules are
loaded by default (marked with `(D)`) when there are multiple versions available and you do
not specify the version when you load. Modules marked with `(L)` are currently loaded.

> ## Licensed software
> Some of the software installed on ARCHER2 requires the user to have their licence validated before they
> can use it on the service. More information on gaining access to licensed software through the SAFE
> is provided below.
{: .callout}

If you want more information on a particular module, you can use the `module help` and
`module whatis` commands. For example, to get more info on the `cray-mpich` module:

```
auser@ln01:~> module help cray-mpich
```
{: .language-bash}
```
----------------- Module Specific Help for "cray-mpich/8.1.9" ------------------


                                  COPYRIGHT

The following is a notice of limited availability of the code, and disclaimer
which must be included in the prologue of the code and in all source listings
of the code.

Copyright Notice
 + 2002 University of Chicago

Permission is hereby granted to use, reproduce, prepare derivative works, and
to redistribute to others.  This software was authored by:

Mathematics and Computer Science Division
Argonne National Laboratory, Argonne IL 60439

(and)

Department of Computer Science
University of Illinois at Urbana-Champaign


                              GOVERNMENT LICENSE

Portions of this material resulted from work developed under a U.S.
Government Contract and are subject to the following license: the Government
is granted for itself and others acting on its behalf a paid-up, nonexclusive,
irrevocable worldwide license in this computer software to reproduce, prepare
derivative works, and perform publicly and display publicly.

                                  DISCLAIMER

This computer code material was prepared, in part, as an account of work
sponsored by an agency of the United States Government.  Neither the United
States, nor the University of Chicago, nor any of their employees, makes any
warranty express or implied, or assumes any legal liability or responsibility
for the accuracy, completeness, or usefulness of any information, apparatus,
product, or process disclosed, or represents that its use would not infringe
privately owned rights.


Cray MPICH 8.1.9:
=======================================

Release Date:
-------------
  August 18, 2021


Purpose:
--------
  Cray MPICH 8.1.9 is based upon ANL MPICH 3.4a2 with support for libfabric
  and is optimized for the Cray Programming Environment.

  Major Differences Cray MPICH 8.1.9 from the XC Cray MPICH include:

      - Uses the new ANL MPICH CH4 code path and libfabric for network
        support.

      - Does not support -default64 mode for Fortran

      - Does not support C++ language bindings

...lots more information... 
```
{: .output}

## Loading and switching modules

Lets look at our environment before we change anything. To see just our loaded modules we use the
`module list` command:

```
auser@ln01:~> module list
```
{: .language-bash}
```
Currently Loaded Modules:
  1) cce/11.0.4                                  8) cray-mpich/8.1.4
  2) craype/2.7.6                                9) cray-libsci/21.04.1.1
  3) craype-x86-rome                            10) PrgEnv-cray/8.0.0
  4) libfabric/1.11.0.4.71                      11) bolt/0.7
  5) craype-network-ofi                         12) epcc-setup-env
  6) perftools-base/21.02.0                     13) load-epcc-module
  7) xpmem/2.2.40-7.0.1.0_2.7__g1d7a24d.shasta
```
{: .output}

You load modules with the `module load` command. For example, to load the default `cray-fftw` module:

```
auser@ln01:~> module load cray-fftw
```
{: .language-bash}

Now, lets list our loaded modules again with `module list`:

```
auser@ln01:~> module list
```
{: .language-bash}
```
Currently Loaded Modules:
  1) cce/11.0.4                                  8) cray-mpich/8.1.4
  2) craype/2.7.6                                9) cray-libsci/21.04.1.1
  3) craype-x86-rome                            10) PrgEnv-cray/8.0.0
  4) libfabric/1.11.0.4.71                      11) bolt/0.7
  5) craype-network-ofi                         12) epcc-setup-env
  6) perftools-base/21.02.0                     13) load-epcc-module
  7) xpmem/2.2.40-7.0.1.0_2.7__g1d7a24d.shasta  14) cray-fftw/3.3.8.9
```
{: .output}

You can see that we now have the default `cray-fftw` module `cray-fftw/3.3.8.9` as we did not specify a version
explicitly.

<!--
  Note: Once Lmod is working, we likely want to choose a module that has dependencies so we 
  can note that Lmod has automatically loaded dependencies. -->

To unload a module, you can handily run the `module unload` command. Since you will generally
only have one version of a given software package loaded at one time, you don't need to specify the
version to be unloaded:

```
auser@ln01:~> module unload cray-fftw
```
{: .language-bash}

If you want to swap two versions of the same module then you use the `module swap` command, which
combines the `module unload` and `module load` commands:

```
auser@ln01:~> module swap cray-fftw cray-fftw/3.3.8.11
```
{: .language-bash}
```
The following have been reloaded with a version change:
  1) cray-fftw/3.3.8.9 => cray-fftw/3.3.8.11
```
{: .output}

> ## The convenience command `ml`
> The command `ml` acts as a shortcut for the `module list` and `module load` commands which are
> typically used very often. It has a different meaning depending on whether or not you
> supply an argument. Without an argument it means `module list`, and with an argument it means
> `module load`.
{: .callout}

## Hierarchical modules

Lmod uses hierarchical modules. This means that modules which depend on others will only become
available to load once the dependencies have been loaded. The reason for this is that some libraries
and applications will be installed multiple times with different dependencies. For example,
`cray-netcdf` has multiple installs, each relying on a different version of `cray-hdf5`. Once we
load a `cray-hdf5` module, Lmod will be able to load a correct `cray-netcdf` module.

Let's try this out. Initially, we won't be able to locate `cray-netcdf` using
`module avail`, for the reasons given above. To search for it we should instead use
`module spider`, which searches for and shows all modules matching the pattern we give, even
if they're currently unavailable.

```
auser@ln01:~> module spider cray-netcdf
```
{: .language-bash}
```
----------------------------------------------------------------------------
  cray-netcdf:
----------------------------------------------------------------------------
     Versions:
        cray-netcdf/4.7.4.3
        cray-netcdf/4.7.4.7
     Other possible modules matches:
        cray-netcdf-hdf5parallel  epcc-cray-netcdf-hdf5parallel

----------------------------------------------------------------------------
  To find other possible module matches execute:

      $ module -r spider '.*cray-netcdf.*'

----------------------------------------------------------------------------
  For detailed information about a specific "cray-netcdf" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider cray-netcdf/4.7.4.7
----------------------------------------------------------------------------
```
{: .output}

We can see that two versions are available, `4.7.4.3` and `4.7.4.7`. Let's find out how to
load `cray-netcdf/4.7.4.3`, providing the full module name as the argument to
`module spider`:

```
auser@ln01:~> module spider cray-netcdf/4.7.4.3
```
{: .language-bash}
```
----------------------------------------------------------------------------
  cray-netcdf: cray-netcdf/4.7.4.3
----------------------------------------------------------------------------

    You will need to load all module(s) on any one of the lines below before the "cray-netcdf/4.7.4.3" module is available to load.

      aocc/2.2.0  cray-hdf5/1.12.0.3
      aocc/2.2.0.1  cray-hdf5/1.12.0.3
      cce/11.0.4  cray-hdf5/1.12.0.3
      cce/11.0.4  cray-hdf5/1.12.0.7
      gcc/10.2.0  cray-hdf5/1.12.0.3
      gcc/10.2.0  cray-hdf5/1.12.0.7
      gcc/10.3.0  cray-hdf5/1.12.0.3
      gcc/10.3.0  cray-hdf5/1.12.0.7
      gcc/11.2.0  cray-hdf5/1.12.0.3
      gcc/11.2.0  cray-hdf5/1.12.0.7
      gcc/9.3.0  cray-hdf5/1.12.0.3
      gcc/9.3.0  cray-hdf5/1.12.0.7

    Help:
      Release info:  /opt/cray/pe/netcdf/4.7.4.3/release_info
```
{: .output}

This output is telling us what we have to load before we can in turn load
`cray-netcdf/4.7.4.3`. We will already have a compiler loaded,
so what we're missing is a version of `cray-hdf5`. If we choose and load one, we
will then find `cray-netcdf` listed by `module avail` and be able to load it:

```
auser@ln01:~> module load cray-hdf5/1.12.0.3
auser@ln01:~> module avail cray-netcdf
```
{: .language-bash}
```
------ /opt/cray/pe/lmod/modulefiles/hdf5/crayclang/10.0/cray-hdf5/1.12.0 ------
   cray-netcdf/4.7.4.3 (D)    cray-netcdf/4.7.4.7

------------------- /work/y07/shared/archer2-lmod/libs/core --------------------
   epcc-cray-netcdf-hdf5parallel/4.7.4.3 (D)
   epcc-cray-netcdf-hdf5parallel/4.7.4.7

  Where:
   D:  Default Module

Use "module spider" to find all possible modules and extensions.
Use "module keyword key1 key2 ..." to search for all possible modules matching
any of the "keys".
```
{: .output}
```
auser@ln01:~> module load cray-netcdf/4.7.4.3
auser@ln01:~> module list
```
{: .language-bash}
```
Currently Loaded Modules:
  1) cce/11.0.4                                  9) cray-libsci/21.04.1.1
  2) craype/2.7.6                               10) PrgEnv-cray/8.0.0
  3) craype-x86-rome                            11) bolt/0.7
  4) libfabric/1.11.0.4.71                      12) epcc-setup-env
  5) craype-network-ofi                         13) load-epcc-module
  6) perftools-base/21.02.0                     14) cray-hdf5/1.12.0.3
  7) xpmem/2.2.40-7.0.1.0_2.7__g1d7a24d.shasta  15) cray-netcdf/4.7.4.3
  8) cray-mpich/8.1.4
```
{: .output}

> ## Finding modules
> If you can't see a module for software you expect to be installed when using `module avail`,
> try searching for it instead with `module spider`.
{: .callout}

## Switching between programming environments

You use the `module load` command to switch between different programming environments available
on ARCHER2, each corresponding to a different compiler suite. When you load a library, Lmod will
load the installation meant for use with the current programming environment. The available
environments are:

* Cray Compiling Environment (CCE): `PrgEnv-cray`
* GNU Compiler Collection (GCC): `PrgEnv-gnu`
* AMD Optimizing Compilers (AOCC): `PrgEnv-aocc`

`PrgEnv-cray` is loaded by default when you log in. You can change to the GCC compiler environment with the
`module load PrgEnv-gnu` command (with `module list` afterwards to show how things have changed):

```
auser@ln01:~> module load PrgEnv-gnu
```
{: .language-bash}
```
Lmod is automatically replacing "cce/11.0.4" with "gcc/10.2.0".

Lmod is automatically replacing "PrgEnv-cray/8.0.0" with "PrgEnv-gnu/8.0.0".

Due to MODULEPATH changes, the following have been reloaded:
  1) cray-mpich/8.1.4
```
{: .output}
```
auser@ln01:~> module list
```
{: .language-bash}
```
Currently Loaded Modules:
  1) gcc/10.2.0                                  8) cray-mpich/8.1.4
  2) craype/2.7.6                                9) cray-libsci/21.04.1.1
  3) craype-x86-rome                            10) bolt/0.7
  4) libfabric/1.11.0.4.71                      11) epcc-setup-env
  5) craype-network-ofi                         12) load-epcc-module
  6) perftools-base/21.02.0                     13) PrgEnv-gnu/8.0.0
  7) xpmem/2.2.40-7.0.1.0_2.7__g1d7a24d.shasta
```
{: .output}

## Switching between different compilers and CPE versions

You can use the `module swap` command to switch between different
versions of the compilers within the programming environments. To do this, you need to know the
names of the compiler modules. These are:

* `cce` for the Cray Compiling Environment
* `gcc` for the GNU Compiler Collection
* `aocc` for the AMD Optimising Compilers

For example, to change the version of GCC you are using you would first load the GNU programming
environment and then switch to a different version of GCC:

```
auser@ln01:~> module load PrgEnv-gnu
auser@ln01:~> module swap gcc gcc/11.2.0
```
{: .language-bash}

Different versions of the compilers are part of coherent Cray Programming Environment (CPE)
releases, each containing the full set of Cray software that is guaranteed to work together
as a whole. As a different version of the compiler is not guaranteed to work, it is generally
preferable to instead change the full CPE version.

The default software on logging in to ARCHER2 is provided by CPE 21.04. A newer version,
CPE 21.09, is also available, and can be loaded by running:

```
auser@ln01:~> module load cpe/21.09
```
{: .language-bash}

This will swap the compiler to the one from the new CPE version and also change the default
version of the HPE Cray modules to those which are part of the CPE. At this point you can load
the programming environment you need followed by any libraries. Finally, you should modify
`LD_LIBRARY_PATH` to include the new non-default libraries:

```
auser@ln01:~> export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
```
{: .language-bash}

## Compiler wrappers

Code compiled on ARCHER2 should use the compiler wrapper scripts to ensure that
all the correct libraries and header files are included in both the compile and
link stages.

> ## Linked libraries
> Loading a module will automatically set the necessary paths and link flags for
> that software, eliminating the need to specify them manually using the
> `LDFLAGS` variable or another mechanism. For example, the `cray-libsci` module
> is loaded by default on login and sets the the wrappers to link in BLAS,
> LAPACK, and ScaLAPACK functionality without the need for the user to provide
> an explicit `-llib_sci`.
{: .callout}

The compiler wrapper scripts are available for Fortran, C, and C++:
* `ftn` -- Fortran compiler
* `cc` -- C compiler
* `CC` -- C++ compiler

The wrapper scripts can be used to compile both sequential and parallel codes.

In practice, this means that your build scripts should always use the compiler
wrappers rather than directly using e.g. `gcc`, `gfortran` or the like. This
may require you to edit your scripts or set environment variables, for example
the `CC` and `MPICC` environment variables:

```
auser@ln01:~> export CC=cc
auser@ln01:~> export MPICC=cc
```
{: .language-bash}

`man` pages are available for the compiler wrappers. This information can be read
no matter which programming environment you are in. Programming environment specific
`man` pages can be read for the CCE and GCC compilers by running with the name
of the desired compilers: so, `craycc` and `crayftn` for CCE, and `gcc` and `gfortran`
for GCC.

## Dynamic linking on ARCHER2

ARCHER2 currently only supports dynamic linking (you will see an error if you 
try to link applications statically). As well as linking dynamically, all 
HPE Cray software libraries are added to the executable *RUNPATH*. This means
that you do not need to load programming environments or particular
Cray modules in your job submission scripts to run on the compute nodes. However,
if you do load different versions of libraries used by your application in your
job submission scripts, these will take precedence over the versions encoded in
the binary RUNPATH.

**Note:** Software libraries provided by the CSE team via modules **are not**
captured in the binary RUNPATH so these **do** need to be loaded in your 
job submission scripts for the dynamic executable to be able to find them
successfully.

## Getting help with software

You can find more information on the software available on ARCHER2 in the ARCHER2 Documentation at:

* [ARCHER2 Documentation](https://docs.archer2.ac.uk)

This includes information on the software provided by Cray and the software provided by the 
ARCHER2 CSE Service at EPCC.

If the software you require is not currently available or you are having trouble with the installed
software then please contact
[the ARCHER2 Service Desk](https://www.archer2.ac.uk/support-access/servicedesk.html) and they
will be able to assist you.

{% include links.md %}

