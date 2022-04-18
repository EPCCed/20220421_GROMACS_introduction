---
title: "Running GROMACS simulations on ARCHER2"
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

## What do we mean by classical molecular dynamics


## The GROMACS molecular dynamics parameters (MDP) file

[GROMACS molecular dynamics parameter file](https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html)


## Running a simulations

[GROMACS `mdrun` page](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-mdrun.html)

## Energy minimisation

You'll need to copy your `5pep-neutral.gro` and `topol.top` files.

```bash
  gmx grommp -f minim.mdp -c 5pep-neutral.gro -p topol.top -o em.tpr
```

```bash
  gmx mdrun -v -deffnm em
```

```bash
  gmx energy -f em.edr -o potential.xvg
```

## NPT simulations

```bash
  ; Intergrator, timestep, and total run time
  integrator               = md
  dt                       = 0.002
  nsteps                   = 500000
  
  ; Logs and outputs
  nstlog                   = 5000
  nstenergy                = 5000
  
  ; Bond constraints
  constraints              = h-bonds
  constraint-algorithm     = lincs
  
  ; Van der Waals interactions
  vdwtype                  = Cut-off
  rvdw                     = 1.0
  cutoff-scheme            = Verlet
  DispCorr                 = EnerPres
  
  ; Coulombic interactions
  coulombtype              = PME
  rcoulomb                 = 1.0
  
  ; Thermostat
  tcoupl                   = V-rescale
  tc-grps                  = Protein  SOL NA
  ref-t                    = 300      300 300
  tau-t                    = 0.1      0.1 0.1
  
  ; Barostat
  pcoupl                   = Parrinello-Rahman
  ref-p                    = 1.0
  tau-p                    = 2.0
  compressibility          = 4.5e-5

```

```bash
  gmx grompp -f npt.mdp -c em.gro -r em.gro -p topol.top -o npt.tpr
```


{% include links.md %}

