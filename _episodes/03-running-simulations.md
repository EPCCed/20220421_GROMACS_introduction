---
title: "Running GROMACS simulations on ARCHER2"
teaching: 25
exercises: 20
questions:
- ""
objectives:
- ""
keypoints:
- ""
---

## What do we mean by classical molecular dynamics



## The GROMACS molecular dynamics parameters (MDP) file

[GROMACS molecular dynamics parameter file](https://manual.gromacs.org/documentation/current/user-guide/mdp-options.html)


## Running a simulations

[GROMACS `mdrun` page](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-mdrun.html)

## Energy minimisation

For this session, you will need a copy of the `5pep-neutral.gro` and 
`topol.top` files generated in the previous session. You can either copy these 
across or use pre-generated ones. To get the pregenerated files, run:

```bash
  svn checkout https://github.com/EPCCed/20220421_GROMACS_introduction/trunk/exercises
```

```bash
  ; minim.mdp - used as input into grompp to generate em.tpr
  ; Parameters describing what to do, when to stop and what to save
  integrator      = steep         ; Algorithm (steep = steepest descent minimization)
  emtol           = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
  emstep          = 0.01          ; Minimization step size
  nsteps          = 50000         ; Maximum number of (minimization) steps to perform

  ; Logs and outputs
  nstlog          = 500
  nstenergy       = 500
  nstxout         = 500
  
  ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
  nstlist         = 1         ; Frequency to update the neighbor list and long range forces
  cutoff-scheme   = Verlet    ; Buffered neighbor searching
  ns_type         = grid      ; Method to determine neighbor list (simple, grid)
  coulombtype     = PME       ; Treatment of long range electrostatic interactions
  rcoulomb        = 1.0       ; Short-range electrostatic cut-off
  rvdw            = 1.0       ; Short-range Van der Waals cut-off
  pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
```

```bash
  gmx grompp -f minim.mdp -c 5pep-neutral.gro -p topol.top -o ener_minim.tpr
```

```bash
  gmx mdrun -v -s ener_minim.tpr
```

```bash
  gmx energy -f ener.edr -o potential.xvg
```

Select option `10` to plot the total potential energy of the system and 
option `0` to exit the program.

## NPT simulations

```bash
  ; Intergrator, timestep, and total run time
  integrator               = md
  dt                       = 0.002
  nsteps                   = 500000
  
  ; Logs and outputs
  nstlog                   = 500
  nstenergy                = 500
  nstxout                  = 500
  
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

