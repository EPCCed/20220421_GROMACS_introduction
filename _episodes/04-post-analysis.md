---
title: "Visualising and analysing simulation data"
teaching: 25
exercises: 20
questions:
- "What output files does GROMACS produce?"
- "What properties of my system can I see after my simulation has run?"
- "How can I visualise my system?"
objectives:
- "Learn how to extract the thermodynamic properties from the output"
- "Learn about the types of analysis we can perform"

keypoints:
- "Use of gmx energy, make_ndx, rms"
- "System visualisation in vmd"
---

Post-processing and analysis tools
==================================

With the simulation complete, we can analyse the simulation trajectory and 
understand what the simulation has demonstrated. GROMACS offers a number of 
post-simulation analysis tools. In this lesson, we will discuss tools that 
can be used to: generate the thermodynamic properties of interest; obtain 
radial distribution functions and correlation functions; 

Thermodynamic properties of the system
--------------------------------------

The GROMACS ``energy`` tool can be used to extract energy components from an 
energy (``.edr``) file. By default, this tool will generate an XMGrace file. 
To use this, run:

```
  gmx energy -f ener.edr -o temp.xvg
```

When running this, you will get a prompt asking which property you would like 
output (*e.g.* potential energy, kinetic energy, pressure, temperature, 
*etc.*). 

```
Select the terms you want from the following list by
selecting either (part of) the name or the number or a combination.
End your selection with an empty line or a zero.
-------------------------------------------------------------------
  1  Bond             2  Angle            3  Proper-Dih.      4  Improper-Dih. 
  5  LJ-14            6  Coulomb-14       7  LJ-(SR)          8  Disper.-corr. 
  9  Coulomb-(SR)    10  Coul.-recip.    11  Potential       12  Kinetic-En.   
 13  Total-Energy    14  Conserved-En.   15  Temperature     16  Pres.-DC      
 17  Pressure        18  Constr.-rmsd    19  Box-X           20  Box-Y         
 21  Box-Z           22  Volume          23  Density         24  pV            
 25  Enthalpy        26  Vir-XX          27  Vir-XY          28  Vir-XZ        
 29  Vir-YX          30  Vir-YY          31  Vir-YZ          32  Vir-ZX        
 33  Vir-ZY          34  Vir-ZZ          35  Pres-XX         36  Pres-XY       
 37  Pres-XZ         38  Pres-YX         39  Pres-YY         40  Pres-YZ       
 41  Pres-ZX         42  Pres-ZY         43  Pres-ZZ         44  #Surf*SurfTen 
 45  Box-Vel-XX      46  Box-Vel-YY      47  Box-Vel-ZZ      48  T-Protein     
 49  T-SOL           50  T-NA            51  Lamb-Protein    52  Lamb-SOL      
 53  Lamb-NA       
```


Enter the `15` to generate an XMGrace file that includes the variation
of the temperature over time. You will need to press Enter twice.


The data can be plotted with xmgrace or other plotting software.
As you can see there are a number of other options for the ``energy`` command,
and these can be found in the GROMACS manual 
[gmx energy](http://manual.gromacs.org/documentation/current/onlinehelp/gmx-energy.html)
page.

Visualising with VMD
--------------------


Now we are going to look at the simulation in vmd. This is installed on 
ARCHER2. You can use it with:

```
module load vmd
vmd
```
It may need a while to open. But you should eventually see the main window.
To load our trajectory:

* Go to File -> New molecule
* Browse and select ``confout.gro``
* Then "Load"

The image of our protein system should appear in the main window.

Then to load the MD trajectory:

* Make sure the following is selected,
   Load files for: 0: confout.gro
* Browse and select ``traj_conf.xtc``
* Then "Load"

The main window will now show the motion of MD of the system.

This can also be done quickly by using the command

```
vmd confout.gro traj_comp.xtc
```


Generating an index file
------------------------


```
  gmx make_ndx -f confout.gro -o 5pep.ndx
```

This takes the coordinates used in the mdrun. It will then
analyse the system, and output the default index groups. 

You should see something like the following:

```

  0 System              : 62149 atoms
  1 Protein             :  4682 atoms
  2 Protein-H           :  2426 atoms
  3 C-alpha             :   326 atoms
  4 Backbone            :   978 atoms
  5 MainChain           :  1305 atoms
  6 MainChain+Cb        :  1596 atoms
  7 MainChain+H         :  1618 atoms
  8 SideChain           :  3064 atoms
  9 SideChain-H         :  1121 atoms
 10 Prot-Masses         :  4682 atoms
 11 non-Protein         : 57467 atoms
 12 Water               : 57429 atoms
 13 SOL                 : 57429 atoms
 14 non-Water           :  4720 atoms
 15 Ion                 :    38 atoms
 16 NA                  :    38 atoms
 17 Water_and_ions      : 57467 atoms
 
  nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index
 ```

This lists all the default groups generated from the coordinate file.

It is possible to create new index groups by using the command prompts listed.
For example you can create a group of the O atoms in water with:

```
> a OW

Found 19143 atoms with name OW
```


>> ### Exercise
>> 
>> Create a group of containing the hydrogen atoms in the water molecules.
>>> ### Solution
>>> ```
>>> a HW1 | a HW2
>>> q
>>> ```
>>> The OR | command can be used to select multiple groups of atoms.
> {: .solution}
{: .challenge}



If you now open the ``5pep.ndx`` file we just generated you should see
the different groups of atoms listed. Each number is an atom number from the 
``confout.gro file``

Index files are useful if you wish to do analysis on specific groups of atoms.
This can be used in calculating the RDF, order parameters, densities, 
displacements etc. The index file can be used in a command with the ``-n``
option.

For more information on ``make_ndx``, please see the
GROMACS manual
[gmx make_ndx](http://manual.gromacs.org/documentation/current/onlinehelp/gmx-make_ndx.html) 
page.

Root mean squared deviation
----------------------------

The following can be run to calculate the root mean squared deviation and
least-squares fit for a group of atoms:

```
 gmx rms -s npt.tpr -f traj.trr -n 5pep.ndx -o rmsd.xvg -tu ns
```

You will be prompted for a group number to do the calculation for.
Select option ``4`` to do the Backbone of the protein. The ``-tu`` option
gives the time unit.


Continuing your simulation
---------------------------

You may wish to restart your simulation from the point you left off
in order to progress the simulation further or if your simulation 
does not complete. GROMACS does checkpointing during a simulation
to allow this. 

A simulation can be restarted with:

```
gmx_mpi mdrun -cpi state
```

providing you have a ``state.cpt`` checkpoint file.

In our case our simulation ran to completion, so we cannot restart it,
however we may extend it. This can be done by generating a new tpr
with the ``convert-tpr`` command:

```
gmx convert-tpr -s ${OLD}.tpr -extend timetoextendby -o ${NEW}.tpr
```

where the timetoextendby is the time in ps

The job can then be run with the following gmx_mpi command:

```
gmx_mpi mdrun -s ${NEW}.tpr -cpi state.cpt
```


>> ### Exercise
>> 
>> Extend our simulation by 100ps.
>>> ### Solution
>>> ```
>>> gmx convert-tpr -s npt.tpr -extend 100 -o npt-new.tpr
>>> gmx_mpi mdrun -s npt-new.tpr -cpi state.cpt
>>> ```
> {: .solution}
{: .challenge}


Other analysis tools
-----------------------

A guide for analysing your trajectory in GROMACS can be found on the 
[website](https://manual.gromacs.org/documentation/2019/reference-manual/analysis.html)


{% include links.md %}


