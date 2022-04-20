---
title: "Preparing a GROMACS system"
teaching: 30
exercises: 15
questions:
- "How do I set up a system ready to be run with GROMACS"
objectives:
- "Learn how to prepare your system prior to running a simulation"
keypoints:
- "Creating Gromacs topology files"
- "Creating your simulation box"
- "Solvating a system"
- "Adding ions to neurtalise your system"
---




# Preparing a system for GROMACS


In this lesson we will describe how to choose and
prepare the input files for a system before to running a Gromacs 
simulation.

As an example we will be looking at pepsin. This is an enzyme used in
digestion. The PDB entry can be 
found [here](https://www.rcsb.org/structure/5pep) and the PDB file can be
downloaded from the protein data bank
using:

```
wget https://files.rcsb.org/download/5pep.pdb
```


This protein is used as an example but this process could apply to a general
system where
you start from a `.pdb` file and end up with a set of input files, ready to
run a GROMACS simulation.

### Splitting up the system

A pdb file downloaded from the database may not be in the format you want
for running your simulation. For example pdb files can contain non-protein
residues and waters. In this case we have some waters already in the pdb
file that we want to remove as we will be solvating the system later on.
You can see these in the pdb file if you open it. 

To just get the protein itself we can use the following grep command to remove
lines containing the water symbol ``'HOH'``.

```
grep -v 'HOH' 5pep.pdb > 5pep_protein.pdb
```

The new ``5pep_protein.pdb`` file contains just the protein itself.

### Creating a Gromacs topology (PDB2GMX)


Now we can create a Gromacs topolopy for the system.
The GROMACS ``pdb2gmx`` command is used to convert a `pdb` coordinate file into a 
 GROMACS topology file and also create a processed structure file 
in the GROMACS format ``.gro``. 

First you will need to load the GROMACS module on ARCHER2:

```
module load gromacs
```
Then you can run:

```
gmx pdb2gmx -f 5pep_protein.pdb
```

You will be prompted to select the forcefield you would like to use. GROMACS comes with 
a number of AMBER and GROMOS forcefields, as well as a CHARMM and an OPLS-AA
option. You will also need to specify your water model (choices included are 
TIP models, and the SPC and SPC/E models). 

---
**NOTE**

We do not currently  have any waters in our system, however the water model 
is needed to ensure the topology is consistent when the system is solvated.

---

```
Command line:
  gmx pdb2gmx -f 5PEP_protein.pdb

Select the Force Field:

From '/work/y07/shared/apps/core/gromacs/2021.3/share/gromacs/top':

 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)

 2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)

 3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)

 4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)

 5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)

 6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)

 7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)

 8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)

 9: GROMOS96 43a1 force field

10: GROMOS96 43a2 force field (improved alkane dihedrals)

11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)

12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)

13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)

14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)

15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)
```

Select the AMBER03 protein (option 1), followed by the TIP3P water model (option 1).

```
ls
5pep.pdb  5pep_protein.pdb  conf.gro  posre.itp  topol.top
```

The code above 
produces three outputs: a system topology ``topol.top``, a 
position restraint file ``posre.itp`` (included in the topology file), and a coordinate file ``conf.gro``. 
Further to these files, ``pdb2gmx`` will output a number of interesting 
details to screen, such as the total mass of the system given the coordinates 
and topology being used as well as the net charge of the system.

```
Total charge -38.000 e
```

The  total charge 
is  important to note and will be used in the `Solvating and 
ionise a system` step of system preparation.

More information about PDB2GMX can be found in 
the GROMACS 
[PDB2GMX manual](http://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html)

#### The topol.top file

This is the topology file for the system. You can open it with a text editor
of your choice.

```
vi topol.top
```

Commented lines begin with a semicolon ;
You will see a number of these at the top of the file. Followed by the following:

```
; Include forcefield parameters
#include "amber03.ff/forcefield.itp"

[ moleculetype ]
; Name            nrexcl
Protein_chain_A     3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
; residue   1 ILE rtp NILE q +1.0
     1         N3      1    ILE      N      1  -0.605757      14.01
     2          H      1    ILE     H1      2   0.451919      1.008
```

The #include statemate includes the amber03 forcefield that we selected
when running pdb2gmx. 

The [ moleculetype ] gives the molecule name and the number of exclusions. The
number of exclusions is the number of bonds away for which non-bonded interactions
between pairs of atoms are excluded.
The [ atoms ] directive lists the atoms in the molecule, grouped by the different
residues. The meaning of the columns are as follows:

* ``nr``: atom number
* ``type``: atom type
* ``resnr``: residue number
* ``atom``: atom name (these are unique to the residue)
* ``cgnr``: charge group number
* ``charge``: the atomic charge
* ``mass``: the atomic mass
* ``typeB``, ``chargeB``, ``massB``: Additional type, charge and mass used for free energy perturbation 

The bonds directive lists pairs of atoms which are bonded:
```
[ bonds ]       
;  ai    aj funct            c0            c1            c2            c3
    1     2     1      
    1     3     1      
```    
 The pairs directive lists LJ pairs of atoms by their atom number:
 ```
 [ pairs ]
;  ai    aj funct            c0            c1            c2            c3
    1     8     1
    1     9     1
 ```   
 The angles directive lists triplets of 3 atoms for angle parameterisation:
 ```
 [ angles ]
;  ai    aj    ak funct            c0            c1            c2            c3
    2     1     3     1
    2     1     4     1
```
The dihedreal angles for 4 atoms:
```
[ dihedrals ]
;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5
    2     1     5     6     9
    2     1     5     7     9
```

These `func` numbers give the function type to be used for the interaction. 
This corresponds to interactions listed in the forcefield field ``.itp`` for 
different atom pair, triplets etc.


The final lines of the file are:
```
; Include Position restraint file
#ifdef POSRES  
#include "posre.itp" 
#endif   
   
; Include water topology    
#include "amber03.ff/tip3p.itp"
   
#ifdef POSRES_WATER 
; Position restraint for each water oxygen
[ position_restraints ]     
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif  
  
; Include topology for ions 
#include "amber03.ff/ions.itp"
  
[ system ]
; Name  
PEPSIN  

[ molecules ] 
; Compound        #mols     
Protein_chain_A     1
```

* This includes the position restriant file (which was created during 
``pdb2gmx``) contains posisition restraits for the atoms. The if 
statement ensures it is only applied when POSRES is true.

* The water forcefield tip3p is included

* The position restraints for waters

* The forcefield for ions

* The name of the system

* The molecules in the system (so far just one protein chain)

#### The conf.gro file

This contains the atomic coordinates of the system

```
PEPSIN
 4682
    1ILE      N    1  -0.269   5.104   9.782
    1ILE     H1    2  -0.361   5.123   9.816
    1ILE     H2    3  -0.209   5.080   9.858
    1ILE     H3    4  -0.234   5.185   9.735
    1ILE     CA    5  -0.275   4.991   9.687
    1ILE     HA    6  -0.335   5.019   9.612
    ...
```
The top lines are the system name and the number of atoms.

The atomic coordinates in x,y,z of each atom are then listed in nm. These
are preceeded by the residue number, residue name, atom name
and atom number.

At the very end of the file is the box dimensions listed as:

```
box[X][X],box[Y][Y],box[Z][Z], box[X][Y],box[X][Z],box[Y][X],box[Y][Z],box[Z][X],box[Z][Y]
```

#### The posre.itp file

This file contains the position restraints for the system. The force
constants in x,y,z are listed for the atom numbers.

More information about how these are defined is given [here](https://manual.gromacs.org/documentation/current/reference-manual/functions/restraints.html#positionrestraint)

More information about file formats can be found in the GROMACS manual
[file formats](http://manual.gromacs.org/documentation/current/reference-manual/file-formats.html).

### Generating your own forcefield file

Occasionally a system will contain non-protein residues which need separate 
parameterisation as there are not forcefields available for them. You will 
see a message when running ``pdb2gmx`` such as:

```
Residue 'XXX' not found in residue topology database
```

In this case you must generate your own forcefield. In GROMACS to do this
you can generate your own forcefield files in a {forcefield}.ff directory.
More about this can be found on the Gromacs manual.

You can also generate an Amber or CHARMM topology by using the   AmberTool 
 ``antechamber`` function or the CHARMM ``cgenff`` function. To do this, you 
should follow the procedures described above, making sure to select an 
appropriate forcefield from the selection GROMACS provides. Then, use a 
parameter-generating tool like ``antechamber`` with ``actype`` (for Amber) 
or ``cgenff`` (for CHARMM). The topologies generated in this way can then be 
added to the GROMACS topology that you generated. 
  

---
**NOTE**

It is advised to use one of the pre-existing GROMACS forcefields where 
possible. Only consider generating your own forcefield if the ones 
available on GROMACS do not fulfill your requirements.

---

For more information on generating your own forcefield, please see the GROMACS
manual pages about 
[adding a residue](http://manual.gromacs.org/documentation/current/how-to/topology.html)
and [force field organisations](http://manual.gromacs.org/documentation/current/reference-manual/topologies/force-field-organization.html)



## Preparing and solvating your simulation box


### Generating a simulation box


Now that a topology has been generated, the next step is to generate a 
simulation box into which to place this topology. For this, use the 
``editconf`` command. This tool has a number of functionalities, including 
generating and orienting a simulation box, and filing it with pre-generated 
topologies. Creating a simulation box with ``editconf`` can
be done with the following:

```
gmx editconf -f ${SYSTEM}.gro -d ${SEPARATION} -bt {BOX_TYPE} -o ${OUTPUT}.gro
```

where ``${SYSTEM}.gro`` is the input coordinate file, 
``${OUTPUT}.gro`` is the chosen output name (the default is ``out.gro``), 
``-d ${SEPARATION}`` defines the minimum 
separation between the input and the edge of the box, and 
``-bt ${BOX_TYPE}`` defines the type of box for the simulation. There 
are a number of other ``editconf`` options that can be found in the GROMACS 
manual 
[gmx editconf page](http://manual.gromacs.org/documentation/current/onlinehelp/gmx-editconf.html)


> ## Exercise
> 
> Use the editconf command to generate a cubic simulation box for our system.
> The molecule should be centred in the box, and a seperation of 1 nm should be used 
> between the edge of the box and the molecule.
> 
> > ## Solution
> > 
> > `gmx editconf -f conf.gro -c -d 1 -bt cubic -o 5pep-box.gro`
> {: .solution}
{: .challenge}


If you now look at the new ``5pep-box.gro`` file you should see the box
dimensions have changed - the box is now cubic.

### Solvating a system


The ``solvate`` tool can be used to create a box of solvent or 
to solvate a pre-existing box. To solvate our system we will run:


```
gmx solvate -cp 5pep-box.gro -cs spc216.gro -o 5pep-solv.gro -p topol.top
``` 
 
where ``5pep-box.gro`` is our simulation box configured using the steps 
described above, ``spc216.gro`` is the solvent configuration file, and ``topol.top`` is the 
topology from earlier. Here ``spc216.gro`` is compatitable 
with the 3 point tip3p model we are using. 

Now we have ``5pep-solv.gro``. You should also see the message when runnning
solvate:

```
Back Off! I just backed up topol.top to ./#topol.top.1#
```

This means that GROMACS generated a new topol.top file and renamed the old
one to ``'#topol.top.1#'``. This is the default behavoiur when GROMACS 
creates files with the same name as a file in your current directory.
The new ``topol.top`` is not very different to the old one only now the
SOL has been added to the ``[ molecules ]`` directive.

```
[ molecules ]
; Compound        #mols
Protein_chain_A     1
SOL             19181
```


For further information on solvate, please see the GROMACS manual 
[gmx solvate](http://manual.gromacs.org/documentation/current/onlinehelp/gmx-solvate.html)

### Adding ions and creating a charge-neutral system


Adding ions to your solvated system can serve two purposes: it can help to 
neutralise any charge in your system; and it allows you to simulate systems 
with similar salt concentrations to their real-world equivalents.

You may have noticed earlier that our system is charged. Here we will
neutralise it by adding ions.

Adding 
ions is done in two parts: first, you need to use the ``grompp`` tool to 
generate a ``.tpr`` file to be used when adding ions, and then you must 
replace some of the recently-added solvent molecules with the necessary 
counterions using ``genion``.

The GROMACS preprocessor tool ``grompp`` reads in coordinate and topology 
files to generate an atomic-level input file (with a ``.tpr`` extension). 
This ``.tpr`` file contains all of the parameters needed for all atoms in 
the system. 


In order to generate a run input ``.tpr`` file, ``grompp`` needs a structure (``.gro``) 
file, a topology (``.top``) file, and a file defining the instructions for 
the simulation run (this is kept in an ``.mdp`` file). This ``.mdp`` file can 
be empty when ionising the system as no actual simulation is run. 
To generate the ``.tpr`` file, run:


```
# create and empty mdp file 
touch mdrun.mdp
gmx grompp -f mdrun.mdp -c 5pep-solv.gro -p topol.top -o ions.tpr 
``` 
 
 
``grompp`` will output a number of notes to screen (one of 
which should be reminding you of the net non-zero charge of your system). In 
this case, these can be ignored.

In some cases you may see warnings at this point. These should not be ignored
unless you are sure that your set up is correct.

Now that the ``.tpr`` has been generated, ``genion`` can be used to make the 
charge of the system neutral. The system charge is decreased by replacing a 
number of parts of the system with anions and cations. This is done by 
running the following: 

```
gmx genion -s ions.tpr -p topol.top -neutral -o 5pep-neutral.gro
```


You are then asked to choose the group within your system that the ions will replace.
Here we will replace some of the solvent `SOL` molecules with ions. This ensures
that the protein itself remains intact.

```
Reading file ions.tpr, VERSION 2021.3 (single precision)
Will try to add 38 NA ions and 0 CL ions.
Select a continuous group of solvent molecules
Group     0 (         System) has 62225 elements
Group     1 (        Protein) has  4682 elements
Group     2 (      Protein-H) has  2426 elements
Group     3 (        C-alpha) has   326 elements
Group     4 (       Backbone) has   978 elements
Group     5 (      MainChain) has  1305 elements
Group     6 (   MainChain+Cb) has  1596 elements
Group     7 (    MainChain+H) has  1618 elements
Group     8 (      SideChain) has  3064 elements
Group     9 (    SideChain-H) has  1121 elements
Group    10 (    Prot-Masses) has  4682 elements
Group    11 (    non-Protein) has 57543 elements
Group    12 (          Water) has 57543 elements
Group    13 (            SOL) has 57543 elements
Group    14 (      non-Water) has  4682 elements
Select a group: 13
```

Once the group is chosen, ``genion`` will replace a number of that 
group with anions and cations until the system is charge neutral.

In our case 38 NA atoms are replaced in the water molecules.

You should see the following, indicating the water molecules that will be
replaced by soduim ions:

```
Replacing solvent molecule 800 (atom 7082) with NA
Replacing solvent molecule 618 (atom 6536) with NA
Replacing solvent molecule 6992 (atom 25658) with NA
Replacing solvent molecule 8244 (atom 29414) with NA
Replacing solvent molecule 3746 (atom 15920) with NA
Replacing solvent molecule 5204 (atom 20294) with NA
Replacing solvent molecule 17626 (atom 57560) with NA
Replacing solvent molecule 17694 (atom 57764) with NA
Replacing solvent molecule 14834 (atom 49184) with NA
Replacing solvent molecule 7080 (atom 25922) with NA
Replacing solvent molecule 8206 (atom 29300) with NA
Replacing solvent molecule 16326 (atom 53660) with NA
Replacing solvent molecule 2039 (atom 10799) with NA
Replacing solvent molecule 11793 (atom 40061) with NA
Replacing solvent molecule 11087 (atom 37943) with NA
Replacing solvent molecule 1214 (atom 8324) with NA
Replacing solvent molecule 822 (atom 7148) with NA
Replacing solvent molecule 11826 (atom 40160) with NA
Replacing solvent molecule 1164 (atom 8174) with NA
Replacing solvent molecule 8729 (atom 30869) with NA
Replacing solvent molecule 2675 (atom 12707) with NA
Replacing solvent molecule 7510 (atom 27212) with NA
Replacing solvent molecule 2319 (atom 11639) with NA
Replacing solvent molecule 5374 (atom 20804) with NA
Replacing solvent molecule 8780 (atom 31022) with NA
Replacing solvent molecule 12110 (atom 41012) with NA
Replacing solvent molecule 4713 (atom 18821) with NA
Replacing solvent molecule 17524 (atom 57254) with NA
Replacing solvent molecule 965 (atom 7577) with NA
Replacing solvent molecule 2713 (atom 12821) with NA
Replacing solvent molecule 8711 (atom 30815) with NA
Replacing solvent molecule 7607 (atom 27503) with NA
Replacing solvent molecule 3728 (atom 15866) with NA
Replacing solvent molecule 15300 (atom 50582) with NA
Replacing solvent molecule 17051 (atom 55835) with NA
Replacing solvent molecule 16208 (atom 53306) with NA
Replacing solvent molecule 14795 (atom 49067) with NA
Replacing solvent molecule 10177 (atom 35213) with NA
```


For further information, please see the GROMACS manual  
[gmx grompp](http://manual.gromacs.org/current/onlinehelp/gmx-grompp.html),
and [gmx genion](http://manual.gromacs.org/documentation/current/onlinehelp/gmx-genion.html) 
pages.

We are now ready to run a simulation.


{% include links.md %}

