---
title: "Preparing a GROMACS system"
teaching: 30
exercises: 30
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

As an example we will be looking at pepsin. The PDB file can be
downloaded from the protein data bank
using

```
wget https://files.rcsb.org/download/5pep.pdb
```


This protein is used as an example but this process could apply to a general
system where
you start from a .pdb file and end up with a set of input files, ready to
run a GROMACS simulation.

### Splitting up the system

A pdb file downloaded from the database may not be in the format you want
for running your simulation. For example pdb files can contain non-protein
residues and waters. In this case we have some waters already in the pdb
file that we want to remove as we will be solvating the system later on.
You can see these in the pdb file if you open it. 

To just get the protein itself we can use the following grep command to remove
lines containing the water symbol 'HOH'.

```
grep -v 'HOH' 5pep.pdb > 5pep_protein.pdb
```

The new 5pep_protein.pdb file contains just the protein itself.

### Creating a Gromacs topology (PDB2GMX)


Now we can create a Gromacs topolopy for the system.
GROMACS ``pdb2gmx`` command is used to convert a coordinate file into a 
set of GROMACS topology files and also create a processed structure file 
in the GROMACS format .gro. 

First you will need to load the GROMACS module on ARCHER2:

```
module load gromacs
```
Then you can run:

```
gmx pdb2gmx -f 5pep_protein.pdb
```

You will be 
prompted to select the forcefield you would like to use. GROMACS comes with 
a number of AMBER and GROMOS forcefields, as well as a CHARMM and an OPLS-AA
option. You will also need to specify your water model (choices included are 
TIP models, and the SPC and SPC/E models). Specifying the water model here 
results in ``pdb2gmx`` to write a complete topology and will ensure that all
topology files are consistent if the system needs to be hydrated.

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
is particularly important to note down and will be used in the `Solvating and 
ionise a system`_ step of system preparation.

More information about the flags and options of this program can be found in 
the GROMACS 
`PDB2GMX manual<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html>`_.

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

The [ moleculetype ] gives the molecule name and the number of exclusions. 
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

The bonds directive lists pairs of atoms which are bonded
```
[ bonds ]       
;  ai    aj funct            c0            c1            c2            c3
    1     2     1      
    1     3     1      
```    
 The pairs directive lists LJ pairs of atoms by their atom number
 ```
 [ pairs ]
;  ai    aj funct            c0            c1            c2            c3
    1     8     1
    1     9     1
 ```   
 The angles directive lists triplets of 3 atoms for angle parameterisation
 ```
 [ angles ]
;  ai    aj    ak funct            c0            c1            c2            c3
    2     1     3     1
    2     1     4     1
```
The dihedreal 
```
[ dihedrals ]
;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5
    2     1     5     6     9
    2     1     5     7     9
```

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
pdb2gmx) contains posisition restraits for the atoms. The if 
statement ensures it is only applied when POSRES is true.

* The water forcefield tip3p is included

* The position restraints for waters

* The forcefield for ions

* The name of the system

* The molecules in the system (so far just one protein chain

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

The atomic coordinates in x,y,z of each atom are then listed. These
are preceeded by the residue number, residue name, atom name
and atom number.

At the very end of the file is the box dimensions listed as:

```
box[X][X],box[Y][Y],box[Z][Z], box[X][Y],box[X][Z],box[Y][X],box[Y][Z],box[Z][X],box[Z][Y]
```

#### The posre.itp file

This file contains the position restraints for the system. The force
constants in x,yz are listed for the atom numbers.

### Generating your own forcefield file

Occasionally a system will contain non-protein residues which need separate 
parameterisation as there are not forcefields available for them. You will 
see a message when running pdb2gmx such as:

```

```

In this case you must generate your own forcefield. In GROMACS to do this
you can generate your own forcefield files in a {forcefield}.ff directory.
More about this can be found on the Gromacs manual.

`file format page<http://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#top>`_.

You can also generate an Amber or CHARMM topology by using the   AmberTool 
 ``antechamber`` function or the CHARMM ``cgenff`` function. To do this, you 
should follow the procedures described above, making sure to select an 
appropriate forcefield from the selection GROMACS provides. Then, use a 
parameter-generating tool like ``antechamber`` with ``actype`` (for Amber) 
or ``cgenff`` (for CHARMM). The topologies generated in this way can then be 
added to the GROMACS topology that you generated. 
  

.. note::

  It is advised to use one of the pre-existing GROMACS forcefields where 
  possible. Only consider generating your own forcefield if the ones 
  available on GROMACS do not fulfill your requirements.

GROMACS comes with a number of forcefields available out of the box. These 
forcefields are stored within the main GROMACS directory in 
``share/gromacs/top``. If the forcefield you want to use is not present, you
will need to generate your own forcefield files. To do this, create a 
directory ``<forcefield>.ff`` in your working directory, with ``<forcefield>``
replaced by a sensible name. Within this directory, create a 
``forcefield.doc`` file and write a simple one-sentence description of your 
forcefield -- this description is what will come up in ``pdb2gmx`` when you 
choose a forcefield. Next, generate a ``forcefield.itp`` included topology 
file. This file is a topology file where you can define the parameters for 
atoms, bond, angles, dihedrals, *etc.*. You can find more information about 
generating topology files from scratch in the GROMACS manual 
`file format page<http://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#top>`_.

Using a ``<forcefield>.ff`` directory has a number of advantages over writing 
out your system topologies directly. For one, this allows for better 
reproducibility in the even that you want to simulate a new system with this 
forcefield. It also has a number of functionalities that can be useful. For 
instance, adding a ``watermodels.dat`` file into the forcefield directory 
makes it easy to keep track of water models available. A line and descripting 
can be added in this file for each water-model included topology file. This 
file is what prompts the choice of water model in ``pdb2gmx``.

Once it is populated, running ``pdb2gmx`` in the directory containing your 
``<forcefield>.ff`` directory will result in your new forcefield being included 
at the top of the list of selectable forcefields. If you are happy with your 
``<forcefield>.ff`` directory and you will use it a lot (and if you have the 
correct permissions to edit parts of the GROMACS directory), you can copy it to 
the ``share/gromacs/top`` of the GROMACS directory (or to ``$GMXLIB`` which 
should be the same directory). In doing so, your forcefield will become a 
permanent part of the forcefields that ``pdb2gmx`` can use.

.. note::

  You can also generate an Amber or CHARMM topology by using the   AmberTool 
  ``antechamber`` function or the CHARMM ``cgenff`` function. To do this, you 
  should follow the procedures described above, making sure to select an 
  appropriate forcefield from the selection GROMACS provides. Then, use a 
  parameter-generating tool like ``antechamber`` with ``actype`` (for Amber) 
  or ``cgenff`` (for CHARMM). The topologies generated in this way can then be 
  added to the GROMACS topology that you generated. This can be done by 
  opening the GROMACS topology file and including the following line at the start:
  

```
    #include "/path/to/forcefield_file.itp"
```


  where the path is to the topology file generated in ``antechamber`` or 
  ``cgenff``.

For more information on generating your own forcefield, please see the GROMACS
manual pages about 
`adding a residue<http://manual.gromacs.org/documentation/current/how-to/topology.html>`_
and `force field organisations<http://manual.gromacs.org/documentation/current/reference-manual/topologies/force-field-organization.html>`_.


## Preparing and solvating your simulation box


### Generating a simulation box


Now that a topology has been generated, the next step is to generate a 
simulation box into which to place this topology. For this, use the 
``editconf`` command. This tool has a number of functionalities, including 
generating and orienting a simulation box, and filing it with pre-generated 
topologies. To create the simulation box with ``editconf``, run the following:

```
gmx editconf -f conf.gro -c -d 1 -bt cubic -o 5pep-box.gro
```

where ``coonf.gro`` is the input forcefield-compliant coordinate file, 
``5pep-box.gro`` is the chosen output name (the default is ``out.gro``), 
the ``-c`` flag will place the system described in ``conf.gro`` into the 
centre of the simulation box, ``-d ${SEPARATION}`` defines the minimum 
separation between the input and the edge of the box (here 1 nm), and 
``-bt ${BOX_TYPE}`` defines the type of box for the simulation (triclinic is 
the default, but other options are cubic, octohedral, or dodecahedral). There 
are a number of other ``editconf`` options, predominantly to have more 
control over defining the simulation box. These can be found in the GROMACS 
manual 
`gmx editconf page<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-editconf.html>`_.

### Solvating a system


The aptly-named ``solvate`` tool can be used to create a box of solvent or 
to solvate a pre-existing box. To use it, run:


```
gmx solvate -cp 5pep-box.gro -cs spc216.gro -o 5pep-solv.gro -p topol.top
``` 
 
where ``5pep-box.gro`` is the simulation box configured using the steps 
described above, ``spc216.gro`` is the solvent configuration file (note 
that GROMACS has a number of pre-defined solvent configuration files but that 
you can also prepare and use your own), and ``topol.top`` is the 
topology obtained when running `GMX2PDB`. Here ``spc216.gro`` is compatitable 
with the 3 point tip3p model we are using. _ If using a GROMACS-provided 
solvent, the addition of this solvent should not alter the net charge of the 
system.

For further information, please see the GROMACS manual 
`gmx solvate<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-solvate.html>`_

### Adding ions and creating a charge-neutral system


Adding ions to your solvated system can serve two purposes: it can help to 
neutralise any charge in your system; and it allows you to simulate systems 
with similar salt concentrations to their real-world equivalents. Adding 
ions is done in two parts: first, you need to use the ``grompp`` tool to 
generate a ``.tpr`` file to be used when adding ions, and then you must 
replace some of the recently-added solvent molecules with the necessary 
counterions using ``genion``.

The GROMACS preprocessor tool ``grompp`` reads in coordinate and topology 
files to generate an atomic-level input file (with a ``.tpr`` extension). 
This ``.tpr`` file contains all of the parameters needed for all atoms in 
the system. We will go into more details about the ``grompp`` tool in the 
`Running a simulation`_ section. For now, the important part is that, to 
generate a run input ``.tpr`` file, ``grompp`` needs a structure (``.gro``) 
file, a topology (``.top``) file, and a file defining the instructions for 
the simulation run (this is kept in an ``.mdp`` file). This ``.mdp`` file can 
be kept empty when ionising the system as no actual simulation is to be run. 
To generate the ``,tpr`` file, run:


```
# create and empty mdp file 
touch mdrun.mdp
gmx grompp -f mdrun.mdp -c 5pep-solv.gro -p topol.top -o ions.tpr 
``` 
 
 
It is likely that ``grompp`` will output a number of notes to screen (one of 
which should be reminding you of the net non-zero charge of your system). In 
this case, these can be ignored (this is an exception and is not usually true).

Now that the ``.tpr`` has been generated, ``genion`` can be used to make the 
charge of the system neutral. The system charge is decreased by replacing a 
number of parts of the system with anions and cations. This is done by 
running the following (note that the ``${INPUT}.tpr`` named below is likely 
to be the ``${OUTPUT.tpr}`` generated in the ``grompp`` step above): 

```
gmx genion -s ions.tpr -p topol.top -neutral -o 5pep-neutral.gro
```

You will be prompted to choose the group within your system (solvents, 
solutes, protein backbones, *etc.*) that you would like ions to replace, with 
the frequency of occurrence of each group also shown. Note that some groups 
may have overlap completely and be different names for the same group. In 
general, it is best to replace solvent molecules with ions (the group named 
``SOL``). Once a group is chosen, ``genion`` will replace a number of that 
group with anions and cations until the system is charge neutral. The default 
anion name is ``CL``, though this name can be changed with the ``-nname`` 
flag, and the default cation name is ``NA``, but this name can be changed with 
the ``nname`` flag. By default, the cation and anion charges are 1 and -1 
respectively, but this can be changed with the ``-pq`` flag for the cation and 
the ``-nq`` flag for the anion.

For further information, please see the GROMACS manual  
`gmx grompp<http://manual.gromacs.org/current/onlinehelp/gmx-grompp.html>`_, 
and `gmx genion<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-genion.html>`_ 
pages.



{: .challenge}

{% include links.md %}

