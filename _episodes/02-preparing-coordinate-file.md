---
title: "Preparing a GROMACS system"
teaching: 30
exercises: 30
questions:
- "How do I set up a system ready to be run with "
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
run a Gromacs simulation.


### Creating a Gromacs topology (PDB2GMX)


The first step is to create a Gromacs topolopy for the system.
GROMACS ``pdb2gmx`` command is used to convert a coordinate file into a 
set of GROMACS topology files (in these examples, we will assume that the 
file is a ``.pdb`` file, but this is not a necessity). To run this:

```
gmx pdb2gmx -f 5pep.pdb
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
  gmx pdb2gmx -f 5PEP.pdb

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
5PEP.pdb  conf.gro  posre.itp  topol.top
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

### Generating your own forcefield file


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


### Generating a system of replicates from a GROMACS structure file


It is possible to populate a simulation box with molecules by replicating the contents 
of a GROMACS structure file (``.gro``) multiple times. This can be achieved 
with the ``insert-molecules`` command. While any structure file can be used 
(including crowded system file), this is particularly useful if you want to 
create a system with a large number of copies of a single molecule (*i.e.* 
as found in a lipid bilayer or a non-aqueous solvent). Furthermore, the 
topology (``.top``) file generated for the system to be replicated will still 
work for the new, larger system, by including the total number of molecules in the directive [molecules].

To generate a system using this command, run:

```
gmx insert-molecules -ci ${INPUT}.gro -o ${OUTPUT}.gro \
                       -nmol ${N} -box ${X_LENGTH} ${Y_LENGTH} ${Z_LENGTH}
 ```
 
 
where ``${INPUT}.gro`` is the structure file of the molecule/system you wish 
to replicate, ``${OUTPUT}.gro`` is the output file, ``${N}`` is the number of 
times that the contents of ``${INPUT}.gro`` will be replicated, and 
``${X_LENGTH}``, ``${Y_LENGTH}``, and ``${Z_LENGTH}`` are the dimensions of 
the cubic box into which these ``${N}`` replicas must be packed.

There are number of further options to help pack your system, including a way 
of defining the default van der Waals distance between atoms in your system, a 
way of inserting new molecules into an existing system, and methods to control 
the amount of random rotation that replicated molecules can undergo. All of 
these options can be found in the 
`gmx insert-molecules<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-insert-molecules.html>`_ page of the GROMACS manual.

### Generating a simulation box


Now that a topology has been generated, the next step is to generate a 
simulation box into which to place this topology. For this, use the 
``editconf`` command. This tool has a number of functionalities, including 
generating and orienting a simulation box, and filing it with pre-generated 
topologies. To create a simulation box with ``editconf``, run:

```
  gmx editconf -f ${INPUT}.gro -center -d ${SEPARATION} -bt ${BOX_TYPE} \
               -o ${OUTPUT}.gro
```

where ``${INPUT}.gro`` is the input forcefield-compliant coordinate file, 
``${OUTPUT}.gro`` is the chosen output name (the default is ``out.gro``), 
the ``-c`` flag will place the system described in ``${INPUT}.gro`` into the 
centre of the simulation box, ``-d ${SEPARATION}`` defines the minimum 
separation between the input and the edge of the box (units are in nm), and 
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
  gmx solvate -cp ${SOLUTE}.gro -cs ${SOLVENT}.gro -p ${TOPOLOGY}.top \
              -o ${OUTPUT}.gro
``` 
 
where ``${SOLUTE}.gro`` is the simulation box configured using the steps 
described above, ``${SOLVENT}.gro`` is the solvent configuration file (note 
that GROMACS has a number of pre-defined solvent configuration files but that 
you can also prepare and use your own), and ``${TOPOLOGY}.top`` is the 
topology obtained when running `GMX2PDB`_. If using a GROMACS-provided 
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
  gmx grompp -f ${RUN_FILE}.mdp -c ${COORDINATES}.gro -p ${TOPOLOGY}.top \
             -o ${OUTPUT}.tpr
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
gmx genion -s {INPUT}.tpr -p ${TOPOLOGY}.top -neutral -o ${OUTPUT}.gro
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

