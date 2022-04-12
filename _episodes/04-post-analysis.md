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
- "Learn how to obtain the radial distribution function"

keypoints:
- "Use of gmx energy, make_ndx, rdf, msd, velacc"
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
  gmx energy -f ener.edr -o total_en.xvg
```

When running this, you will get a prompt asking which property you would like 
output (*e.g.* potential energy, kinetic energy, pressure, temperature, 
*etc.*). 

```
End your selection with an empty line or a zero.
-------------------------------------------------------------------
  1  Angle            2  Proper-Dih.      3  Improper-Dih.    4  LJ-14         
  5  Coulomb-14       6  LJ-(SR)          7  Disper.-corr.    8  Coulomb-(SR)  
  9  Coul.-recip.    10  Potential       11  Kinetic-En.     12  Total-Energy  
 13  Conserved-En.   14  Temperature     15  Pres.-DC        16  Pressure      
 17  Constr.-rmsd    18  Box-X           19  Box-Y           20  Box-Z         
 21  Volume          22  Density         23  pV              24  Enthalpy      
 25  Vir-XX          26  Vir-XY          27  Vir-XZ          28  Vir-YX        
 29  Vir-YY          30  Vir-YZ          31  Vir-ZX          32  Vir-ZY        
 33  Vir-ZZ          34  Pres-XX         35  Pres-XY         36  Pres-XZ       
 37  Pres-YX         38  Pres-YY         39  Pres-YZ         40  Pres-ZX       
 41  Pres-ZY         42  Pres-ZZ         43  #Surf*SurfTen   44  Box-Vel-XX    
 45  Box-Vel-YY      46  Box-Vel-ZZ      47  T-Protein       48  T-SOL         
 49  T-NA            50  Lamb-Protein    51  Lamb-SOL        52  Lamb-NA       
```


Enter the `12` to generate an XMGrace file that includes the variation
of the Total energy over time. This can be plotted with xmgrace or other
plotting software. 
There are a number of other options for the ``energy`` command, and these 
can be found in the GROMACS manual 
`gmx energy<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-energy.html#gmx-.. energy>`_
page.

Generating an index file
------------------------

GROMACS has a post-analysis tool for generating radial distribution functions 
(RDFs). Before generating an RDF, we will need to create a GROMACS index 
(``.ndx``) file to categorise the various parts that compose the simulation 
into indexed groups. This can be done with the ``gmx make_ndx`` command. To 
use it, run:

```
  gmx make_ndx -f confout.gro -o 5pep.ndx
```

 Provided you used the default names in 
your ``mdrun``, you can simply use ``confout.gro``. The ``make_ndx`` command 
will analyse the system, and output the default index groups. 

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

It is possible to create new index groups by using the command prompts listed.
For now we will just save and quit with option ``q``.

If you now open the ``5pep.ndx`` file we just generated you should see
the different groups of atoms listed. Each number is an atom number from the 
``confout.gro file``


For more information on ``make_ndx``, please see the
GROMACS manual
`gmx make_ndx<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-make_ndx.html>`_ 
page.


For more complex manipulations than selecting all of one group of atoms, 
GROMACS provides the ``gmx select`` option. This will allow you to define 
the exact time or particles or regions of interest within your simulation. 
You can find more information on how to use this in the GROMACS manual
`Groups and Selections<https://manual.gromacs.org/documentation/2019/reference-manual/analysis/using-groups.html#selections>`_
page.


Radial distribution function
----------------------------

Once an appropriate index file is generated, with the atoms for which an RDF 
is to be calculated indexed into appropriate groups, we can use the 
``gmx rdf`` command to generate the RDFs. This is done by running:

```
  gmx rdf -f ${TRAJECTORY_INPUT}.trr -n ${INDEX_INPUT}.ndx  \
          -ref ${REFERENCE_GROUP} -sel ${SELECTED_GROUP} -bin ${BIN_WIDTH}
          -o ${OUTPUT}.xvg
```
  
where ``${TRAJECTORY_INPUT}.trr`` is the trajectory file for which you would 
like to generate an RDF, and ``${INDEX_INPUT}.ndx`` is the index file that you 
produced using ``make_ndx``. ``${REFERENCE_GROUP}`` should be replaced with 
the name of the principal group to be used in the RDF as it appears in the 
``${INDEX_INPUT}.ndx`` file. Likewise, ``${SELECTED_GROUP}`` should be 
replaced with the name of the atom group(s) for which you want to calculate 
the RDF against the position of the reference group (*e.g.* if you want to 
calculate the RDF between sodium ions and chloride ions, your reference 
group would be one of ``NA`` or ``CL``, and your selected group would be the 
one not chosen as reference). Note that it is possible for your reference and 
selected groups to be the same group.

Mean squared displacement and velocity autocorrelation functions
----------------------------------------------------------------

Gromacs offers a number of tools to calculate correlation and autocorrelation 
functions. Here, we will look at two specific example: the mean-squared 
displacement (MSD) and velocity autocorrelation function (VACF). We will focus 
on how to generate these functions within GROMACS but you can use these links 
to find an overview of the theory behind the 
`MSD<http://manual.gromacs.org/documentation/current/reference-manual/analysis/mean-square-displacement.html>`_
and the 
`VACF<http://manual.gromacs.org/documentation/2019/reference-manual/analysis/correlation-function.html>`_.

Calculating the MSD of parts of a system can be done using the ``gmx msd``. 
This can be run using:

```
  gmx msd -f ${INPUT_TRAJECTORY}.trr -s ${INPUT_TOPOLOGY}.tpr -o ${OUTPUT}.xvg
```


where ``${INPUT_TRAJECTORY}.trr`` is the trajectory file of the simulation for 
which the MSD is being calculated, and ``${INPUT_TOPOLOGY}.tpr`` can be the 
input file used to obtain this trajectory (note that it is possible to use 
the final topology ``confout.gro`` file here instead to obtain the same 
results). Running this command will prompt you to choose the group for which 
you would like the MSD. Note that, if the group you are looking for is not 
present in the list, you can generate an index file (see 
`Generating an index file`_) where you can define this new group. To include 
this index file, add the option ``-n ${INDEX_FILE}.ndx`` to the command above.
For more information and options, please look at the GROMACS manual page on 
the `gmx msd command<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-msd.html#gmx-msd>`_.

VACFs can be generated using the ``gmx velacc`` command:

```
  gmx velacc -f ${INPUT_TRAJECTORY}.trr -o ${OUTPUT}.xvg
```


where ``${INPUT_TRAJECTORY}.trr`` is the trajectory file of the simulation 
for which the VACF is being produced. You will get a prompt asking for which 
group of atoms the VACF should be calculated. If the group you want is not 
present, you may need to create it by following the instructions in the 
`Generating an index file`_ section of the manual. To include your index file, 
add it with the ``-n ${INPUT_INDEX}.ndx`` option. You can find more options 
and information on the GROMACS manual 
`gmx velacc<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-velacc.html#gmx-velacc>`_ page.



{% include links.md %}


