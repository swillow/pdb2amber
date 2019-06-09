# pdb2amber
Generate an Amber prmtop file from a pdb file.

# Introduction
OpenMM provides the outstanding concept that the potential energy of the system can be calculated with the pdb file (the forcefields are provided internally).
Sometimes the potential energy calculated with the pdb file differs from the one obtained with the amber prmtop and inpcrd files.
Here, in order to figure out why, I made a code which generates an amber prmtop file from a pdb file.

# Note
_topology.py and _pdbfile.py are topology.py and pdbfile.py at the directory of (...)/site-packages/simtk/openmm/app.
The key modification in _topology.py and _pdbfile.py is that those files read the local data files (./data/residues_new.xml ./data/pdbNames.xml).
I added the bonds of 'CYX' into './data/residues_new.xml'. 

_forcefield.py was made based on forcefield.py at (...)/site-packages/simtk/openmm/app.

# Howto 
python pdb2amber -i receptor.pdb -o receptor.prmtop
