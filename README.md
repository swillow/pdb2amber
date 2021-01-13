# pdb2amber

Generate an AMBER prmtop file from a pdb file.

# Introduction

OpenMM offers an excellent concept that the potential energy of the system can be calculated with the pdb file since the **standard** forcefields are provided internally. 
In special cases where non-standard amino acids and ligands are in the pdb file, we need to provide additional forcefields.
Here, 'build_forcefield' introduces how to build an OpenMM forcefield file (in a xml file) from an AMBER prmtop file. With these additional OpenMM forcefield files, we generate an AMBER prmtop file from a PDB file. 

# Note

_topology.py, _pdbfile.py, and _forcefield.py were modified codes of topology.py, pdbfile.py, and forcefield.py at the directory of (...)/site-packages/simtk/openmm/app for this goal.


# Howto

python pdb2amber.py -i input.json

---`input.json'---

```
{
    "fname_pdb" : "abc.pdb",
    "fname_prmtop" : "abc.prmtop",
    "fname_ff" : [
        "./data/protein.ff14SB.xml",
        "./data/wat_opc3.xml"
    ]
}
```
