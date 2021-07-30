##### How to generate a ligand force fields (using Amber Tools)

1. Add hydrogen atoms.
   
   > reduce ligand.pdb > ligand_H.pdb

2. Check a molecular structure of a ligand and then identify its total charge.
   
   

3. Estimate atomic point charges as follows:
   
   > antechamber -fi pdb -fo mol2 -i ligand_H.pdb -o ligand_H.mol2 -c bcc -pf y -nc 3
   
   in case that the total charge of the ligand is 3.

4. Generate a frcmod file.
   
   > parmchk2 -i ligand_H.mol2 -o ligand_H.frcmod -f mol2

5. Generate the parameter and topology file using a General Amber Force Field.
   
   tleap -s -f ligand.in
   
   ---ligand.in---
   
   ```textile
   source leaprc.gaff2
   mol = loadmol2 ligand_H.mol2
   loadamberparams ligand_H.frcmod
   saveamberparam mol ligand.prmtop ligand.inpcrd
   quit
   
   ```
   (Note, there is a empty line after 'quit'.)

6. Build a force field file for OpenMM from an AMBER prmtop file.
   
   > python write_xml_pretty.py -i input_ff.json

---input_ff.json----

```
{
    "fname_prmtop" : "ligand.prmtop",
    "fname_xml" : "ligand.ff.xml",
    "ff_prefix" : "lig"
}
```


