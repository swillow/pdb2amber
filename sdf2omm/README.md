##### Using openforcefield, estimate an energy of a small molecule (ligand).

1. Build a sdf file from a pdb file.
   
   > obabel -ipdb ligand.pdb -osdf > ligand.sdf

2. Run sdf2omm.py

   > python sdf2omm.py

   [Energy: -461.7389 kJ/mol]
   

3. To Estimate an Amber potential energy
   
   > python test_amber.py
   
   [Energy: -722.55 kJ/mol]
