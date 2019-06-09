import parmed as pmd

amber = pmd.load_file ('receptor.prmtop', 'receptor.inpcrd')
amber.save ('receptor.pdb')

