natom = 5
nres  = 5
x     = 1.0
y     = 2.0
z     = 3.0
line = "HETATM " + "%4d"%natom + "  O   HOH  " + "%4d    "%nres +\
       "%8.3f%8.3f%8.3f"%(x, y, z) + "  1.00  0.00           O"
print line

