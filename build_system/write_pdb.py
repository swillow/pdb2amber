import sys
from decimal import Decimal

f_pdb = open ('dppe1.pdb', 'r')
l_pdb = f_pdb.read().split('\n')
f_pdb.close()

pdb_line = []
for line in l_pdb:
    str_key = line[0:4]
    l_atom = False
    if str_key == 'ATOM':
        l_atom = True
    if str_key == 'HETA':
        l_atom = True
       
    if l_atom:
        pdb_line.append(line)

#for line in pdb_line:
#    print line
    
#sys.exit(1)

f_rst = open ('bomdrr.sav.box75', 'r')
l_rst = f_rst.read().split('\n')
f_rst.close()

natom = int(l_rst[0].split()[0])
l_rst.pop(0) # nmol
l_rst.pop(0) # nstep

for ires in range (8):
    
    for ii in range (natom):
    iatom = ii*2
    words = l_rst[iatom].split()
    x     = Decimal(words[0])
    y     = Decimal(words[1])
    z     = Decimal(words[2])
    

    str_xyz = "%8.3f%8.3f%8.3f"%(x, y, z)
    line_old  = pdb_line[ii]
#    print 'OLD ', line_old
    line_new  = line_old[0:30] + str_xyz + line_old[54:78]

    print line_new
    
