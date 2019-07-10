import os
from decimal import Decimal

# box = [-60:60, -60:60, -75:75]
# cell_list:120/5
cell_list = {}

for ic in range (-12, 12):
    for jc in range (-12, 12):
        for kc in range (-15, 15):
            cell_list[ (ic, jc, kc) ] = []
            
neigh_list = []
for ic in range (-1,2):
    for jc in range (-1,2):
        for kc in range (-1,2):
            neigh_list.append ([ic, jc, kc])
            
f_pdb = open ('PaSDH_opt_H_DPP.pdb')
l_pdb = f_pdb.read().split('\n')
f_pdb.close()

# xrange : -30 ~ 70 : 100
# yrange : -22 ~ 66 : 88
# zrange : -80 ~ 56 : 136
import math
heavy_atoms = []
for line in l_pdb[:-1]:
    key   = line[0:4]
    l_atom = False
    if key == 'ATOM':
        l_atom = True
    if key == 'HETA':
        l_atom = True

    if l_atom:
        atName = line[12:16].strip()
        if atName[0] != 'H':
            words = line[30:54].split()
            x     = Decimal(words[0]) 
            y     = Decimal(words[1]) 
            z     = Decimal(words[2]) 

            if x < Decimal(-60.0):
                x += Decimal(120.0)
            if x >= Decimal(60.0):
                x -= Decimal(120.0)
                
            if y < Decimal(-60.0):
                y += Decimal(120.0)
            if y >= Decimal(60.0):
                y -= Decimal(120.0)
                
#            if z < Decimal(-75.0):
#                z += Decimal(150.0)
#            if z >= Decimal(75.0):
#                z -= Decimal(150.0)
                
            ic = int(math.floor(x/5))
            jc = int(math.floor(y/5))
            kc = int(math.floor(z/5))
            
#            print '(x,y,z) ', x, y, z, ' (ic, jc, kc) ', ic, jc, kc
            
            cell_list[(ic,jc,kc)].append( (x,y,z) )
            
#            heavy_atoms.append ( (x,y,z) )

f_wat = open ('water.xyz')
l_wat = f_wat.read().split('\n')
f_wat.close()

l_wat.pop(0) # natom
l_wat.pop(0) # comment

nwat = 902
wo_list = []
wh1_list = []
wh2_list = []
for iw in range (nwat):
    io  = 3*iw
    ih1 = io+1
    ih2 = io+2

    words = l_wat[io].split()
    x = Decimal(words[1])
    y = Decimal(words[2])
    z = Decimal(words[3])
    wo_list.append ( (x, y, z) )
    
    words = l_wat[ih1].split()
    x = Decimal(words[1])
    y = Decimal(words[2])
    z = Decimal(words[3])
    wh1_list.append ( (x, y, z) )
    
    words = l_wat[ih2].split()
    x = Decimal(words[1])
    y = Decimal(words[2])
    z = Decimal(words[3])
    wh2_list.append ( (x, y, z) )

f_sol= open ('water_cut.pdb', 'w', 1)
# 
# xrange : -30 ~ 70 : 100
# yrange : -22 ~ 66 : 88
# zrange : -80 ~ 56 : 136
nres  = 0
natom = 0
for ix in range (-2, 2): # -60 ~ 60
    cx = Decimal(30.0*ix)
    
    for iy in range (-2, 2): # -60 ~ 60
        cy = Decimal(30.0*iy)
        
        for iz in range (-3, 3): # -90 ~ 90 : cut (-75:75)
            cz = Decimal(30.0*iz)

            print ix, iy, iz, cx, cy, cz
        
            for iw in range(nwat):

                xi, yi, zi  = wo_list[iw]

                xo = xi + cx
                yo = yi + cy
                zo = zi + cz

                xi, yi, zi  = wh1_list[iw]
                xh1 = xi + cx
                yh1 = yi + cy
                zh1 = zi + cz
                
                xi, yi, zi  = wh2_list[iw]
                xh2 = xi + cx
                yh2 = yi + cy
                zh2 = zi + cz
                
                ic = int (math.floor (xo/5))
                jc = int (math.floor (yo/5))
                kc = int (math.floor (zo/5))

                
                if zo <= -75.0:
                    continue
                if zo >= 75.0:
                    continue
                
                if -50.0 < zo and zo < -7.0:
                   continue
               
                l_overlap = 0
                for ic1, jc1, kc1 in neigh_list:
                   ic2 = ic + ic1
                   jc2 = jc + jc1
                   kc2 = kc + kc1
                   
                   if ic2 < -12: ic2 =  11
                   if ic2 > 11:  ic2 = -12
                       
                   if jc2 < -12: jc2 =  11
                   if jc2 > 11:  jc2 = -12
                       
                   if kc2 < -12: kc2 =  11
                   if kc2 > 11:  kc2 = -12

                                
                   for xj, yj, zj in cell_list[ (ic2, jc2, kc2) ]:
                    
                       dx = xo - xj
                       dy = yo - yj
                       dz = zo - zj
                    
                       r2 = dx*dx + dy*dy + dz*dz
                    
                       if (r2 < 16.0):
                           l_overlap = 1
                           break

                if (l_overlap == 0):
                    natom += 1
                    nres  += 1
                    if natom == 100000: natom = 0
                    line = "HETATM" + "%5d"%natom + "  O   HOH " + \
                      "%5d    "%nres + \
                      "%8.3f%8.3f%8.3f"%(xo, yo, zo) + \
                      "  1.00  0.00           O"
                    
                    f_sol.write(line + '\n')
                    
                    natom += 1
                    if natom == 100000: natom = 0
                    line = "HETATM" + "%5d"%natom + "  H1  HOH " + \
                      "%5d    "%nres + \
                      "%8.3f%8.3f%8.3f"%(xh1, yh1, zh1) + \
                      "  1.00  0.00           H"
                    f_sol.write(line + '\n')
                    
                    natom += 1
                    if natom == 100000: natom = 0
                    line = "HETATM" + "%5d"%natom + "  H2  HOH " + \
                      "%5d    "%nres + \
                      "%8.3f%8.3f%8.3f"%(xh2, yh2, zh2) + \
                      "  1.00  0.00           H"
                    f_sol.write(line + '\n')
                    
