import os
import sys
import math
from decimal import Decimal

f_prt = open ('PaSDH_opt_H.pqr')
l_prt = f_prt.read().split('\n')
f_prt.close()

max_prt = [-100, -100, -100]
min_prt = [ 100,  100,  100]

# xrange : -30 ~ 70 : 100
# yrange : -22 ~ 66 : 88
# zrange : -80 ~ 56 : 136
f_pdb = open ('membrane.pdb', 'w', 1)
a_len = 120.0
b_len = 120.0
c_len = 150.0
alpha = 90.0
beta  = 90.0
gamma = 90.0

line_head = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 " % (
    a_len, b_len, c_len, alpha, beta, gamma)
f_pdb.write(line_head + '\n')

heavy_atoms = []
for line in l_prt[:-1]:
    key   = line[0:4]
    l_atom = False
    if key == 'ATOM':
        l_atom = True
    if key == 'HETA':
        l_atom = True

    if l_atom:
        words = line[30:54].split()
        x     = Decimal(words[0]) - Decimal(43.0+35.0)
        y     = Decimal(words[1]) - Decimal(50.6+35.0)
        z     = Decimal(words[2]) - Decimal(65.5+55.0)

        str_xyz = "%8.3f%8.3f%8.3f"%(x, y, z)

        line_new = line[0:30]+str_xyz+line[54:]

        #print line_new
        
    
        f_pdb.write(line_new + '\n')

        atName = line[12:16].strip()
        if atName[0] != 'H':
            heavy_atoms.append ( (x,y,z) )
            if max_prt[2] < z: max_prt[2] = z
            if min_prt[2] > z: min_prt[2] = z

            if z < -5:
                if max_prt[0] < x: max_prt[0] = x
                if max_prt[1] < y:
                    max_prt[1] = y
    
                if min_prt[0] > x: min_prt[0] = x
                if min_prt[1] > y: min_prt[1] = y
        
#f_xyz.close()

print 'min ', min_prt[0],  min_prt[1], min_prt[2]
print 'max ', max_prt[0],  max_prt[1], max_prt[2]
print 'heavy atoms ', len (heavy_atoms)

#sys.exit(1)

f_mem = open ('dppe1.pdb')
l_mem = f_mem.read().split('\n')
f_mem.close()

f_rst = open ('bomdrr.sav.box75')
l_rst = f_rst.read().split('\n')
f_rst.close()
l_rst.pop(0) # natom
l_rst.pop(0) # time step

max_dpp = [-100, -100, -100]
min_dpp = [ 100,  100,  100]
tatom = 0
res_atoms_list = []
for ires in range (8):
    res_atoms = []
    for iatom in range (121):
        words  = l_rst[2*tatom].split()
        x = Decimal(words[0]) 
        y = Decimal(words[1]) 
        z = Decimal(words[2]) - Decimal(30.)
        res_atoms.append ( (x,y,z) )
        tatom += 1
        
        if max_dpp[0] < x: max_dpp[0] = x
        if max_dpp[1] < y: max_dpp[1] = y
        if max_dpp[2] < z: max_dpp[2] = z
    
        if min_dpp[0] > x: min_dpp[0] = x
        if min_dpp[1] > y: min_dpp[1] = y
        if min_dpp[2] > z: min_dpp[2] = z
            
    res_atoms_list.append (res_atoms)

print 'min ', min_dpp[0],  min_dpp[1], min_dpp[2]
print 'max ', max_dpp[0],  max_dpp[1], max_dpp[2]


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

for ia in range (len(heavy_atoms)):
    x, y, z = heavy_atoms[ia]
    
    if x < Decimal(-60.0):
        x += Decimal(120.0)
    if x >= Decimal(60.0):
        x -= Decimal(120.0)
                
    if y < Decimal(-60.0):
        y += Decimal(120.0)
    if y >= Decimal(60.0):
        y -= Decimal(120.0)
        
    ic = int(math.floor(x/5))
    jc = int(math.floor(y/5))
    kc = int(math.floor(z/5))

    cell_list[(ic,jc,kc)].append( (x,y,z) )
    
#f_pdb = open ('build.xyz', 'w', 1)
# 
# xrange : -30 ~ 70 : 100
# yrange : -22 ~ 66 : 88
# zrange : -80 ~ 56 : 136
lbox = [ [0.0,0.0], [0.0,0.0], \
         [0.0,7.5], [0.0,7.5], \
         [7.5,0.0], [7.5,0.0], \
         [7.5,7.5], [7.5,7.5] ] 
             
n_mem = 0
for ix in range (-4, 4):
    cx = Decimal(15.0*ix)
    for iy in range (-4, 4):
        cy = Decimal(15*iy)

        print ix, iy, cx, cy

        

        for ires in range (8):
            print 'ires ---', ires
            res_atoms = res_atoms_list[ires]
            
            l_overlap = 0
            
            for iatom in range (121):
                xi, yi, zi = res_atoms[iatom]
                
                xn = xi + cx 
                yn = yi + cy 
                zn = zi 

                ic = int (math.floor (xn/5))
                jc = int (math.floor (yn/5))
                kc = int (math.floor (zn/5))
                
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
            
                        dx = xn - xj
                        dy = yn - yj
                        dz = zn - zj
                        
                        r2 = dx*dx + dy*dy + dz*dz
                        
                        if (r2 < 10.0):
                            l_overlap = 1
                            break

                
            print 'ires is overlap ', l_overlap
           
            if (l_overlap == 0):
                   
                n_mem += 1
                str_resNM = "%4d"%n_mem
                
                for iatom in range (121):
                    line        = l_mem[iatom]
                    xi, yi, zi  = res_atoms[iatom]
                    xi         += cx
                    yi         += cy
                    
                    str_xyz = "%8.3f%8.3f%8.3f"%(xi, yi, zi)

                    line_new = line[0:22]+str_resNM+'    '+str_xyz+line[54:]
                    f_pdb.write(line_new + '\n')
                    
