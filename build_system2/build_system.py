import os
import sys
import numpy as np


def read_prt_heavy_atoms(fname_prt, box_min, fout):
    """
    fname_prt : string, file name containing coordinates of a protein
    fout: output file
    return:
    numpy.array([x,y,z])
    """

    f_pdb = open(fname_prt, 'r')
    l_pdb = f_pdb.read().split('\n')
    f_pdb.close()

    prt_heavy_atoms = []
    for line in l_pdb[:-1]:
        if line[:3] == 'END':
            continue
        if line[:6] in ['ATOM  ', 'HETATM']:
            words = line[30:].split()
            sym = words[-1]
            x = float(words[0]) - box_min[0]
            y = float(words[1]) - box_min[1]
            z = float(words[2]) - box_min[2]
            sx = ("%8.3f" % x)[:8]
            sy = ("%8.3f" % y)[:8]
            sz = ("%8.3f" % z)[:8]

            line_new = line[0:30]+sx+sy+sz+line[54:]
            print(line_new, file=fout)

            if sym != 'H':
                prt_heavy_atoms.append([x, y, z])
        else:
            print(line, file=fout)

    prt_heavy_atoms = np.array(prt_heavy_atoms)

    return prt_heavy_atoms


def build_membrane(prt_heavy_atoms, fname_mem, box, shift_z, fout):
    """
    box: np.array (3) the system box 
    """

    prt_min = prt_heavy_atoms.min(axis=0)
    prt_max = prt_heavy_atoms.max(axis=0)

    ndim_min = np.array(prt_min//5, dtype=np.int)
    ndim_max = np.array(prt_max//5, dtype=np.int) + 1

    cell_list = {}

    for ic in range(ndim_min[0], ndim_max[0]):  # min_box[0]/5, max_box[0]/5
        for jc in range(ndim_min[1], ndim_max[1]):
            for kc in range(ndim_min[2], ndim_max[2]):
                cell_list[(ic, jc, kc)] = []

    neigh_list = []
    for ic in range(-1, 2):
        for jc in range(-1, 2):
            for kc in range(-1, 2):
                neigh_list.append([ic, jc, kc])

    neigh_list = np.array(neigh_list, dtype=np.int)

    for pi in prt_heavy_atoms:
        ic, jc, kc = np.array(pi//5, dtype=np.int)
        cell_list[(ic, jc, kc)].append(pi)
    # READ MEMBRANE
    f_mem = open(fname_mem)
    l_mem = f_mem.read().split('\n')
    f_mem.close()

    tatom = 0
    res_atoms_list = []
    #res_atoms = np.zeros((121, 3), dtype=np.float)
    minx = 100.0
    maxx = 0.0
    miny = 100.0
    maxy = 0.0
    minz = 100.0
    maxz = 0.0
    for ires in range(8):
        res_atoms = []
        for iatom in range(121):
            line = l_mem[tatom]
            words = line[30:].split()

            x = float(words[0])
            minx = min(x, minx)
            maxx = max(x, maxx)

            y = float(words[1])
            miny = min(y, miny)
            maxy = max(y, maxy)

            z = float(words[2])
            minz = min(z, minz)
            maxz = max(z, maxz)

            res_atoms.append([x, y, z])
            tatom += 1

        res_atoms_list.append(res_atoms)

    res_atoms_list = np.array(res_atoms_list)

    ndim_box = np.array(box//15, dtype=np.int)

    mem_heavy_atoms = []
    n_mem = 0
    # membrane is in a dimension [0:15, 0:15, -24:24]
    for ix in range(ndim_box[0]):
        cx = 15.0*ix
        for iy in range(ndim_box[1]):
            cy = 15.0*iy

            pc = np.array([cx, cy, shift_z])

            for ires in range(8):
                res_atoms = res_atoms_list[ires]

                l_overlap = 0

                for pi in res_atoms:

                    if l_overlap == 1:
                        break

                    pn = pi + pc

                    icel = np.array(pn//5, dtype=np.int)

                    for ncel in neigh_list:

                        if l_overlap == 1:
                            break

                        ic, jc, kc = icel + ncel

                        if (ic, jc, kc) in cell_list:
                            # crds of prt_heavy_atoms
                            for pj in cell_list[(ic, jc, kc)]:
                                pij = pn - pj
                                rij2 = (np.einsum('i,i', pij, pij))

                                if (rij2 < 10.0):
                                    l_overlap = 1
                                    break

                if l_overlap == 0:

                    n_mem += 1

                    for iatom, pi in enumerate(res_atoms):

                        line = l_mem[iatom]
                        xi, yi, zi = pi + pc

                        str_xyz = "%4d    %8.3f%8.3f%8.3f" % (
                            n_mem, xi, yi, zi)

                        line_new = line[0:22] + str_xyz+line[54:]
                        print(line_new, file=fout)

                        words = line[54:].split()
                        if words[-1] != 'H':
                            mem_heavy_atoms.append([xi, yi, zi])

    return np.array(mem_heavy_atoms)


def build_solution(prt_heavy_atoms, mem_heavy_atoms,
                   fname_wat, prt_chg, molarity, box, shift_z, fout):

    #
    perWat = int(1.0/(molarity*0.03345))

    heavy_atoms = np.concatenate((prt_heavy_atoms, mem_heavy_atoms), axis=0)
    ndim_min = np.array(heavy_atoms.min(axis=0)//5, dtype=np.int)
    ndim_max = np.array(heavy_atoms.max(axis=0)//5, dtype=np.int)+1

    cell_list = {}
    for ic in range(ndim_min[0], ndim_max[0]):
        for jc in range(ndim_min[1], ndim_max[1]):
            for kc in range(ndim_min[2], ndim_max[2]):
                cell_list[(ic, jc, kc)] = []

    for pi in heavy_atoms:
        ic, jc, kc = np.array(pi//5, dtype=np.int)
        if (ic, jc, kc) in cell_list:
            cell_list[(ic, jc, kc)].append(pi)
        else:
            print('missing cell_list in water build', ic, jc, kc)
    del heavy_atoms

    neigh_list = []
    for ic in range(-1, 2):
        for jc in range(-1, 2):
            for kc in range(-1, 2):
                neigh_list.append([ic, jc, kc])
    neigh_list = np.array(neigh_list, dtype=np.int)

    f_xyz = open(fname_wat)
    l_xyz = f_xyz.read().split('\n')
    f_xyz.close()

    nwat = int(l_xyz[0].split()[0])//3
    l_xyz.pop(0)  # natom
    l_xyz.pop(0)  # comment

    wo_list = np.zeros((nwat, 3))
    wh1_list = np.zeros((nwat, 3))
    wh2_list = np.zeros((nwat, 3))

    for iw in range(nwat):
        io = 3*iw
        ih1 = io+1
        ih2 = io+2

        words = l_xyz[io].split()
        wo_list[iw] = [float(words[1]), float(words[2]), float(words[3])]

        words = l_xyz[ih1].split()
        wh1_list[iw] = [float(words[1]), float(words[2]), float(words[3])]

        words = l_xyz[ih2].split()
        wh2_list[iw] = [float(words[1]), float(words[2]), float(words[3])]

    ###
    nres = 0
    natom = 0

    ndim_box = np.array(box//15, dtype=np.int)
    ion_list = []
    l_add_ion = True

    for ix in range(ndim_box[0]):
        cx = 15.0*ix
        for iy in range(ndim_box[1]):
            cy = 15.0*iy

            for iz in range(ndim_box[2]):
                cz = 15.0*iz

                cell0 = np.array([cx, cy, cz])

                for iw in range(nwat):

                    po = wo_list[iw] + cell0
                    if -18 + shift_z < po[2] and po[2] < 18 + shift_z:
                        # where the mebrane is placed
                        continue

                    icel = np.array(po//5, dtype=np.int)

                    l_overlap = 0
                    for ncel in neigh_list:
                        if l_overlap == 1:
                            break

                        ic2, jc2, kc2 = icel + ncel

                        if (ic2, jc2, kc2) not in cell_list:
                            continue

                        for pj in cell_list[(ic2, jc2, kc2)]:
                            poj = po - pj
                            r2 = np.einsum('i,i', poj, poj)

                            if (r2 < 16.0):
                                l_overlap = 1
                                break

                    if l_overlap == 0:
                        xo, yo, zo = po
                        if l_add_ion and nres % perWat == 0:
                            ion_list.append([xo, yo, zo])
                            l_add_ion = False
                        else:
                            natom += 1
                            nres += 1
                            l_add_ion = True

                            line = "HETATM" + "%5d" % (natom % 100000) + "  O   HOH " + \
                                " %4d    " % (nres % 10000) + \
                                "%8.3f%8.3f%8.3f" % (xo, yo, zo) + \
                                "  1.00  0.00           O"
                            print(line, file=fout)

                            natom += 1
                            xh, yh, zh = wh1_list[iw] + cell0
                            line = "HETATM" + "%5d" % (natom % 100000) + "  H1  HOH " + \
                                " %4d    " % (nres % 10000) + \
                                "%8.3f%8.3f%8.3f" % (xh, yh, zh) + \
                                "  1.00  0.00           H"
                            print(line, file=fout)

                            natom += 1
                            xh, yh, zh = wh2_list[iw] + cell0
                            line = "HETATM" + "%5d" % (natom % 100000) + "  H2  HOH " + \
                                " %4d    " % (nres % 10000) + \
                                "%8.3f%8.3f%8.3f" % (xh, yh, zh) + \
                                "  1.00  0.00           H"
                            print(line, file=fout)

    ion_list = np.array(ion_list)
    numIon = (len(ion_list)-abs(prt_chg))//2
    idx_Cl = [ii for ii in range(0, 2*numIon, 2)]
    idx_Na = [ii for ii in range(1, 2*numIon, 2)]

    if prt_chg < 0:
        # to be neutral
        # Na- > Cl-
        idx_Na += [ii for ii in range(2*numIon, 2*numIon+abs(prt_chg))]
    else:
        # to be neutral of the system
        #  Cl- > Na+
        idx_Cl += [ii for ii in range(2*numIon, 2*numIon+prt_chg)]

    natom = 0
    nres = 0
    for ii in idx_Na:
        natom += 1
        nres += 1
        x, y, z = ion_list[ii]
        line = "HETATM" + "%5d" % natom + "  NA   NA " + \
            " %4d    " % nres + \
            "%8.3f%8.3f%8.3f" % (x, y, z) + \
            "  1.00  0.00          Na"
        print(line, file=fout)
    natom = 0
    nres = 0
    for ii in idx_Cl:
        natom += 1
        nres += 1
        x, y, z = ion_list[ii]
        line = "HETATM" + "%5d" % natom + "  CL   CL " + \
            " %4d    " % nres + \
            "%8.3f%8.3f%8.3f" % (x, y, z) + \
            "  1.00  0.00          Cl"
        print(line, file=fout)


if __name__ == '__main__':

    import getopt
    """
    1. Build membrane with DPPE
    2. Build Solution with 0.1 M NaCl
    (Molarity = # of NaCl/ (0.03345* # of Water), 
     where 0.03345 is the number density of water (1 g/cm^{-3}))
    """

    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "hi:b:c:s:",
                               ["help=", "input=", "box=", "charge=", "save="])

    if len(opts) == 0:
        print('python build_system.py -i <opm pdb file> -b "box lengths" -c <charge> -s <save_pdb_file>')
        sys.exit(2)
    #    box_min = np.array([-60.0, -60.0, -90.0])
    #    box_max = np.array([60.0, 60.0, 75.0])
    box_min = np.array([-75.0, -75.0, -75.0])
    box_max = np.array([75.0, 75.0, 75.0])
    box = box_max - box_min
    prt_chg = 0  # -27
    molarity = 0.1
    fname_prt = "opm.pdb"
    fname_sys = "system.pdb"
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(
                'python build_system.py -i <opm pdb file> -b "box lengths" -c <charge> -s <save_pdb_file>')
            sys.exit(1)
        elif opt in ('-i', '--input'):
            fname_prt = arg
        elif opt in ('-b', '--box'):
            bl = arg.split()
            box = np.array([float(bl[0]), float(bl[1]), float(bl[2])])
            box_min = -0.5*box
            box_max = 0.5*box
        elif opt in ('-c', '--charge'):
            prt_chg = int(arg)
        elif opt in ('-s', '--save'):
            fname_sys = arg

    fout = open(fname_sys, 'w', 1)

    alpha = 90.0
    beta = 90.0
    gamma = 90.0
    print("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1 1 " % (
        box[0], box[1], box[2], alpha, beta, gamma), file=fout)
    fname_prt = 'rib_nqr_opm.pdb'
    prt_heavy_atoms = read_prt_heavy_atoms(fname_prt, box_min, fout)

    #
    print("Build Membrane")
    fname_mem = 'dppe_box15.pdb'
    shift_z = -box_min[2]
    mem_heavy_atoms = build_membrane(prt_heavy_atoms, fname_mem,
                                     box, shift_z, fout)

    print("Build Solution (0.1 M NaCl)")
    fname_wat = 'water_box15.xyz'
    build_solution(prt_heavy_atoms, mem_heavy_atoms, fname_wat,
                   prt_chg, molarity,
                   box, shift_z, fout)

    fout.close()
