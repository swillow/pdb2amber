import numpy as np
import sys
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd


def get_CA_positions(to_pdb_fname, from_pdb_fname):

    l_to = open(to_pdb_fname).readlines()
    l_from = open(from_pdb_fname).readlines()

    sel_to_ca = {}
    sel_from_ca = {}

    for line in l_to:
        if line[:6] in ['ATOM  ']:
            atName = line[12:16].strip()
            if atName in ['CA']:
                resID = line[17:26]

                words = line[30:].split()
                x = float(words[0])
                y = float(words[1])
                z = float(words[2])
                sel_to_ca[resID] = [x, y, z]

    sel_to_pos = []
    sel_from_pos = []

    for line in l_from:
        if line[:6] in ['ATOM  ']:
            atName = line[12:16].strip()
            if atName in ['CA']:
                resID = line[17:26]

                words = line[30:].split()
                x = float(words[0])
                y = float(words[1])
                z = float(words[2])

                if resID in sel_to_ca:
                    sel_to_pos.append(sel_to_ca[resID])
                    sel_from_pos.append([x, y, z])

    return np.array(sel_to_pos), np.array(sel_from_pos)


if __name__ == "__main__":
    import getopt

    argv = sys.argv[1:]
    opts, args = getopt.getopt(
        argv, "ht:f:s:", ["help=", "to_file=", "from_file=", "save_file="])

    if len(opts) == 0:
        print('python align_two_proteins.py -t <to_pdb_file> -f <from_pdb_file> -s <save_pdb_file>')
        sys.exit(2)

    for opt, arg in opts:

        if opt in ("-t", "--to_file"):
            to_pdb_file = arg
        elif opt in ("-f", "--from_file"):
            from_pdb_file = arg
        elif opt in ("-s", "--save_file"):
            save_pdb_file = arg
        elif opt in ('-h', "--help"):
            print(
                'python align_two_proteins.py -t <to_pdb_file> -f <from_pdb_file> -s <save_pdb_file>')
            sys.exit(1)

    to_pos, from_pos = get_CA_positions(to_pdb_file, from_pdb_file)

    to_com = to_pos.sum(axis=0)/to_pos.shape[0]
    from_com = from_pos.sum(axis=0)/from_pos.shape[0]

    to_pos -= to_com
    from_pos -= from_com

    R, rmsd = align.rotation_matrix(to_pos, from_pos)

    prt = mda.Universe(from_pdb_file)
    prt.atoms.translate(-from_com)
    prt.atoms.rotate(R)
    prt.atoms.translate(to_com)
    prt.atoms.write(save_pdb_file)