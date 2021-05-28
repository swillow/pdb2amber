import getopt
import sys


def pqr2pdb(fname_pqr, fname_pdb):

    lines = open(fname_pqr, 'r').readlines()
    fout = open(fname_pdb, 'w')

    chain_id = ['A', 'B', 'C', 'D', 'E', 'F',
                'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']
    ich = 0
    for line in lines:
        if line[:6] in ['ATOM  ', 'HETATM']:
            atnm = line[13:15].upper()
            resName = line[17:20].upper()
            symbol = atnm[0]
            if atnm[:2] in ['CL', 'NA', 'MG', 'BE', 'LI', 'ZN']:
                symbol = atnm[:2]
            elif atnm[:2] == 'CA' and resName[:2] == 'CA':
                symbol = 'CA'
            sline = line[:21] + chain_id[ich] + line[22:54] + \
                "  1.00  0.00          %2s  " % symbol
            print(sline, file=fout)
        elif line[:3] == 'END':
            print('END', file=fout)
        else:
            print(line[:-1], file=fout)

        if line[:3] == 'TER':
            ich += 1

    fout.close()


if __name__ == "__main__":

    argv = sys.argv[1:]
    opts, args = getopt.getopt(argv, "h:i:o:", ["help=", "input=", "output="])

    if len(opts) == 0:
        print('python pqr2pdb.py -i <pqr file> -o <pdb file>')
        sys.exit(2)

    fname_pdb = 'output.pdb'
    fname_pqr = 'input.pqr'
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('python pqr2pdb.py -i <pqr file> -o <pdb file>')
            sys.exit(1)
        if opt in ('-i', '--input'):
            fname_pqr = arg
        if opt in ('-o', '--output'):
            fname_pdb = arg

    pqr2pdb(fname_pqr, fname_pdb)
