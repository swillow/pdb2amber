import sys, getopt
import simtk.unit as units

from _pdbfile import *

def print_inpcrd(fname, crd):

    fout = open (fname, 'w')

    title = "default_name\n"
    fout.write (title)
    natom = len(crd)
    str_natom = '%5d'%natom +'\n'
    fout.write (str_natom)

    for ii in range (0, natom, 2):
        if ii + 2 < natom:
            str = ''
            x, y, z = crd[ii]
            str += '%12.7f'%x
            str += '%12.7f'%y
            str += '%12.7f'%z
            x, y, z = crd[ii+1]
            str += '%12.7f'%x
            str += '%12.7f'%y
            str += '%12.7f'%z
            str += '\n'
            fout.write (str)
        else:
            str = ''
            x, y, z = crd[ii]
            str += '%12.7f'%x
            str += '%12.7f'%y
            str += '%12.7f'%z

            if natom - ii == 2:
                x, y, z = crd[ii+1]
                str += '%12.7f'%x
                str += '%12.7f'%y
                str += '%12.7f'%z
                
            str += '\n'
            fout.write(str)
                
    fout.close()


if __name__ == "__main__":

    pdb_fname = ''
    inpcrd_fname = ''

    argv = sys.argv[1:]

    opts, args = getopt.getopt (argv, "hi:o:", ["help=", "ifile=", "ofile="])

    if (len (opts) == 0):
        print ("pdb2inpcrd.py -i <input.pdb> -o <output.inpcrd>")
        sys.exit(1)
        
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ("pdb2inpcrd.py -i <input.pdb> -o <output.inpcrd>")
            sys.exit(2)
            
        elif opt in ("-i", "--ifile"):
            pdb_fname = arg
        elif opt in ("-o", "--ofile"):
            inpcrd_fname = arg

    pdb = MyPDBFile (pdb_fname)
    
    print_inpcrd (inpcrd_fname, pdb.positions.value_in_unit (units.angstroms))
