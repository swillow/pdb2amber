import sys


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
    
