import sys


def __get_line(data, iatom):
    """
    Generate the string based on data[iatom] and data[iatom+1]

    Parameters:
    data  : list (natom, 3), float
    iatom: intger
    """

    line = ""
    x, y, z = data[iatom]
    line += '%12.7f' % x
    line += '%12.7f' % y
    line += '%12.7f' % z

    if iatom+1 < len(data):
        x, y, z = data[iatom+1]
        line += '%12.7f' % x
        line += '%12.7f' % y
        line += '%12.7f' % z
    line += '\n'

    return line


def print_inpcrd(fname, crd, vel=None, box=None, time=None):
    """
    Save Amber Inpcrd format into fname

    Parameters
    ----------
    crd : list (natom, 3), float, coordinates
    vel : list (natom, 3), float, velocities
    box : list (3), float, box dimensions
    time: float, time in picosecond
    """

    fout = open(fname, 'w')

    title = "default_name\n"
    fout.write(title)
    natom = len(crd)
    if vel is None:
        line_natom = '%5d' % natom + '\n'
    else:
        line_natom = '%5d' % natom + '%15.7f' % time + '\n'
    fout.write(line_natom)

    for ii in range(0, natom, 2):
        line = __get_line(crd, ii)
        fout.write(line)

    if vel is not None:
        for ii in range(0, natom, 2):
            line = __get_line(vel, ii)
            fout.write(line)

    # Print Box Information
    if box is not None:
        angle = 90.0
        line = ""
        line += '%12.7f' % box[0]
        line += '%12.7f' % box[1]
        line += '%12.7f' % box[2]
        line += '%12.7f' % angle
        line += '%12.7f' % angle
        line += '%12.7f' % angle
        line += '\n'
        fout.write(line)

    fout.close()
