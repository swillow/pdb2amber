import sys
import numpy as np

def print_solvent_pointers  (data, fout):
    key_words = 'SOLVENT_POINTERS'
    str_flag   = '%-80s'%('%FLAG '+key_words)+'\n'
    str_format = '%-80s'%('%FORMAT(3I8)')+'\n'
    fout.write (str_flag)
    fout.write (str_format)
    line = '%8d%8d%8d'%(data[0], data[1], data[2]) + '\n'
    fout.write (line)

def print_string_radius_set (fout):
    key_words = 'RADIUS_SET'
    str_flag   = '%-80s'%('%FLAG '+key_words)+'\n'
    str_format = '%-80s'%('%FORMAT(1a80)')+'\n'
    fout.write (str_flag)
    fout.write (str_format)
    line = '%-80s'%('H(N)-modified Bondi radii (mbondi2)') + '\n'
    fout.write (line)

def print_string_data (key_words, data, fout):

    str_flag   = '%-80s'%('%FLAG '+key_words)+'\n'
    str_format = '%-80s'%('%FORMAT(20a4)')+'\n'

    fout.write (str_flag)
    fout.write (str_format)

    num = len(data)
    for ii in range (0, num, 20):
        if ii + 20 < num:
            str = ''
            for jj in range(20):
                val = data[ii+jj]
                str += '%-4s'%val
            line = str + '\n'
            fout.write(line)
        else:
            str = ''
            len0 = num - ii
            for jj in range (len0):
                val = data[ii+jj]
                str += '%-4s'%val
            line='%-80s'%str+'\n'
            fout.write(line)


def print_integer_data (key_words, data, fout):

    str_flag   = '%-80s'%('%FLAG '+key_words)+'\n'
    str_format = '%-80s'%('%FORMAT(10I8)')+'\n'

    fout.write (str_flag)
    fout.write (str_format)

    num = len(data)
    for ii in range (0, num, 10):
        if ii + 10 < num:
            str = ''
            for jj in range(10):
                val = data[ii+jj]
                str += '%8d'%val
            line = str+'\n'
            fout.write(line)
        else:
            str = ''
            len0 = num - ii
            for jj in range (len0):
                val = data[ii+jj]
                str += '%8d'%val
            line = '%-80s'%str+'\n'
            fout.write(line)


def print_double_data (key_words, data, fout):

    str_flag   = '%-80s'%('%FLAG '+key_words)+'\n'
    str_format = '%-80s'%('%FORMAT(5E16.8)')+'\n'

    fout.write (str_flag)
    fout.write (str_format)

    num = len(data)
    if num == 0:
        line = '\n'
        fout.write(line)
    else:
        for ii in range (0, num, 5):
            if ii + 5 < num:
                str = ''
                for jj in range(5):
                    val = data[ii+jj]
                    str += '%16.8E'%val
                line = str+'\n'
                fout.write(line)
            else:
                str = ''
                len0 = num - ii
                for jj in range (len0):
                    val = data[ii+jj]
                    str += '%16.8E'%val
                line = '%-80s'%str+'\n'
                fout.write(line)



class PrmTop (object):

    def __init__ (self):

        self.atom_name_list = []
        self.atom_type_list = []
        self.atomic_number_list = []
        self.atom_type_index_list = []
        self.amber_atom_type_list = []
        self.res_name_list  = []
        self.res_ptr_list   = []

        self.chg_list  = []
        self.mass_list = []
        self.radii_list = []

        self.bond_k_list  = []
        self.bond_length_list = []
        self.bond_wH_list = []
        self.bond_nH_list = []

        self.angle_k_list  = []
        self.angle_length_list = []
        self.angle_wH_list = []
        self.angle_nH_list = []

        self.proper_wH_list = []
        self.proper_nH_list = []
        self.proper_periodicity_list = []
        self.proper_phase_list = []
        self.proper_k_list     = []
#        self.scee_list = []
#        self.scnb_list = []

        self.solty_list = [0.0, 0.0]
        self.lj_acoef_list = []
        self.lj_bcoef_list = []
        self.nb_idx_list   = []

        self.num_excluded_atoms  = []
        self.excluded_atoms_list = []

        self.is_box   = 0
        self.box_info = []
        self.sol_ptr  = [0, 1, 2] # The last residue number of the solute, total number of molecules, the first solvent "molecule"
        self.atoms_per_molecule = []

        self.max_res_natom = 100

    def write_pointers(self, key_words, fout):

        str_flag   = '%-80s'%('%FLAG '+key_words)+'\n'
        str_format = '%-80s'%('%FORMAT(10I8)')+'\n'

        fout.write(str_flag)
        fout.write(str_format)

        data = []
        data.append(len(self.atom_name_list))
        data.append(len(self.atom_type_list))
        data.append(len(self.bond_wH_list)/3)
        data.append(len(self.bond_nH_list)/3)
        data.append(len(self.angle_wH_list)/4)
        data.append(len(self.angle_nH_list)/4)
        data.append(len(self.proper_wH_list)/5)
        data.append(len(self.proper_nH_list)/5)
        data.append(0) # NHPARM: NOT USED
        data.append(0) # NPARM:  NOT USED
        nnb = sum(self.num_excluded_atoms)
        data.append(nnb)
        data.append(len(self.res_name_list))
        data.append(len(self.bond_nH_list)/3)
        data.append(len(self.angle_nH_list)/4)
        data.append(len(self.proper_nH_list)/5)
        data.append(len(self.bond_k_list))
        data.append(len(self.angle_k_list))
        data.append(len(self.proper_k_list))
        data.append(len(self.solty_list)) # NATYP: NOT USED
        data.append(0) # NPHB: number of distinct 10-12 HB pair types
        data.append(0) # IFPERT: set to 1 if perturbation info is to be read in
        data.append(0) # NBPER
        data.append(0) # NGPER
        data.append(0) # NDPER
        data.append(0) # MBPER
        data.append(0) # MGPER
        data.append(0) # MDPER
        data.append(self.is_box) # IFBXO: 1 if standard periodic box
        data.append(self.max_res_natom) #
        data.append(0) # ifcap : CAP option
        data.append(0) # numextra
        data.append(0) # ncopy

        num = len(data)
        for ii in range (0, num, 10):
            if ii + 10 < num:
                str = ''
                for jj in range(10):
                    val = data[ii+jj]
                    str += '%8d'%val
                line = str+'\n'
                fout.write(line)
            else:
                str = ''
                len0 = num - ii
                for jj in range (len0):
                    val = data[ii+jj]
                    str += '%8d'%val
                line = '%-80s'%str+'\n'
                fout.write(line)



    def write (self, fname):

        fout = open (fname, 'w')

        head = '%-80s'%'%VERSION  VERSION_STAMP = V0001.000  DATE = 02/26/19  15:36:56'+'\n'
        head += '%-80s'%'%FLAG TITLE'   + '\n'
        head += '%-80s'%'%FORMAT(20a4)' + '\n'
        head += '%-80s'%'default_name'  + '\n'
        fout.write(head)


        self.write_pointers('POINTERS', fout)

        #f_C.write('%FLAG POINTERS\n')
        #f_C.write('%FORMAT(10I8) \n')

        print_string_data ('ATOM_NAME', self.atom_name_list, fout)
        print_double_data ('CHARGE', self.chg_list, fout)
        tot_chg = 0.0
        for chg in self.chg_list:
            tot_chg += chg
        print ('tot_chg ', tot_chg/18.2223)

        print_integer_data ('ATOMIC_NUMBER', self.atomic_number_list, fout)
        print_double_data ('MASS',   self.mass_list, fout)
        print_integer_data ('ATOM_TYPE_INDEX', self.atom_type_index_list, fout)
        print_integer_data ('NUMBER_EXCLUDED_ATOMS', self.num_excluded_atoms, fout)

        print_integer_data ('NONBONDED_PARM_INDEX', self.nb_idx_list, fout)
        print_string_data ('RESIDUE_LABEL', self.res_name_list, fout)
        print_integer_data('RESIDUE_POINTER', self.res_ptr_list, fout)
        print_double_data ('BOND_FORCE_CONSTANT', self.bond_k_list, fout)
        print_double_data ('BOND_EQUIL_VALUE', self.bond_length_list, fout)
        print_double_data ('ANGLE_FORCE_CONSTANT', self.angle_k_list, fout)
        print_double_data ('ANGLE_EQUIL_VALUE', self.angle_length_list, fout)
        print_double_data ('DIHEDRAL_FORCE_CONSTANT', self.proper_k_list, fout)
        print_double_data ('DIHEDRAL_PERIODICITY', self.proper_periodicity_list, fout)
        print_double_data ('DIHEDRAL_PHASE', self.proper_phase_list, fout)
        scee_list = np.ones ( (len(self.proper_phase_list)) )*1.2
        print_double_data ('SCEE_SCALE_FACTOR', scee_list, fout)
        scnb_list = np.ones ( (len(self.proper_phase_list)) )*2.0
        print_double_data ('SCNB_SCALE_FACTOR', scnb_list, fout)
        print_double_data ('SOLTY', self.solty_list, fout)

        print_double_data ('LENNARD_JONES_ACOEF', self.lj_acoef_list, fout)
        print_double_data ('LENNARD_JONES_BCOEF', self.lj_bcoef_list, fout)

        print_integer_data('BONDS_INC_HYDROGEN', self.bond_wH_list, fout)
        print_integer_data('BONDS_WITHOUT_HYDROGEN', self.bond_nH_list, fout)
        print_integer_data('ANGLES_INC_HYDROGEN', self.angle_wH_list, fout)
        print_integer_data('ANGLES_WITHOUT_HYDROGEN', self.angle_nH_list, fout)
        print_integer_data('DIHEDRALS_INC_HYDROGEN', self.proper_wH_list, fout)
        print_integer_data('DIHEDRALS_WITHOUT_HYDROGEN', self.proper_nH_list, fout)

        print_integer_data('EXCLUDED_ATOMS_LIST', self.excluded_atoms_list, fout)
        data = []
        print_double_data ('HBOND_ACOEF',   data, fout)
        print_double_data ('HBOND_BCOEF',   data, fout)
        print_double_data ('HBCUT',   data, fout)
        print_string_data ('AMBER_ATOM_TYPE', self.amber_atom_type_list, fout)

        data = []
        for ii in range(len(self.atomic_number_list)):
            data.append('BLA')
        print_string_data('TREE_CHAIN_CLASSIFICATION', data, fout)

        data = []
        for ii in range(len(self.atomic_number_list)):
            data.append(0)
        print_integer_data('JOIN_ARRAY', data, fout)
        print_integer_data('IROTAT', data, fout)

        if self.is_box == 1:
            print_solvent_pointers (self.sol_ptr, fout)
            if len (self.atoms_per_molecule) == 0:
                self.atoms_per_molecule.append (len (self.atom_name_list))
            print_integer_data('ATOMS_PER_MOLECULE', self.atoms_per_molecule, fout)
            print_double_data('BOX_DIMENSIONS', self.box_info, fout)

        print_string_radius_set(fout)
        print_double_data('RADII', self.radii_list, fout)

        data = []
        for znum in self.atomic_number_list:
            if znum == 1:
                data.append(0.85)
            elif znum == 6:
                data.append(0.72)
            elif znum == 7:
                data.append(0.79)
            elif znum == 8:
                data.append(0.85)
            elif znum == 9:
                data.append(0.88)
            elif znum == 15:
                data.append(0.86)
            elif znum == 16:
                data.append(0.96)
            elif znum == 26:
                data.append(0.88)
            else:
                data.append(0.80)


        print_double_data('SCREEN', data, fout)


        data = [0]
        print_integer_data('IPOL', data, fout)

        fout.close()
