import os
import itertools
import sys
import math
import numpy as np
import xml.etree.ElementTree as etree
from copy import copy
from collections import defaultdict
import simtk.openmm as mm
from datetime import date
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.openmm.app.internal import amber_file_parser
from simtk.openmm.app.internal.unitcell import computeLengthsAndAngles

import simtk.unit as units
from simtk.openmm.app import element as elem
from decimal import Decimal
import sys, getopt

from _pdbfile import *
from _prmtop import *
from _forcefield import *
from _inpcrd import *

       
class MyTorsion (object):

    def __init__ (self, index, idx_torsion):
        # torsion: class PeriodicTorsion
       self.index       = index 
       self.idx_torsion = idx_torsion

       
def get_proper_typeID(proper_database, type1, type2, type3, type4):

    typeID = 'PR_' + type1 + '_' + type2 + '_' + type3 + '_' + type4

    if typeID in proper_database:
        return typeID
    
    typeID = 'PR_' + type4 + '_' + type3 + '_' + type2 + '_' + type1
    if typeID in proper_database:
        return typeID

    typeID = 'PR_X_' + type2 + '_' + type3 + '_X'
    if typeID in proper_database:
        return typeID
    
    typeID = 'PR_X_' + type3 + '_' + type2 + '_X'
    if typeID in proper_database:
        return typeID

    return None


def get_improper_typeID(improper_database, type1, type2, type3, type4):

    typeID = 'IM_' + type1 + '_' + type2 + '_' + type3 + '_' + type4
    if typeID in improper_database:
        return typeID
    
    typeID = 'IM_' + type2 + '_' + type1 + '_' + type3 + '_' + type4
    if typeID in improper_database:
        return typeID
#debug 
#    typeID = 'IM_' + type2 + '_' + type4 + '_' + type3 + '_' + type1
#    if typeID in improper_database:
#        return typeID
    
#    typeID = 'IM_' + type4 + '_' + type2 + '_' + type3 + '_' + type1
#    if typeID in improper_database:
#        return typeID
    
#    typeID = 'IM_' + type1 + '_' + type4 + '_' + type3 + '_' + type2
#    if typeID in improper_database:
#        return typeID
    
#    typeID = 'IM_' + type4 + '_' + type1 + '_' + type3 + '_' + type2
#    if typeID in improper_database:
#        return typeID
    
    typeID = 'IM_X_' + type1 + '_' + type3 + '_' + type4
    if typeID in improper_database:
        return typeID

    typeID = 'IM_X_' + type2 + '_' + type3 + '_' + type4
    if typeID in improper_database:
        return typeID

    typeID = 'IM_X_X_' + type3 + '_' + type4
    if typeID in improper_database:
        return typeID

    typeID = 'IM_X_X_' + type3 + '_' + type2
    if typeID in improper_database:
        return typeID

    typeID = 'IM_X_X_' + type3 + '_' + type1
    if typeID in improper_database:
        return typeID

    return None


def pdb2amber (pdb_fname, prmtop_fname, inpcrd_fname, ff_fnames):

    pdb   = MyPDBFile(pdb_fname)
    my_ff = MyForceFields (ff_fnames)

    _atomType  = {}
    _excludedAtomWith = []
    _bonds     = []
    _angles    = []
    _propers   = []
    _impropers = []
    _bondedToAtom = []

    _atoms = list(pdb.topology.atoms())

    prm = PrmTop()

    vectors = pdb.topology.getPeriodicBoxVectors()
    if vectors is not None:
        prm.is_box = 1
        a, b, c, alpha, beta, gamma = computeLengthsAndAngles(vectors)
        # radian-->degree, nm --> A
        prm.box_info = [ beta*180.0/np.pi, a*10.0, b*10.0, c*10.0] 
        
    if inpcrd_fname != '':
        print_inpcrd (inpcrd_fname,
                      pdb.positions.value_in_unit (units.angstroms))
    #for debuging
    #prm_amber = amber_file_parser.PrmtopLoader('receptor.prmtop')
    
    # Gather Atom Info
    for atom in _atoms:
        _excludedAtomWith.append([])
        atName = atom.name
        atElem = atom.element
        prm.atom_name_list.append(atom.name)
        prm.atomic_number_list.append(atElem._atomic_number)

    # Find the Residue Template maching each Residue
    unmatchedResidues = []
    at_type_list  = []
    max_res_natom = 0
    is_solvent    = ["DPPE", "DPP", "HOH", "WAT"]
    
    for chain in pdb.topology.chains():
        nres = chain._residues[0];
        cres = chain._residues[-1];
        
        for res in chain.residues():
            prm.res_name_list.append(res.name)

            resName = res.name


                
#            if resName in MyPDBFile._residueNameReplacements:
#                resName = MyPDBFile._residueNameReplacements[resName]
                
            atomReplacements = {}
            if nres == res: 
                res.name = "N"+resName
                if res.name not in my_ff._residues:
                    res.name = resName
                atomReplacements = MyPDBFile._atomNameReplacements["Protein"]
                
            if cres == res:
                res.name = "C"+resName
                if res.name not in my_ff._residues:
                    res.name = resName
                atomReplacements = MyPDBFile._atomNameReplacements["Protein"]
                
            if resName in MyPDBFile._atomNameReplacements:
                ff_res_atomNameRepl = MyPDBFile._atomNameReplacements[resName]
                for ff_atName in ff_res_atomNameRepl:
                    atomReplacements[ff_atName] = ff_res_atomNameRepl[ff_atName]
            else: 
                print res.name, 'has no atomNameReplacement'

            
            if res.name in my_ff._residues:
                ff_resData = my_ff._residues[res.name]
                atoms = list(res.atoms())
                natom = len(atoms)
                natom_ff = len(ff_resData.atoms)
                if max_res_natom < natom:
                    max_res_natom = natom
                if natom != natom_ff:
                    print ("%-5s"%res.name + "has %d"%natom + " but Forcefield has %d"%natom_ff)
                    
                prm.res_ptr_list.append (atoms[0].index+1)
                
                if resName in is_solvent:
                    if prm.sol_ptr[1] == 1:
                        prm.sol_ptr[0] = len(prm.res_name_list) - 1
                        prm.atoms_per_molecule.append (atoms[0].index)
                    prm.sol_ptr[1] += 1
                    prm.atoms_per_molecule.append(natom)
                    
                for ii in range (natom):
                    atName = atoms[ii].name
                    imatch = -1
                    for jj in range (natom_ff):
                        ff_at_name = ff_resData.atoms[jj].at_name
                        
                        if atName == ff_at_name:
                            imatch = jj
                            break
                        
                        elif ff_at_name in atomReplacements:
                            ff_at_name = atomReplacements[ff_at_name]
                            if atName == ff_at_name:
                                imatch = jj
                                break
                    if imatch == -1:
#                        raise Exception("Could not identify atom '%s'"%atName + " at residue %s "%res.name)
                        print("Could not identify atom '%s'"%atName + " at residue %s "%res.name)
                        
                        for ii in range (natom):
                            print ("%5d"%(ii+1) + "%5s"%atoms[ii].name)
                        for jj in range (natom_ff):
                            print ("%5d"%(jj+1) + "%5s"%ff_resData.atoms[jj].at_name)
                                       
                    
                    typeName = ff_resData.atoms[imatch].at_type
                    if typeName not in at_type_list:
                        at_type_list.append(typeName)
                        
                    mass     = my_ff._atomTypes[typeName].mass 
                    _atomType[atoms[ii]] = typeName
                    at_chg = Decimal(ff_resData.atoms[imatch].at_chg)*Decimal(18.2223)
                    
                    prm.chg_list.append(at_chg)
                    prm.mass_list.append(mass)
            else: 
                raise Exception("No Residue2 '%s'."%res.name)


    prm.max_res_natom = max_res_natom        
    at_type_list = sorted (at_type_list)
    prm.atom_type_list = at_type_list
    numTypes = len(at_type_list)
                    
    for atom in _atoms:
        typeName = _atomType[atom]
        atomClass = my_ff._atomTypes[typeName].atomClass
        
        index = -1
        for ii in range(len (at_type_list)):
            if typeName == at_type_list[ii]:
                index = ii+1
                break
        prm.atom_type_index_list.append(index)
        prm.amber_atom_type_list.append(atomClass)

        
    for ii in range (numTypes):
        for jj in range (numTypes):
            prm.nb_idx_list.append(0)
            
    nbIdx = 0
    for ii in range(numTypes):
        typeName       = at_type_list[ii]
        eps_i, sigma_i = my_ff._nonbonds[typeName]
        for jj in range (ii+1):
            nbIdx = nbIdx+1
            
            idx = numTypes*ii+jj
            prm.nb_idx_list[idx] = nbIdx
            idx = numTypes*jj+ii
            prm.nb_idx_list[idx] = nbIdx
            
            typeName       = at_type_list[jj]
            eps_j, sigma_j = my_ff._nonbonds[typeName]
            
            eps_ij = Decimal(math.sqrt(eps_i*eps_j))
            sig_ij = Decimal(0.5)*(sigma_i + sigma_j)
            
            lj_B   = Decimal(4.0)*eps_ij*sig_ij**6
            lj_A   = Decimal(lj_B*sig_ij**6)
            prm.lj_acoef_list.append (lj_A)
            prm.lj_bcoef_list.append (lj_B)
            

    # BOND LIST
    _type_list = []
    
    for bond in pdb.topology.bonds():
        _bonds.append (BondData (bond[0].index, bond[1].index) )
        
        iatom = bond[0].index
        jatom = bond[1].index
            
        iatElem = bond[0].element
        jatElem = bond[1].element
        
        type1 = _atomType[bond[0]]
        type2 = _atomType[bond[1]]

        itype = -1
        if type1 < type2:
            if (type1, type2) not in _type_list:
                _type_list.append ((type1, type2))
                itype = len(_type_list)
            else:
                for ii in range (len(_type_list)):
                    if (type1, type2) == _type_list[ii]:
                        itype = ii+1
                        break
        else:        
            if (type2, type1) not in _type_list:
                _type_list.append ((type2, type1))
                itype = len(_type_list)
            else:
               for ii in range (len(_type_list)):
                   if (type2, type1) == _type_list[ii]:
                       itype = ii+1
                       break
                   
        if iatElem._atomic_number == 1 or \
           jatElem._atomic_number == 1:
            if iatom < jatom:
                prm.bond_wH_list.append(iatom*3)
                prm.bond_wH_list.append(jatom*3)
            else: 
                prm.bond_wH_list.append(jatom*3)
                prm.bond_wH_list.append(iatom*3)
            prm.bond_wH_list.append(itype)
        else:
            if iatom < jatom:
                prm.bond_nH_list.append(iatom*3)
                prm.bond_nH_list.append(jatom*3)
            else: 
                prm.bond_nH_list.append(jatom*3)
                prm.bond_nH_list.append(iatom*3)
            prm.bond_nH_list.append(itype)
        
    
    for (type1, type2) in _type_list:
        bond_k      = 0.0
        bond_length = 0.0
        l_found    = False 
        for jj in range (len(my_ff._harmonicBonds)):
            types1 = my_ff._harmonicBonds[jj].types1
            types2 = my_ff._harmonicBonds[jj].types2
            if (type1 == types1 and type2 == types2):
                bond_k      = my_ff._harmonicBonds[jj].k
                bond_length = my_ff._harmonicBonds[jj].length
                l_found = True
                break
        if not l_found:
            print 'Error: No defined Bond ', type1, ' ', type2
            
        prm.bond_k_list.append(bond_k)
        prm.bond_length_list.append(bond_length)

    ###### 
    for ii in range(len(_atoms)):
        _bondedToAtom.append(set())

    _bonds = sorted (_bonds)  # add willow
    for ii in range(len(_bonds)):
        bond = _bonds[ii]
        _bondedToAtom[bond.atom1].add(bond.atom2)
        _bondedToAtom[bond.atom2].add(bond.atom1)
        if bond.atom2 not in _excludedAtomWith[bond.atom1]:
            if bond.atom1 not in _excludedAtomWith[bond.atom2]:
                _excludedAtomWith[bond.atom1].append(bond.atom2)
        

    
    # Make a list of all unique angles
    uniqueAngles = set()
    for bond in _bonds:
        iatom = bond.atom1
        jatom = bond.atom2
        
        for katom in _bondedToAtom[iatom]:
            if katom != jatom:
                if katom < jatom:
                    uniqueAngles.add ( (katom, iatom, jatom) )
                else:
                    uniqueAngles.add ( (jatom, iatom, katom) )

        for katom in _bondedToAtom[jatom]:
            if katom != iatom:
                if katom < iatom:
                    uniqueAngles.add ( (katom, jatom, iatom) )
                else:
                    uniqueAngles.add ( (iatom, jatom, katom) )


    _angles = sorted (list(uniqueAngles))


    _type_list = []
    for angle in _angles:
        iatom = angle[0]
        jatom = angle[1]
        katom = angle[2]

        iatElem = _atoms[iatom].element
        katElem = _atoms[katom].element
        
        type1 = _atomType[_atoms[iatom]]
        type2 = _atomType[_atoms[jatom]]
        type3 = _atomType[_atoms[katom]]

        if katom not in _excludedAtomWith[iatom]:
            if iatom not in _excludedAtomWith[katom]:
                _excludedAtomWith[iatom].append(katom)
        
        itype = -1
        if type1 < type3:
            if (type1, type2, type3) not in _type_list:
                _type_list.append ((type1, type2, type3))
                itype = len (_type_list)
            else:
                for ii in range (len(_type_list)):
                    if (type1, type2, type3) == _type_list[ii]:
                        itype = ii+1
                        break
        else:
            if (type3, type2, type1) not in _type_list:
                _type_list.append ((type3, type2, type1))
                itype = len (_type_list)
            else:
                for ii in range (len(_type_list)):
                    if (type3, type2, type1) == _type_list[ii]:
                        itype = ii+1
                        
        if iatElem._atomic_number == 1 or \
           katElem._atomic_number == 1:
            if iatom < katom:
                prm.angle_wH_list.append(iatom*3)
                prm.angle_wH_list.append(jatom*3)
                prm.angle_wH_list.append(katom*3)
            else: 
                prm.angle_wH_list.append(katom*3)
                prm.angle_wH_list.append(jatom*3)
                prm.angle_wH_list.append(iatom*3)
            prm.angle_wH_list.append(itype)
        else:
            if iatom < katom:
                prm.angle_nH_list.append(iatom*3)
                prm.angle_nH_list.append(jatom*3)
                prm.angle_nH_list.append(katom*3)
            else: 
                prm.angle_nH_list.append(katom*3)
                prm.angle_nH_list.append(jatom*3)
                prm.angle_nH_list.append(iatom*3)
            prm.angle_nH_list.append(itype)

    # angle energy conversion factor

    for (type1, type2, type3) in _type_list:
        ang_k      = 0.0
        ang_length = 0.0
        l_found    = False 

        for jj in range (len (my_ff._harmonicAngles) ):
            types1 = my_ff._harmonicAngles[jj].types1
            types2 = my_ff._harmonicAngles[jj].types2
            types3 = my_ff._harmonicAngles[jj].types3

            if type1 == types1 and type2 == types2 and type3 == types3:
                ang_k      = my_ff._harmonicAngles[jj].k
                ang_length = my_ff._harmonicAngles[jj].angle
                l_found    = True
                break
        if not l_found:
            print 'Error: No defined Angle ', type1, ' ', type2, ' ', type3
            
        prm.angle_k_list.append(ang_k)
        prm.angle_length_list.append(ang_length)

        
    # Make a list of all unique proper torsions

    uniquePropers = set()
    
    for angle in _angles:
            
        for atom in _bondedToAtom[angle[0]]:
            if atom not in angle:
                if atom < angle[2]:
                    uniquePropers.add( (atom, angle[0], angle[1], angle[2]) )
                else:
                    uniquePropers.add( (angle[2], angle[1], angle[0], atom) )

        for atom in _bondedToAtom[angle[2]]:
            if atom not in angle:
                if atom < angle[0]:
                    uniquePropers.add( (atom, angle[2], angle[1], angle[0]) )
                else:
                    uniquePropers.add( (angle[0], angle[1], angle[2], atom) )

    _propers = sorted (list(uniquePropers))
    _type_dict = {}
    _unique_torsion_dict = {}
    _dihedral_index = 0
    
    for proper in _propers:
        iatom = proper[0]
        jatom = proper[1]
        katom = proper[2]
        latom = proper[3]
        
        type1 = _atomType[_atoms[iatom]]
        type2 = _atomType[_atoms[jatom]]
        type3 = _atomType[_atoms[katom]]
        type4 = _atomType[_atoms[latom]]

        vsign  = -1.0
        if latom not in _excludedAtomWith[iatom]:
            if iatom not in _excludedAtomWith[latom]:
                _excludedAtomWith[iatom].append(latom)
                vsign = 1.0

        typeID = get_proper_typeID (my_ff._propers, type1, type2, type3, type4)

        if typeID == None:
            print 'Error: No defined Proper Torsion ', type1, ' ', type2, ' ', type3, ' ', type4
            continue
        
        iatElem = _atoms[iatom].element
        latElem = _atoms[latom].element


        if typeID not in _type_dict:
            
            idx         = my_ff._propers[typeID].index_torsion
            periodicity = my_ff._unique_torsion_list[idx].periodicity
            phase       = my_ff._unique_torsion_list[idx].phase
            kval        = my_ff._unique_torsion_list[idx].k

            if idx not in _unique_torsion_dict:
                _unique_torsion_dict[idx]    = _dihedral_index 
                prm.proper_periodicity_list += periodicity
                prm.proper_phase_list       += phase
                prm.proper_k_list           += kval
                _dihedral_index             += len (phase)
                
            _type_dict[typeID] = MyTorsion(_unique_torsion_dict[idx], idx)
            
            for ii in range(len(phase)):
                prm.scee_list.append (Decimal(1.2))
                prm.scnb_list.append (Decimal(2.0))
            
        itype     = _type_dict[typeID].index + 1
        idx       = _type_dict[typeID].idx_torsion
        
        if iatElem._atomic_number == 1 or \
           latElem._atomic_number == 1:

            for jj in range (len (my_ff._unique_torsion_list[idx].phase)):
                
                prm.proper_wH_list.append(iatom*3)
                prm.proper_wH_list.append(jatom*3)
                if jj == 0:
                    prm.proper_wH_list.append(vsign*katom*3)
                else:
                    prm.proper_wH_list.append(-katom*3)
                prm.proper_wH_list.append(latom*3)
                prm.proper_wH_list.append(itype+jj)
        else:
            
            for jj in range (len (my_ff._unique_torsion_list[idx].phase)):
                prm.proper_nH_list.append(iatom*3)
                prm.proper_nH_list.append(jatom*3)
                if jj == 0:
                    prm.proper_nH_list.append(vsign*katom*3)
                else:
                    prm.proper_nH_list.append(-katom*3)
                prm.proper_nH_list.append(latom*3)
                prm.proper_nH_list.append(itype+jj)
                
    # ---- PASS proper_list ----


    #
    # Make a list of all unique improper torsions
    for iatom in range (len(_bondedToAtom)):
        bondedTo = _bondedToAtom[iatom]
        
        if len(bondedTo) == 3:
            subset = sorted(list(bondedTo))
            _impropers.append ((subset[0], subset[1], iatom, subset[2]))

        
    debug = [19, 23, 21, 22] 
    _impropers = sorted(_impropers)
    _dbg_improper_list = []
    for improper in _impropers:
        iatom = improper[0]
        jatom = improper[1]
        katom = improper[2]
        latom = improper[3]
        
        iatType = _atomType[_atoms[iatom]]
        jatType = _atomType[_atoms[jatom]]
        katType = _atomType[_atoms[katom]]
        latType = _atomType[_atoms[latom]]

#        if  katType == 'protein-N' and _atoms[katom].residue.name != 'PRO':
#            continue
        
        iatElem = _atoms[iatom].element
        jatElem = _atoms[jatom].element
        latElem = _atoms[latom].element

        typeID = get_improper_typeID (my_ff._impropers,
                                      iatType, jatType, katType, latType)

        if typeID == None:
            if jatType != latType:
                (jatType, latType) = (latType, jatType)
                (jatom,   latom)   = (latom,   jatom)
                typeID = get_improper_typeID (my_ff._impropers,
                                              iatType, jatType, katType, latType)
                
        if typeID == None:
            if iatType != latType:
                (iatType, latType) = (latType, iatType)
                (iatom,   latom)   = (latom,   iatom)
                typeID = get_improper_typeID (my_ff._impropers,
                                              iatType, jatType, katType, latType)
            
        if typeID == None:
            continue
    
        if typeID not in _type_list:

            idx         = my_ff._impropers[typeID].index_torsion
            
            periodicity = my_ff._unique_torsion_list[idx].periodicity
            phase       = my_ff._unique_torsion_list[idx].phase
            kval        = my_ff._unique_torsion_list[idx].k

            if idx not in _unique_torsion_dict:
                _unique_torsion_dict[idx] = _dihedral_index
                prm.proper_periodicity_list += periodicity
                prm.proper_phase_list += phase
                prm.proper_k_list += kval
                _dihedral_index += len(phase)
                
            _type_dict[typeID] = MyTorsion (_unique_torsion_dict[idx], idx)
            
            prm.scee_list.append (Decimal(0.0))
            prm.scnb_list.append (Decimal(0.0))
            
        itype     = _type_dict[typeID].index + 1

        if iatElem._atomic_number == 1 or \
           jatElem._atomic_number == 1 or \
           latElem._atomic_number == 1:
#            print iatom, jatom, katom, latom
            prm.proper_wH_list.append(iatom*3)
            prm.proper_wH_list.append(jatom*3)
            prm.proper_wH_list.append(-katom*3)
            prm.proper_wH_list.append(-latom*3)
            prm.proper_wH_list.append(itype)
        else:
            prm.proper_nH_list.append(iatom*3)
            prm.proper_nH_list.append(jatom*3)
            prm.proper_nH_list.append(-katom*3)
            prm.proper_nH_list.append(-latom*3)
            prm.proper_nH_list.append(itype)


    
    for ii in range (len (_excludedAtomWith)):
        excluded_atoms_list = sorted(_excludedAtomWith[ii])
        numExAtom     = len(excluded_atoms_list)
        if numExAtom == 0:
            numExAtom = 1
            prm.num_excluded_atoms.append(numExAtom)
            prm.excluded_atoms_list.append(0)
        else:
            prm.num_excluded_atoms.append(numExAtom)
            for iatom in excluded_atoms_list:
                prm.excluded_atoms_list.append(iatom+1)
            
        
    prm.write (prmtop_fname) 


if __name__ == "__main__":

    pdb_fname = ''
    inpcrd_fname = ''
    prmtop_fname = ''
    argv = sys.argv[1:]
    
    opts, args = getopt.getopt (argv, "hi:o:c:", ["help=","ifile=", "ofile=","cfile="])

    if (len(opts) == 0):
        print ("pdb2amber.py -i <pdbfile.pdb> -o <prmtopfile.prmtop> -c <inpcrdfile.inpcrd>")
        sys.exit(2)
        
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ("pdb2amber.py -i <pdbfile.pdb> -o <prmtopfile.prmtop> -c <inpcrdfile.inpcrd>")
            sys.exit(1)
        elif opt in ("-i", "--ifile"):
            pdb_fname = arg;
        elif opt in ("-o", "--ofile"):
            prmtop_fname = arg;
        elif opt in ("-c", "--cfile"):
            inpcrd_fname = arg;
            
    ff_fnames = ['./data/protein.ff14SB.xml',\
            './data/lipid17.xml', \
            './data/fad.ff.xml', \
            './data/heme.ff.xml', \
            './data/uq2.ff.xml', \
            './data/ironsulfur.ff.xml', \
            './data/tip3pfb.xml'] 

    pdb2amber (pdb_fname, prmtop_fname, inpcrd_fname, ff_fnames)
