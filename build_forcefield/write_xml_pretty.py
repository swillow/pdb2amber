#!/usr/bin/python
# -*- coding: utf-8 -*-
try:
    from simtk.openmm.app.internal import amber_file_parser
    from simtk.openmm.app import *
    from simtk.openmm.vec3 import Vec3
    from simtk.openmm import *
    from simtk.unit import *
except:
    from openmm.app.internal import amber_file_parser
    from openmm.app import *
    from openmm.vec3 import Vec3
    from openmm import *
    from openmm.unit import *

import os
import sys
from collections import namedtuple
import numpy as np
import xml.etree.ElementTree as ET
import datetime
from decimal import Decimal


# import read_crdprm

def xml_indent(elem, level=0):
    i = '\n' + level * '  '
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + '  '
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            xml_indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def xml_ff_info(data):

    ff_info = ET.SubElement(data, 'Info')
    ff_info_date = ET.SubElement(ff_info, 'DateGenerated')
    ff_info_date.text = datetime.datetime.today().strftime('%Y-%m-%d')


def xml_ff_atom_types(data, prm, ff_prefix):

    element_list = [
        '',
        'H',
        'He',
        'Li',
        'Be',
        'B',
        'C',
        'N',
        'O',
        'F',
        'Ne',
    ]
    element_list += [
        'Na',
        'Mg',
        'Al',
        'Si',
        'P',
        'S',
        'Cl',
        'Ar',
    ]
    element_list += [
        'K',
        'Ca',
        'Sc',
        'Ti',
        'V',
        'Cr',
        'Mn',
        'Fe',
        'Co',
        'Ni',
        'Cu',
        'Zn',
    ]

    prm_atom_type = prm._raw_data['AMBER_ATOM_TYPE']
    prm_atom_znum = prm._raw_data['ATOMIC_NUMBER']
    prm_atom_mass = prm._raw_data['MASS']

    ff_atom_types = ET.SubElement(data, 'AtomTypes')
    atom_type_list = []
    for iat in range(len(prm_atom_type)):
        at_type = prm_atom_type[iat].upper()
        at_mass = '%8.4f' % float(prm_atom_mass[iat])
        at_name = ff_prefix + '-' + prm_atom_type[iat].upper()
        iznum = int(prm_atom_znum[iat])
        at_elem = element_list[iznum]

        if not at_type in atom_type_list:
            ff_atom_type = ET.SubElement(ff_atom_types, 'Type')
            ff_atom_type.set('name', at_name)
            ff_atom_type.set('class', at_type)
            ff_atom_type.set('element', at_elem)
            ff_atom_type.set('mass', at_mass.strip())
            atom_type_list.append(at_type)


def xml_ff_residues(data, prm, ff_prefix):

    res_name = prm._raw_data['RESIDUE_LABEL']
    prm_atom_name = prm._raw_data['ATOM_NAME']
    prm_atom_type = prm._raw_data['AMBER_ATOM_TYPE']
    prm_atom_chg = prm.getCharges()

    prm_bonds = prm._raw_data['BONDS_WITHOUT_HYDROGEN'] \
        + prm._raw_data['BONDS_INC_HYDROGEN']

    ff_residues = ET.SubElement(data, 'Residues')
    ff_residue = ET.SubElement(ff_residues, 'Residue')
    ff_residue.set('name', res_name[0].upper())

    natom = len(prm_atom_name)
    nheavy = 0
    tot_chg = 0.0
    for iat in range(natom):
        at_name = prm_atom_name[iat].upper()
        chg = float(prm_atom_chg[iat])
        tot_chg += chg
#        if at_name[0] != 'H':
#            nheavy += 1

    tot_chg_round = round(tot_chg)
    ave_chg = (tot_chg_round - tot_chg)/natom
    print(tot_chg, tot_chg_round, ave_chg)

    tot_chg = 0.0
    for iat in range(len(prm_atom_type)):
        at_name = prm_atom_name[iat].upper()
        chg = float(prm_atom_chg[iat]) + ave_chg

        at_chg = '%10.6f' % chg
        at_type = ff_prefix + '-' + prm_atom_type[iat].upper()

        tot_chg += float(at_chg)
        ff_res_atom = ET.SubElement(ff_residue, 'Atom')
        ff_res_atom.set('charge', at_chg)
        ff_res_atom.set('name', at_name)
        ff_res_atom.set('type', at_type)

    print('tot_chg', tot_chg)

    for ii in range(0, len(prm_bonds), 3):
        ibnd = int(prm_bonds[ii]) // 3
        jbnd = int(prm_bonds[ii + 1]) // 3

        iatnm = prm_atom_name[ibnd]
        jatnm = prm_atom_name[jbnd]

#    if iatnm == "O4P" or jatnm == "O4P": continue

        ff_res_bond = ET.SubElement(ff_residue, 'Bond')
        ff_res_bond.set('atomName1', iatnm)
        ff_res_bond.set('atomName2', jatnm)


def xml_ff_bond_force(data, prm, ff_prefix):

    bnd_k0 = prm._raw_data['BOND_FORCE_CONSTANT']
    bnd_r0 = prm._raw_data['BOND_EQUIL_VALUE']
    prm_atom_type = prm._raw_data['AMBER_ATOM_TYPE']

    forceConstConversionFactor = (kilocalorie_per_mole / (angstrom
                                                          * angstrom)).conversion_factor_to(kilojoule_per_mole
                                                                                            / (nanometer * nanometer))
    lengthConversionFactor = angstrom.conversion_factor_to(nanometer)

    ff_bnd_frc = ET.SubElement(data, 'HarmonicBondForce')
    raw_data = prm._raw_data['BONDS_WITHOUT_HYDROGEN'] \
        + prm._raw_data['BONDS_INC_HYDROGEN']
    typ_list = []
    for ii in range(0, len(raw_data), 3):
        ibnd = int(raw_data[ii]) // 3
        jbnd = int(raw_data[ii + 1]) // 3
        ityp = int(raw_data[ii + 2]) - 1
        r0 = '%10.6f' % (float(bnd_r0[ityp]) * lengthConversionFactor)

    # amber (k (x-x0)^2) ---> openmm (0.5 k' (x-x0)^2) : k' = 2k

        k0 = '%20.6f' % (float(bnd_k0[ityp])
                         * forceConstConversionFactor * 2.0)

        if not ityp in typ_list:
            typ_list.append(ityp)
            ff_bnd = ET.SubElement(ff_bnd_frc, 'Bond')
            itypeName = prm_atom_type[ibnd].upper()
            jtypeName = prm_atom_type[jbnd].upper()
            if itypeName > jtypeName:
                (itypeName, jtypeName) = (jtypeName, itypeName)

            ff_bnd.set('type1', ff_prefix + '-' + itypeName)
            ff_bnd.set('type2', ff_prefix + '-' + jtypeName)

            ff_bnd.set('length', r0.strip())
            ff_bnd.set('k', k0.strip())


def xml_ff_angle_force(data, prm, ff_prefix):

    forceConstConversionFactor = (kilocalorie_per_mole / (radian
                                                          * radian)).conversion_factor_to(kilojoule_per_mole
                                                                                          / (radian * radian))

    ang_k0 = prm._raw_data['ANGLE_FORCE_CONSTANT']
    ang_r0 = prm._raw_data['ANGLE_EQUIL_VALUE']
    prm_atom_type = prm._raw_data['AMBER_ATOM_TYPE']

    ff_ang_frc = ET.SubElement(data, 'HarmonicAngleForce')
    raw_data = prm._raw_data['ANGLES_WITHOUT_HYDROGEN'] \
        + prm._raw_data['ANGLES_INC_HYDROGEN']
    typ_list = []
    for ii in range(0, len(raw_data), 4):
        iang = int(raw_data[ii]) // 3
        jang = int(raw_data[ii + 1]) // 3
        kang = int(raw_data[ii + 2]) // 3
        ityp = int(raw_data[ii + 3]) - 1
        r0 = '%10.6f' % float(ang_r0[ityp])

    # amber (k (x-x0)^2) ---> openmm (0.5 k' (x-x0)^2) : k' = 2k

        k0 = '%20.6f' % (float(ang_k0[ityp])
                         * forceConstConversionFactor * 2.0)

        if not ityp in typ_list:
            typ_list.append(ityp)
            ff_ang = ET.SubElement(ff_ang_frc, 'Angle')
            itypeName = prm_atom_type[iang].upper()
            jtypeName = prm_atom_type[jang].upper()
            ktypeName = prm_atom_type[kang].upper()
            if itypeName > ktypeName:
                (itypeName, ktypeName) = (ktypeName, itypeName)

            ff_ang.set('type1', ff_prefix + '-' + itypeName)
            ff_ang.set('type2', ff_prefix + '-' + jtypeName)
            ff_ang.set('type3', ff_prefix + '-' + ktypeName)

            ff_ang.set('angle', r0.strip())
            ff_ang.set('k', k0.strip())


def xml_ff_torsion_force(data, prm, ff_prefix):
    forceConstConversionFactor = \
        kilocalorie_per_mole.conversion_factor_to(kilojoule_per_mole)
    forceConstant = prm._raw_data['DIHEDRAL_FORCE_CONSTANT']
    periodicity = prm._raw_data['DIHEDRAL_PERIODICITY']
    phase = prm._raw_data['DIHEDRAL_PHASE']
    prm_atom_type = prm._raw_data['AMBER_ATOM_TYPE']

    raw_data = prm._raw_data['DIHEDRALS_INC_HYDROGEN'] \
        + prm._raw_data['DIHEDRALS_WITHOUT_HYDROGEN']

    ff_tor_frc = ET.SubElement(data, 'PeriodicTorsionForce')

    proper_typ_list = []

    for ii in range(0, len(raw_data), 5):
        it = int(raw_data[ii]) // 3
        jt = int(raw_data[ii + 1]) // 3
        kt = abs(int(raw_data[ii + 2])) // 3
        lt = abs(int(raw_data[ii + 3])) // 3
        ityp = int(raw_data[ii + 4]) - 1

        fk0 = '%12.6f' % (float(forceConstant[ityp])
                          * forceConstConversionFactor)
        ph0 = '%20.16f' % Decimal(0.0)
        if Decimal(phase[ityp]) != Decimal(0.0):
            ph0 = '%20.16f' % Decimal(np.pi)
        pn0 = '%2d' % int(0.5 + float(periodicity[ityp]))

        if int(raw_data[ii + 3]) > 0:  # Proper
            typeID = 'PR_'
            typeID += prm_atom_type[it] + '_'
            typeID += prm_atom_type[jt] + '_'
            typeID += prm_atom_type[kt] + '_'
            typeID += prm_atom_type[lt]
            if not typeID in proper_typ_list:
                proper_typ_list.append(typeID)

                ff_proper = ET.SubElement(ff_tor_frc, 'Proper')
                ff_proper.set('k1', fk0.strip())
                ff_proper.set('periodicity1', pn0.strip())
                ff_proper.set('phase1', ph0.strip())
                ff_proper.set('type1', ff_prefix + '-'
                              + prm_atom_type[it].upper())
                ff_proper.set('type2', ff_prefix + '-'
                              + prm_atom_type[jt].upper())
                ff_proper.set('type3', ff_prefix + '-'
                              + prm_atom_type[kt].upper())
                ff_proper.set('type4', ff_prefix + '-'
                              + prm_atom_type[lt].upper())

    improp_typ_list = []

    for ii in range(0, len(raw_data), 5):
        it = int(raw_data[ii]) // 3
        jt = int(raw_data[ii + 1]) // 3
        kt = abs(int(raw_data[ii + 2])) // 3
        lt = abs(int(raw_data[ii + 3])) // 3
        ityp = int(raw_data[ii + 4]) - 1
        fk0 = '%12.6f' % (float(forceConstant[ityp])
                          * forceConstConversionFactor)
        ph0 = '%20.16f' % Decimal(0.0)
        if Decimal(phase[ityp]) != Decimal(0.0):
            ph0 = '%20.16f' % Decimal(np.pi)
        pn0 = '%2d' % int(0.5 + float(periodicity[ityp]))

        if int(raw_data[ii + 3]) < 0:  # improper
            typeID = 'IM_'
            typeID += prm_atom_type[it] + '_'
            typeID += prm_atom_type[jt] + '_'
            typeID += prm_atom_type[kt] + '_'
            typeID += prm_atom_type[lt]
            if not typeID in improp_typ_list:
                improp_typ_list.append(typeID)

                ff_improp = ET.SubElement(ff_tor_frc, 'Improper')
                ff_improp.set('k1', fk0.strip())
                ff_improp.set('periodicity1', pn0.strip())
                ff_improp.set('phase1', ph0.strip())
                ff_improp.set('type1', ff_prefix + '-'
                              + prm_atom_type[it].upper())
                ff_improp.set('type2', ff_prefix + '-'
                              + prm_atom_type[jt].upper())
                ff_improp.set('type3', ff_prefix + '-'
                              + prm_atom_type[kt].upper())
                ff_improp.set('type4', ff_prefix + '-'
                              + prm_atom_type[lt].upper())


def xml_ff_nonbond_force(data, prm, ff_prefix):

    rmin2sig = 1.0 / 2.0 ** (1.0 / 6.0)
    lengthConversionFactor = angstrom.conversion_factor_to(nanometer)
    energyConversionFactor = \
        kilocalorie_per_mole.conversion_factor_to(kilojoule_per_mole)

    prm_atom_type = prm._raw_data['AMBER_ATOM_TYPE']
    prm_atom_type_ind = prm._raw_data['ATOM_TYPE_INDEX']

    ff_nonbnd_frc = ET.SubElement(data, 'NonbondedForce')
    ff_nonbnd_frc.set('coulomb14scale', '0.8333333333333334')
    ff_nonbnd_frc.set('lj14scale', '0.5')
    ff_nonbnd_att = ET.SubElement(ff_nonbnd_frc,
                                  'UseAttributeFromResidue')
    ff_nonbnd_att.set('name', 'charge')

    numTypes = prm.getNumTypes()  # len(prm_atom_type)

    atom_type_list = []
    for ii in range(len(prm_atom_type)):
        at_type = ff_prefix + '-' + prm_atom_type[ii].upper()
        ityp = int(prm_atom_type_ind[ii]) - 1

        if not at_type in atom_type_list:
            atom_type_list.append(at_type)
            index = int(prm._raw_data['NONBONDED_PARM_INDEX'][numTypes
                                                              * ityp + ityp]) - 1
            if index < 0:
                continue

            acoef = float(prm._raw_data['LENNARD_JONES_ACOEF'][index])
            bcoef = float(prm._raw_data['LENNARD_JONES_BCOEF'][index])

            try:
                rmin_val = (2.0 * acoef / bcoef) ** (1.0 / 6.0)
                eps_val = 0.25 * bcoef * bcoef / acoef
            except ZeroDivisionError:
                rmin_val = 1.0
                eps_val = 0.0

            eps = '%10.7f' % Decimal(eps_val * energyConversionFactor)
            sig = '%20.16f' % Decimal(rmin_val * rmin2sig
                                      * lengthConversionFactor)

            ff_nonbnd_atm = ET.SubElement(ff_nonbnd_frc, 'Atom')
            ff_nonbnd_atm.set('epsilon', eps.strip())
            ff_nonbnd_atm.set('sigma', sig.strip())
            ff_nonbnd_atm.set('type', at_type)


if __name__ == '__main__':
    import json
    import getopt

    argv = sys.argv[1:]

    (opts, args) = getopt.getopt(argv, 'hi:', ['help=', 'input='])

    if len(opts) == 0:
        print('write_xml_pretty.py -i <input.json>')
        sys.exit(2)

    fname_json = 'input.json'
    for (opt, arg) in opts:
        if opt in ('-h', '--help'):
            print('write_xml_pretty.py -i <input.json>')
            sys.exit(1)
        elif opt in ('-i', '--input'):
            fname_json = arg

    with open(fname_json) as f:
        data = json.load(f)

        prmtop_fname = data['fname_prmtop']
        xml_fname = data['fname_xml']
        ff_prefix = data['ff_prefix']

        prm = amber_file_parser.PrmtopLoader(prmtop_fname)

        data = ET.Element('ForceField')

        xml_ff_info(data)

        xml_ff_atom_types(data, prm, ff_prefix)

        xml_ff_residues(data, prm, ff_prefix)

        xml_ff_bond_force(data, prm, ff_prefix)

        xml_ff_angle_force(data, prm, ff_prefix)

        xml_ff_torsion_force(data, prm, ff_prefix)

        xml_ff_nonbond_force(data, prm, ff_prefix)

        xml_indent(data)

        str_data = ET.tostring(data, method='xml').decode()
        xml_file = open(xml_fname, 'w')
        xml_file.write(str_data)
