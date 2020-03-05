import os
import itertools
import sys
import math
import numpy
import xml.etree.ElementTree as etree
from copy import copy
from collections import defaultdict
from datetime import date
from simtk.openmm import *
from simtk.openmm.app import *
import simtk.unit as units
from simtk.openmm.app import element as elem

"""
This is a modified 'simtk.openmm.app.forcefield.py' 
to generate the amber prmtop file for the non-standard Residue.
For example, in case that the amino acid is covalented to the co-factor or ligand, 
we cannot use the standard protein database (data/protein.ff14SB.xml).
"""

class AtomType (object):

    def __init__ (self, name, atomClass, mass, element):
        self.name      = name
        self.atomClass = atomClass
        self.mass      = mass
        self.element   = element

class AtomData (object):

    def __init__ (self, str_name, str_type, str_element, str_chg):

        self.at_name = str_name
        self.at_type = str_type
        self.at_element = str_element
        self.at_chg  = float(str_chg)
        self.bondedTo = []
        self.externalBonds = 0

class BondData (object):
    
    def __init__ (self, iatom, jatom):

        self.atom1 = iatom
        self.atom2 = jatom
        self.length = 0.0

        if iatom > jatom:
            self.atom2 = iatom
            self.atom1 = jatom


class ResidueData (object):
    """
     AMBER PRMTOP did not support virtual sites with massless.
    """
    def __init__ (self, res_name):

        self.name  = res_name
        self.atoms = []  # AtomData
        self.atomIndices = {}
        self.bonds = []
        self.externalBonds = []
        

    def addBond (self, atom1, atom2):
        if atom1 < atom2:
            self.bonds.append ( (atom1, atom2) )
        else:
            self.bonds.append ( (atom2, atom1) )
        self.atoms[atom1].bondedTo.append (atom2)
        self.atoms[atom2].bondedTo.append (atom1)


    def addBondByName (self, atom1_name, atom2_name):
        atom1 = self.atomIndices[atom1_name]
        atom2 = self.atomIndices[atom2_name]
        self.addBond (atom1, atom2)


    def addExternalBond (self, atom1):
        self.externalBonds.append (atom1)
        self.atoms[atom1].externalBonds += 1

    def addExternalBondByName (self, atom1_name):
        atom1 = self.atomIndices[atom1_name]
        self.addExternalBond (atom1)


class HarmonicBondData(object):
    """ a HarmonicBondForce. """
    def __init__ (self, type1, type2, length, k):
        self.types1 = type1
        self.types2 = type2
        if type1 > type2:
            self.types1 = type2
            self.types2 = type1
        self.length = length
        self.k      = k


class HarmonicAngleData(object):
    """ a HarmonicAngleForce. """
    def __init__ (self, type1, type2, type3, angle, k):
        self.types1 = type1
        self.types2 = type2
        self.types3 = type3
        if type1 > type3:
            self.types1 = type3
            self.types3 = type1
        self.angle  = angle
        self.k      = k


class PeriodicTorsionData (object):

    def __init__ (self):
        self.periodicity = []
        self.phase       = []
        self.k           = []


class PeriodicTorsion(object):

    def __init__(self, types):

        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.index_torsion = -1
        self.ordering    = 'default'



class MyForceFields (object):


    def __init__ (self, file_names): #pdb_file_name):
        # read filename

        self._atomTypes = {}
        self._residues  = {}
        self._atomClasses = {'':set()}
        self._harmonicBonds = []
        self._harmonicAngles = []
        self._propers   = {}
        self._impropers = {}
        self._nonbonds  = {}
        self._nonbonds_coulomb14scale = 0.8333333333333334
        self._nonbonds_lj14scale      = 0.5
        self._bondsForAtomType = defaultdict(set)
        self._unique_torsion_list = []
        self._unique_lj_list = []

        trees = []

        for xml_file_name in file_names:
            print ('file_names', xml_file_name)
            tree = etree.parse (xml_file_name)
            trees.append(tree)

        # Load the atom types
        for tree in trees:
            element = tree.getroot().find('AtomTypes')
            if element is not None:
                for at_type in element.findall('Type'):
                    typeName    = at_type.attrib['name']
                    atomClass   = at_type.attrib['class']
                    atomMass    = float(at_type.attrib['mass'])
                    atomElement = None
                    if 'element' in at_type.attrib:
                        atomElement = at_type.attrib['element']
                    self._atomTypes[typeName] = AtomType (typeName, atomClass,
                                                          atomMass, atomElement)

                    if atomClass not in self._atomClasses:
                        self._atomClasses[atomClass] = set()

                    self._atomClasses[atomClass].add(typeName)
                    self._atomClasses[''].add(typeName)



        # Load the Residue templates
        for tree in trees:
            element = tree.getroot().find('Residues')
            if element is not None:
                for residue in element.findall('Residue'):
                    resName     = residue.attrib['name']
                    resData     = ResidueData (resName)

                    #sum_charge = 0.0

                    for ia, atom in enumerate (residue.findall('Atom')):
                        atomName = atom.attrib['name']
                        typeName = atom.attrib['type']
                        charge   = atom.attrib['charge']

                     #   sum_charge += float(charge)

                        resData.atomIndices[atomName] = ia

                        atomData = AtomData(atomName,
                                            typeName,
                                            self._atomTypes[typeName].element,
                                            charge)

                        resData.atoms.append (atomData)

                    #if abs(sum_charge) > 0.001:
                    #    print (resName, sum_charge)

                    for bond in residue.findall('Bond'):
                        resData.addBondByName(bond.attrib['atomName1'],
                                              bond.attrib['atomName2'])

                    for bond in residue.findall('ExternalBond'):
                        resData.addExternalBondByName(bond.attrib['atomName'])


                    # Register Template
                    self._residues[resData.name] = resData

        # Load the HarmonicBondForce
        # OpenMM = 1/2 k (x - x0)**2 : E(kJ/mol/nm**2)
        # Amber  = k' (x-x0)**2      : E(kcal/mol/A**2)
        # k' = 1/2 k
        ene_conv = (units.kilojoule_per_mole/(units.nanometer*units.nanometer)
                             ).conversion_factor_to( units.kilocalorie_per_mole/(units.angstrom*units.angstrom))
        len_conv = units.nanometer.conversion_factor_to(units.angstrom)

        for tree in trees:
            element = tree.getroot().find('HarmonicBondForce')
            if element is not None:
                for bond in element.findall('Bond'):
                    type1 = bond.attrib['type1']
                    type2 = bond.attrib['type2']
                    bond_length = float(bond.attrib['length'])*len_conv
                    bond_k      = 0.5*float(bond.attrib['k'])*ene_conv
                    harmBond    = HarmonicBondData (type1, type2, bond_length, bond_k)
                    self._harmonicBonds.append ( harmBond )


        ene_conv = units.kilojoule_per_mole.conversion_factor_to (units.kilocalorie_per_mole)
        for tree in trees:
            element = tree.getroot().find('HarmonicAngleForce')
            if element is not None:
                for angle in element.findall('Angle'):
                    type1 = angle.attrib['type1']
                    type2 = angle.attrib['type2']
                    type3 = angle.attrib['type3']
                    ang_length = float(angle.attrib['angle'])
                    ang_k      = 0.5*float(angle.attrib['k'])*ene_conv
                    harmAngle  = HarmonicAngleData (type1, type2, type3,
                                                    ang_length, ang_k)
                    self._harmonicAngles.append ( harmAngle )


        ene_conv = units.kilojoule_per_mole.conversion_factor_to (units.kilocalorie_per_mole)
        for tree in trees:
            element = tree.getroot().find('PeriodicTorsionForce')
            if element is not None:

                ordering = 'default'
                if 'ordering' in element.attrib:
                    ordering = element.attrib['ordering']

                for proper in element.findall('Proper'):
                    types = []
                    for i in range(4):
                        suffix = str(i+1)
                        typeAttrib = 'type'+suffix
                        typeName   = proper.attrib[typeAttrib]
                        if typeName == '':
                            types.append('X') #self._atomClasses[''])
                        elif typeName not in self._atomTypes:
                            types.append(None)
                            print ('Proper ', typeName)
                        else:
                            types.append(typeName)

                    torsion = PeriodicTorsion(types)
                    index  = 1
                    torsionData = PeriodicTorsionData()

                    while 'phase%d'%index in proper.attrib:
                        periodicity = int(proper.attrib['periodicity%d'%index])
                        phase = float(proper.attrib['phase%d'%index])
                        k     = float(proper.attrib['k%d'%index])*ene_conv
                        torsionData.periodicity.append(periodicity)
                        torsionData.phase.append(phase)
                        torsionData.k.append(k)
                        index += 1

                    ladd = 1
                    for uniqueTorsion in self._unique_torsion_list:
                        if torsionData.periodicity == uniqueTorsion.periodicity and \
                           torsionData.phase       == uniqueTorsion.phase and \
                           torsionData.k           == uniqueTorsion.k:
                            ladd = 0
                            break

                    if ladd == 1:
                        self._unique_torsion_list.append (torsionData)

                    for ii in range (len(self._unique_torsion_list)):
                        uniqueTorsion = self._unique_torsion_list[ii]
                        if torsionData.periodicity == uniqueTorsion.periodicity and \
                           torsionData.phase       == uniqueTorsion.phase and \
                           torsionData.k           == uniqueTorsion.k:
                            torsion.index_torsion = ii

                    typeID  = 'PR_'
                    if types[0] == types[3]:
                        if types[1] < types[2]:
                            typeID += types[0]+'_'
                            typeID += types[1]+'_'
                            typeID += types[2]+'_'
                            typeID += types[3]
                        else:
                            typeID += types[3]+'_'
                            typeID += types[2]+'_'
                            typeID += types[1]+'_'
                            typeID += types[0]

                    elif types[0] < types[3]:
                        typeID += types[0]+'_'
                        typeID += types[1]+'_'
                        typeID += types[2]+'_'
                        typeID += types[3]
                    else:
                        typeID += types[3]+'_'
                        typeID += types[2]+'_'
                        typeID += types[1]+'_'
                        typeID += types[0]

                    self._propers[typeID] = torsion

                for improper in element.findall('Improper'):
                    types = []
                    improper_openmm2amber = [2,3,1,4]
                    for i in improper_openmm2amber:
                        suffix = str(i)
                        typeAttrib = 'type'+suffix
                        typeName   = improper.attrib[typeAttrib]
                        if typeName == '':
                            types.append('X')#self._atomClasses[''])
                        elif typeName not in self._atomTypes:
                            types.append(None)
                            print ('ImProper ', typeName)
                        else:
                            types.append(typeName)

                    torsion = PeriodicTorsion(types)
                    torsionData = PeriodicTorsionData()

                    torsion.ordering = ordering
                    periodicity = int(improper.attrib['periodicity1'])
                    phase       = float(improper.attrib['phase1'])
                    k           = float(improper.attrib['k1'])*ene_conv
                    torsionData.periodicity.append(periodicity)
                    torsionData.phase.append(phase)
                    torsionData.k.append(k)

                    ladd = 1
                    for uniqueTorsion in self._unique_torsion_list:
                        if torsionData.periodicity == uniqueTorsion.periodicity and \
                           torsionData.phase       == uniqueTorsion.phase and \
                           torsionData.k           == uniqueTorsion.k:
                            ladd = 0
                            break

                    if ladd == 1:
                        self._unique_torsion_list.append (torsionData)

                    for ii in range (len(self._unique_torsion_list)):
                        uniqueTorsion = self._unique_torsion_list[ii]
                        if torsionData.periodicity == uniqueTorsion.periodicity and \
                           torsionData.phase       == uniqueTorsion.phase and \
                           torsionData.k           == uniqueTorsion.k:
                            torsion.index_torsion = ii


                    if types[0] > types[1]:
                        (types[0], types[1]) = (types[1], types[0])

                    typeID  = 'IM_'
                    typeID += types[0]+'_'
                    typeID += types[1]+'_'
                    typeID += types[2]+'_'
                    typeID += types[3]

                    self._impropers[typeID] = torsion

        # Eunit: kJ/mol --> kcal/mol
        # Lunit: nm --> A
        ene_conv = units.kilojoules_per_mole.conversion_factor_to(units.kilocalories_per_mole) 
        len_conv = units.nanometers.conversion_factor_to(units.angstrom) 
        # Load NonbondedForce
        for tree in trees:
            element = tree.getroot().find('NonbondedForce')
            if element is not None:

                if 'coulomb14scale' in element.attrib:
                    self._nonbonds_coulomb14scale = float(element.attrib['coulomb14scale'])
                if 'lj14scale' in element.attrib:
                    self._nonbonds_lj14scale      = float(element.attrib['lj14scale'])
                for atom in element.findall('Atom'):
                    typeName = atom.attrib['type']

                    epsilon  = float(atom.attrib['epsilon'])*ene_conv
                    sigma    = float(atom.attrib['sigma'])*len_conv

                    self._nonbonds[typeName] = [epsilon, sigma]
