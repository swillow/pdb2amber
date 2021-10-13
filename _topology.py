"""
topology.py: Used for storing topological information about a system.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2018 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import
__author__ = "Peter Eastman"  # modified by Soohaeng Yoo Willow
__version__ = "1.0"

from collections import namedtuple
import os
import xml.etree.ElementTree as etree
try:
    from simtk.openmm.vec3 import Vec3
    from simtk.openmm.app.internal.singleton import Singleton
    from simtk.unit import nanometers, sqrt, is_quantity
except:
    from openmm.vec3 import Vec3
    from openmm.app.internal.singleton import Singleton
    from openmm.unit import nanometers, sqrt, is_quantity
    
from copy import deepcopy
import sys
# Enumerated values for bond type

_standardResidues = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX', 'CYG', 'CYF', 'GLH', 'GLN', 'GLU',
                     'GLY', 'HIS', 'HID', 'HIE', 'HIP', 'HYP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE',
                     'PRO', 'SER', 'THO', 'THR', 'TRP', 'TYR', 'VAL', 'HI5']

_atomNameReplacements = {'HN': 'H', 'H1': 'H', '1H': 'H', 'HN1': 'H', 'HT1': 'H',
                         '2H': 'H2', 'HN2': 'H2', 'HT2': 'H2',
                         '3H': 'H3', 'HN3': 'H3', 'HT3': 'H3',
                         'O1': 'O', 'OT1': 'O', 'OCT1': 'O', 'OC1': 'O',
                         'O2': 'OXT', 'OT2': 'OXT', 'OCT2': 'OXT', 'OC2': 'OXT', 'OT': 'OXT'}


class Single(Singleton):
    def __repr__(self):
        return 'Single'


Single = Single()


class Double(Singleton):
    def __repr__(self):
        return 'Double'


Double = Double()


class Triple(Singleton):
    def __repr__(self):
        return 'Triple'


Triple = Triple()


class Aromatic(Singleton):
    def __repr__(self):
        return 'Aromatic'


Aromatic = Aromatic()


class Amide(Singleton):
    def __repr__(self):
        return 'Amide'


Amide = Amide()


class MyTopology(object):
    """Topology stores the topological information about a system.

    The structure of a Topology object is similar to that of a PDB file.  It consists of a set of Chains
    (often but not always corresponding to polymer chains).  Each Chain contains a set of Residues,
    and each Residue contains a set of Atoms.  In addition, the Topology stores a list of which atom
    pairs are bonded to each other, and the dimensions of the crystallographic unit cell.

    Atom and residue names should follow the PDB 3.0 nomenclature for all molecules for which one exists.
    """

    _standardBonds = {}
    _hasLoadedStandardBonds = False

    def __init__(self):
        """Create a new Topology object"""
        self._chains = []
        self._numResidues = 0
        self._numAtoms = 0
        self._bonds = []
        self._periodicBoxVectors = None

    def __repr__(self):
        nchains = len(self._chains)
        nres = self._numResidues
        natom = self._numAtoms
        nbond = len(self._bonds)
        return '<%s; %d chains, %d residues, %d atoms, %d bonds>' % (
            type(self).__name__, nchains, nres, natom, nbond)

    def getNumAtoms(self):
        """Return the number of atoms in the Topology.
        """
        return self._numAtoms

    def getNumResidues(self):
        """Return the number of residues in the Topology.
        """
        return self._numResidues

    def getNumChains(self):
        """Return the number of chains in the Topology.
        """
        return len(self._chains)

    def getNumBonds(self):
        """Return the number of bonds in the Topology.
        """
        return len(self._bonds)

    def addChain(self, id=None):
        """Create a new Chain and add it to the Topology.

        Parameters
        ----------
        id : string=None
            An optional identifier for the chain.  If this is omitted, an id is
            generated based on the chain index.

        Returns
        -------
        Chain
             the newly created Chain
        """
        if id is None:
            id = str(len(self._chains)+1)
        chain = Chain(len(self._chains), self, id)
        self._chains.append(chain)
        return chain

    def addResidue(self, name, chain, id=None, insertionCode=''):
        """Create a new Residue and add it to the Topology.

        Parameters
        ----------
        name : string
            The name of the residue to add
        chain : Chain
            The Chain to add it to
        id : string=None
            An optional identifier for the residue.  If this is omitted, an id
            is generated based on the residue index.
        insertionCode: string=''
            An optional insertion code for the residue.

        Returns
        -------
        Residue
             the newly created Residue
        """
        if len(chain._residues) > 0 and self._numResidues != chain._residues[-1].index+1:
            raise ValueError('All residues within a chain must be contiguous')
        if id is None:
            id = str(self._numResidues+1)
        residue = Residue(name, self._numResidues, chain, id, insertionCode)
        self._numResidues += 1
        chain._residues.append(residue)
        return residue

    def addAtom(self, name, element, residue, id=None):
        """Create a new Atom and add it to the Topology.

        Parameters
        ----------
        name : string
            The name of the atom to add
        element : Element
            The element of the atom to add
        residue : Residue
            The Residue to add it to
        id : string=None
            An optional identifier for the atom.  If this is omitted, an id is
            generated based on the atom index.

        Returns
        -------
        Atom
             the newly created Atom
        """
        if len(residue._atoms) > 0 and self._numAtoms != residue._atoms[-1].index+1:
            raise ValueError('All atoms within a residue must be contiguous')
        if id is None:
            id = str(self._numAtoms+1)
        atom = Atom(name, element, self._numAtoms, residue, id)
        self._numAtoms += 1
        residue._atoms.append(atom)
        return atom

    def addBond(self, atom1, atom2, type=None, order=None):
        """Create a new bond and add it to the Topology.

        Parameters
        ----------
        atom1 : Atom
            The first Atom connected by the bond
        atom2 : Atom
            The second Atom connected by the bond
        type : object=None
            The type of bond to add.  Allowed values are None, Single, Double, Triple,
            Aromatic, or Amide.
        order : int=None
            The bond order, or None if it is not specified
        """
        self._bonds.append(Bond(atom1, atom2, type, order))

    def chains(self):
        """Iterate over all Chains in the Topology."""
        return iter(self._chains)

    def residues(self):
        """Iterate over all Residues in the Topology."""
        for chain in self._chains:
            for residue in chain._residues:
                yield residue

    def atoms(self):
        """Iterate over all Atoms in the Topology."""
        for chain in self._chains:
            for residue in chain._residues:
                for atom in residue._atoms:
                    yield atom

    def bonds(self):
        """Iterate over all bonds (each represented as a tuple of two Atoms) in the Topology."""
        return iter(self._bonds)

    def getPeriodicBoxVectors(self):
        """Get the vectors defining the periodic box.

        The return value may be None if this Topology does not represent a periodic structure."""
        return self._periodicBoxVectors

    def setPeriodicBoxVectors(self, vectors):
        """Set the vectors defining the periodic box."""
        if vectors is not None:
            if not is_quantity(vectors[0][0]):
                vectors = vectors*nanometers
            if vectors[0][1] != 0*nanometers or vectors[0][2] != 0*nanometers:
                raise ValueError(
                    "First periodic box vector must be parallel to x.")
            if vectors[1][2] != 0*nanometers:
                raise ValueError(
                    "Second periodic box vector must be in the x-y plane.")
            if vectors[0][0] <= 0*nanometers or vectors[1][1] <= 0*nanometers or vectors[2][2] <= 0*nanometers or vectors[0][0] < 2*abs(vectors[1][0]) or vectors[0][0] < 2*abs(vectors[2][0]) or vectors[1][1] < 2*abs(vectors[2][1]):
                raise ValueError(
                    "Periodic box vectors must be in reduced form.")
        self._periodicBoxVectors = deepcopy(vectors)

    def getUnitCellDimensions(self):
        """Get the dimensions of the crystallographic unit cell.

        The return value may be None if this Topology does not represent a periodic structure.
        """
        if self._periodicBoxVectors is None:
            return None
        xsize = self._periodicBoxVectors[0][0].value_in_unit(nanometers)
        ysize = self._periodicBoxVectors[1][1].value_in_unit(nanometers)
        zsize = self._periodicBoxVectors[2][2].value_in_unit(nanometers)
        return Vec3(xsize, ysize, zsize)*nanometers

    def setUnitCellDimensions(self, dimensions):
        """Set the dimensions of the crystallographic unit cell.

        This method is an alternative to setPeriodicBoxVectors() for the case of a rectangular box.  It sets
        the box vectors to be orthogonal to each other and to have the specified lengths."""
        if dimensions is None:
            self._periodicBoxVectors = None
        else:
            if is_quantity(dimensions):
                dimensions = dimensions.value_in_unit(nanometers)
            self._periodicBoxVectors = (Vec3(dimensions[0], 0, 0), Vec3(
                0, dimensions[1], 0), Vec3(0, 0, dimensions[2]))*nanometers

    @staticmethod
    def loadBondDefinitions(ff_file):
        """Load an XML file containing definitions of bonds that should be used by createStandardBonds().

        The built in residues.xml file containing definitions for standard amino acids and nucleotides is loaded automatically.
        This method can be used to load additional definitions for other residue types.  They will then be used in subsequent
        calls to createStandardBonds().  This is a static method, so it affects subsequent calls on all Topology objects.
        Also note that PDBFile calls createStandardBonds() automatically when a file is loaded, so the newly loaded definitions
        will be used for any PDB file loaded after this is called.
        """

        tree = etree.parse(ff_file)
        element = tree.getroot().find('Residues')
        for residue in element.findall('Residue'):
            bonds = []
            resname = residue.attrib['name']
            l_pdb = False
            if len(resname) >= 3:
                if resname[-3:] in _standardResidues:
                    l_pdb = True
            MyTopology._standardBonds[resname] = bonds
            for bond in residue.findall('Bond'):
                atnm1 = bond.attrib['atomName1']
                atnm2 = bond.attrib['atomName2']
                if l_pdb:
                    if atnm1 in _atomNameReplacements:
                        atnm1 = _atomNameReplacements[atnm1]
                    if atnm2 in _atomNameReplacements:
                        atnm2 = _atomNameReplacements[atnm2]

                #bonds.append((bond.attrib['from'], bond.attrib['to']))
                bonds.append((atnm1, atnm2))

            exBonds = []
            for bond in residue.findall('ExternalBond'):
                exBonds.append(bond.attrib['atomName'])
            if l_pdb:
                # Protein Backbone
                bonds.append(('C', '+N'))
            elif len(exBonds) == 2:
                # Nucleic Acids (DNA and RNA) Backbone
                if exBonds[0] in ['O3', 'P'] and exBonds[1] in ['O3', 'P']:
                    bonds.append(('-O3', 'P'))

    def createStandardBonds(self, positions, ff_fileNames):
        """Create bonds based on the atom and residue names for all standard residue types.

        Definitions for standard amino acids and nucleotides are built in.  You can call loadBondDefinitions() to load
        additional definitions for other residue types.
        """
        if not MyTopology._hasLoadedStandardBonds:
            # Load the standard bond definitions.

            for ff_fname in ff_fileNames:
                MyTopology.loadBondDefinitions(ff_fname)
            MyTopology._hasLoadedStandardBonds = True

        for chain in self._chains:
            # First build a map of atom names to atoms.
            nres = chain._residues[0]
            cres = chain._residues[-1]

            atomMaps = []

            for residue in chain._residues:
                atomMap = {}
                atomMaps.append(atomMap)

                for atom in residue._atoms:
                    # if atom.name in _atomNameReplacements:
                    #    atom.name = _atomNameReplacements[atom.name]
                    atomMap[atom.name] = atom

            # Loop over residues and construct bonds.

            for i in range(len(chain._residues)):
                res = chain._residues[i]
                name = res.name
                if name in _standardResidues:
                    if res == nres:
                        name = "N" + res.name
                    elif res == cres:
                        name = "C" + res.name

                    if name in ['HIE', 'HIS']:
                        """
                        HIE: no HD1, yes HE2
                        HID: yes HD1, no HE2
                        HIP: yes HD1, yes HE2
                        """
                        atomNames = []
                        name = 'HIE'
                        for atom in res._atoms:
                            atomNames.append(atom.name)

                        if 'HD1' in atomNames:
                            if 'HE2' in atomNames:
                                name = 'HIP'
                            else:
                                name = 'HID'
                        res.name = name

                if name in MyTopology._standardBonds:
                    for bond in MyTopology._standardBonds[name]:

                        if bond[0].startswith('-') and i > 0:
                            fromResidue = i-1
                            fromAtom = bond[0][1:]
                        elif bond[0].startswith('+') and i < len(chain._residues)-1:
                            fromResidue = i+1
                            fromAtom = bond[0][1:]
                        else:
                            fromResidue = i
                            fromAtom = bond[0]
                        if bond[1].startswith('-') and i > 0:
                            toResidue = i-1
                            toAtom = bond[1][1:]
                        elif bond[1].startswith('+') and i < len(chain._residues)-1:
                            toResidue = i+1
                            toAtom = bond[1][1:]
                        else:
                            toResidue = i
                            toAtom = bond[1]

                        l_add_Bond = False
                        if fromAtom in atomMaps[fromResidue] and toAtom in atomMaps[toResidue]:
                            l_add_Bond = True

                        if fromResidue != toResidue and l_add_Bond:
                            atom1 = atomMaps[fromResidue][fromAtom]
                            atom2 = atomMaps[toResidue][toAtom]
                            pos1 = positions[atom1.index]
                            pos2 = positions[atom2.index]
                            delta = [x - y for (x, y) in zip(pos1, pos2)]
                            distance = sqrt(
                                delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2])
                            if distance > 0.20*nanometers:
                                l_add_Bond = False

                        if l_add_Bond:
                            self.addBond(
                                atomMaps[fromResidue][fromAtom], atomMaps[toResidue][toAtom])

        # sys.exit()

    def createDisulfideBonds(self, positions):
        """Identify disulfide bonds based on proximity and add them to the
        Topology.

        Parameters
        ----------
        positions : list
            The list of atomic positions based on which to identify bonded atoms
        """
        def isCyx(res):
            names = [atom.name for atom in res._atoms]
            return 'SG' in names and 'HG' not in names

        cyx_sg_list = []
        for res in self.residues():
            is_cyx = False
            if res.name == 'CYX':
                is_cyx = True
            elif res.name == 'CYS' and isCyx(res):
                is_cyx = True

            if is_cyx:
                atomNames = [atom.name for atom in res._atoms]
                sg = res._atoms[atomNames.index('SG')]
                cyx_sg_list.append(sg)

        for i in range(len(cyx_sg_list)):
            sg1 = cyx_sg_list[i]
            pos1 = positions[sg1.index]
            for j in range(i):
                sg2 = cyx_sg_list[j]
                pos2 = positions[sg2.index]
                delta = [x-y for (x, y) in zip(pos1, pos2)]
                distance = sqrt(delta[0]*delta[0] + delta[1]
                                * delta[1] + delta[2]*delta[2])
                if distance < 0.23*nanometers:
                    self.addBond(sg1, sg2)

    def createIronSulfurBonds(self, positions):
        import sys

        print('createIronSulfureBonds ')

        def isCyx(res):
            names = [atom.name for atom in res._atoms]
            return 'SG' in names and 'HG' not in names

        cyx_sg_list = []
        for res in self.residues():
            is_cyx = False
            if res.name in ['CYX', 'CYF', 'CYG', 'CY4', 'CY3']:
                is_cyx = True
            elif res.name == 'CYS' and isCyx(res):
                is_cyx = True

            if is_cyx:
                atomNames = [atom.name for atom in res._atoms]
                sg = res._atoms[atomNames.index('SG')]
                cyx_sg_list.append(sg)

        fe_list = []
        for res in self.residues():
            if res.name in ['FES', 'FET']:
                atomNames = [atom.name for atom in res._atoms]
                fe = res._atoms[atomNames.index('FE1')]
                fe_list.append(fe)
                fe = res._atoms[atomNames.index('FE2')]
                fe_list.append(fe)

            if res.name == 'SF4':
                atomNames = [atom.name for atom in res._atoms]
                fe = res._atoms[atomNames.index('FE1')]
                fe_list.append(fe)
                fe = res._atoms[atomNames.index('FE2')]
                fe_list.append(fe)
                fe = res._atoms[atomNames.index('FE3')]
                fe_list.append(fe)
                fe = res._atoms[atomNames.index('FE4')]
                fe_list.append(fe)

            if res.name == 'F3S':
                atomNames = [atom.name for atom in res._atoms]
                fe = res._atoms[atomNames.index('FE1')]
                fe_list.append(fe)
                fe = res._atoms[atomNames.index('FE3')]
                fe_list.append(fe)
                fe = res._atoms[atomNames.index('FE4')]
                fe_list.append(fe)

            if res.name == 'FE':
                atomNames = [atom.name for atom in res._atoms]
                fe = res._atoms[atomNames.index('FE')]
                fe_list.append(fe)

        for fe in fe_list:
            pos1 = positions[fe.index]
            print('---FE---  ', fe.residue.name, fe.index, pos1)

            for sg in cyx_sg_list:
                pos2 = positions[sg.index]

                delta = [x1 - x2 for (x1, x2) in zip(pos1, pos2)]
                distance = sqrt(delta[0]*delta[0] + delta[1]
                                * delta[1] + delta[2]*delta[2])

                if distance < 0.3*nanometers:
                    print('addBond ', fe.index, sg.index, distance)
                    self.addBond(fe, sg)

    # resName1, atomName1, resName2, atomName2):
    def createUserDefined(self, positions, link_residues):
        import sys

        res1_atoms = {}
        res2_atoms = {}

        for res1, res2 in link_residues:
            atomName, resName = res1.split()
            if resName not in res1_atoms:
                res1_atoms[resName] = []
            if atomName not in res1_atoms[resName]:
                res1_atoms[resName].append(atomName)

            atomName, resName = res2.split()
            if resName not in res2_atoms:
                res2_atoms[resName] = []
            if atomName not in res2_atoms[resName]:
                res2_atoms[resName].append(atomName)

        atom1_list = []
        atom2_list = []
        for res in self.residues():
            if res.name in res1_atoms:
                for atom1 in res._atoms:
                    if atom1.name in res1_atoms[res.name]:
                        atom1_list.append(atom1)

            if res.name in res2_atoms:
                for atom2 in res._atoms:
                    if atom2.name in res2_atoms[res.name]:
                        atom2_list.append(atom2)

        for atom1 in atom1_list:
            pos1 = positions[atom1.index]

            for atom2 in atom2_list:
                pos2 = positions[atom2.index]

                delta = [x1 - x2 for (x1, x2) in zip(pos1, pos2)]
                distance = sqrt(delta[0]*delta[0] + delta[1]
                                * delta[1] + delta[2]*delta[2])

                if distance < 0.3*nanometers:
                    print('addBond ', atom1.index, atom2.index, distance)
                    self.addBond(atom1, atom2)


class Chain(object):
    """A Chain object represents a chain within a Topology."""

    def __init__(self, index, topology, id):
        """Construct a new Chain.  You should call addChain() on the Topology instead of calling this directly."""
        # The index of the Chain within its Topology
        self.index = index
        # The Topology this Chain belongs to
        self.topology = topology
        # A user defined identifier for this Chain
        self.id = id
        self._residues = []

    def residues(self):
        """Iterate over all Residues in the Chain."""
        return iter(self._residues)

    def atoms(self):
        """Iterate over all Atoms in the Chain."""
        for residue in self._residues:
            for atom in residue._atoms:
                yield atom

    def __len__(self):
        return len(self._residues)

    def __repr__(self):
        return "<Chain %d>" % self.index


class Residue(object):
    """A Residue object represents a residue within a Topology."""

    def __init__(self, name, index, chain, id, insertionCode):
        """Construct a new Residue.  You should call addResidue() on the Topology instead of calling this directly."""
        # The name of the Residue
        self.name = name
        # The index of the Residue within its Topology
        self.index = index
        # The Chain this Residue belongs to
        self.chain = chain
        # A user defined identifier for this Residue
        self.id = id
        # A user defined insertion code for this Residue
        self.insertionCode = insertionCode
        self._atoms = []

    def atoms(self):
        """Iterate over all Atoms in the Residue."""
        return iter(self._atoms)

    def bonds(self):
        """Iterate over all Bonds involving any atom in this residue."""
        return (bond for bond in self.chain.topology.bonds() if ((bond[0] in self._atoms) or (bond[1] in self._atoms)))

    def internal_bonds(self):
        """Iterate over all internal Bonds."""
        return (bond for bond in self.chain.topology.bonds() if ((bond[0] in self._atoms) and (bond[1] in self._atoms)))

    def external_bonds(self):
        """Iterate over all Bonds to external atoms."""
        return (bond for bond in self.chain.topology.bonds() if ((bond[0] in self._atoms) != (bond[1] in self._atoms)))

    def __len__(self):
        return len(self._atoms)

    def __repr__(self):
        return "<Residue %d (%s) of chain %d>" % (self.index, self.name, self.chain.index)


class Atom(object):
    """An Atom object represents an atom within a Topology."""

    def __init__(self, name, element, index, residue, id):
        """Construct a new Atom.  You should call addAtom() on the Topology instead of calling this directly."""
        # The name of the Atom
        self.name = name
        # That Atom's element
        self.element = element
        # The index of the Atom within its Topology
        self.index = index
        # The Residue this Atom belongs to
        self.residue = residue
        # A user defined identifier for this Atom
        self.id = id

    def __repr__(self):
        return "<Atom %d (%s) of chain %d residue %d (%s)>" % (self.index, self.name, self.residue.chain.index, self.residue.index, self.residue.name)


class Bond(namedtuple('Bond', ['atom1', 'atom2'])):
    """A Bond object represents a bond between two Atoms within a Topology.

    This class extends tuple, and may be interpreted as a 2 element tuple of Atom objects.
    It also has fields that can optionally be used to describe the bond order and type of bond."""

    def __new__(cls, atom1, atom2, type=None, order=None):
        """Create a new Bond.  You should call addBond() on the Topology instead of calling this directly."""
        bond = super(Bond, cls).__new__(cls, atom1, atom2)
        bond.type = type
        bond.order = order
        return bond

    def __getnewargs__(self):
        "Support for pickle protocol 2: http://docs.python.org/2/library/pickle.html#pickling-and-unpickling-normal-class-instances"
        return self[0], self[1], self.type, self.order

    def __getstate__(self):
        """
        Additional support for pickle since parent class implements its own __getstate__
        so pickle does not store or restore the type and order, python 2 problem only
        https://www.python.org/dev/peps/pep-0307/#case-3-pickling-new-style-class-instances-using-protocol-2
        """
        return self.__dict__

    def __deepcopy__(self, memo):
        return Bond(self[0], self[1], self.type, self.order)

    def __repr__(self):
        s = "Bond(%s, %s" % (self[0], self[1])
        if self.type is not None:
            s = "%s, type=%s" % (s, self.type)
        if self.order is not None:
            s = "%s, order=%d" % (s, self.order)
        s += ")"
        return s
