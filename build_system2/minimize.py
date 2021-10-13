import numpy as np
import sys
try:
    from simtk.openmm.app.internal import amber_file_parser
    from simtk.openmm.app import *
    from simtk.openmm.vec3 import Vec3
    from simtk.openmm import *
    import simtk.unit as unit
except:
    from openmm.app.internal import amber_file_parser
    from openmm.app import *
    from openmm.vec3 import Vec3
    from openmm import *
    import openmm.unit as unit
    
from _inpcrd import *


boltz_kcal = 8.31446261815324/4184  # kcal/(Kmol)


def _init_neigh_list(n_neigh):
    '''
    n_neigh: int
    generate neighboring cells with rcut_cell = rcut/n_neigh,
    rcut is the cutoff distance

    Returns
    -------
    np.array ([-n_neigh:n_neigh, -n_neigh:n_neigh, -n_neigh:n_neigh])
    '''
    neigh_list = []

    for ic in range(-n_neigh, n_neigh+1):
        for jc in range(-n_neigh, n_neigh+1):
            for kc in range(-n_neigh, n_neigh+1):
                neigh_list.append([ic, jc, kc])

    return np.array(neigh_list)


def _init_cell_list(ndim):
    '''
    ndim: np.array (3)
     the number of cell in x, y, z directions

    using a cell (linked) list method,

    Returns
    -------
    cell_list: dictionary
      each cell will store the atom lists
    '''

    cell_list = {}
    for ic in range(ndim[0]):
        for jc in range(ndim[1]):
            for kc in range(ndim[2]):
                cell_list[(ic, jc, kc)] = []

    return cell_list


def pot_grd_bond(pos, bonds):
    """
    pos : np.array (natm, 3)
    bonds: list [(iatom, jatom, k, dist)]
    """

    en = 0.0
    grd = np.zeros(pos.shape, dtype=np.float32)

    for (atom1, atom2, k, r0) in bonds:

        dij = pos[atom1] - pos[atom2]
        r = np.sqrt(np.einsum('i,i', dij, dij))
        en += k*(r-r0)*(r-r0)
        dedr = 2.0*k*(r-r0)
        de = dedr/r
        gij = de*dij
        grd[atom1] += gij
        grd[atom2] -= gij

    return en, grd


def pot_grd_angle(pos, angles):
    """
    pos : np.array (natm, 3)
    angles: list [(atom1, atom2, atom3, k, dist)]
    """

    en = 0.0
    grd = np.zeros(pos.shape, dtype=np.float32)

    for (atom1, atom2, atom3, k, r0) in angles:

        pos_ij = pos[atom1] - pos[atom2]
        pos_kj = pos[atom3] - pos[atom2]
        pos_p = np.cross(pos_kj, pos_ij)

        rij2 = np.einsum('i,i', pos_ij, pos_ij)
        rkj2 = np.einsum('i,i', pos_kj, pos_kj)
        rp = np.sqrt(np.einsum('i,i', pos_p, pos_p))
        cs = np.einsum('i,i', pos_ij, pos_kj)/np.sqrt(rij2*rkj2)
        cs = min(1.0, max(-1.0, cs))
        theta = np.arccos(cs)
        delta = theta - r0
        en += k*delta*delta

        dedt = 2.0*k*delta
        termi = -dedt/(rij2*rp)
        termk = dedt/(rkj2*rp)

        gi = termi*np.cross(pos_ij, pos_p)
        gk = termk*np.cross(pos_kj, pos_p)

        grd[atom1] += gi
        grd[atom3] += gk
        grd[atom2] -= (gi+gk)

    return en, grd


def pot_grd_torsion(pos, torsions_cos_phase, torsions_kvals):

    en = 0.0
    grd = np.zeros(pos.shape)

    for (ia, ib, ic, id) in torsions_kvals:

        r_ba = pos[ib] - pos[ia]
        r_cb = pos[ic] - pos[ib]
        r_dc = pos[id] - pos[ic]
        r_t = np.cross(r_ba, r_cb)
        r_u = np.cross(r_cb, r_dc)
        r_tu = np.cross(r_t, r_u)
        rt2 = np.einsum('i,i', r_t, r_t)
        ru2 = np.einsum('i,i', r_u, r_u)
        rcb = np.einsum('i,i', r_cb, r_cb)
        rtru = np.sqrt(rt2*ru2)

        kval = torsions_kvals[(ia, ib, ic, id)]
        cs0 = torsions_cos_phase[(ia, ib, ic, id)]

        cs = np.empty(6, dtype=np.float32)
        sn = np.empty(6, dtype=np.float32)
        phi = np.empty(6, dtype=np.float32)
        dphi = np.empty(6, dtype=np.float32)
        cs[0] = np.einsum('i,i', r_t, r_u)/rtru
        sn[0] = np.einsum('i,i', r_cb, r_tu)/(rcb*rtru)
        phi[0] = 1.0 + cs[0]*cs0[0]
        dphi[0] - sn[0]*cs0[0]
        for kk in range(1, 6):
            cs[kk] = cs[0]*cs[kk-1] - sn[0]*sn[kk-1]
            sn[kk] = cs[0]*sn[kk-1] + sn[0]*cs[kk-1]
            phi[kk] = 1.0 + cs[kk]*cs0[kk]
            dphi[kk] = -sn[kk]*cs0[kk]

        en += np.einsum('i,i', kval, phi)
        dedphi = np.einsum('i,i', kval, dphi)

        dedt = dedphi*np.cross(r_t, r_cb)/(rt2*rcb)
        dedu = -dedphi*np.cross(r_u, r_cb)/(ru2*rcb)

        r_ac = pos[ia] - pos[ic]
        r_bd = pos[ib] - pos[id]

        grd[ia] += np.cross(dedt, r_cb)
        grd[ib] += np.cross(dedt, r_ac) + np.cross(dedu, r_dc)
        grd[ic] += np.cross(dedt, r_ba) + np.cross(dedu, r_bd)
        grd[id] += np.cross(dedu, r_cb)

    return en, grd


def pot_grd_nonbond(pos, nb_parms, nb_pairs):
    """
    nonbonding interactions between only overlapped atoms (rij < 2.0 or 1.5)
    """
    grd = np.zeros(pos.shape)
    cutoff = 10.0  # A
    en = 0.0

    for ia, ja in nb_pairs:
        rmin2_i, eps_i, q_i, m_i = nb_parms[ia]
        rmin2_j, eps_j, q_j, m_j = nb_parms[ja]

        pij = pos[ia] - pos[ja]
        rij = np.sqrt(np.einsum('i,i', pij, pij))

        if rij < cutoff:
            rmin = rmin2_i + rmin2_j
            eps = np.sqrt(eps_i*eps_j)
            qij = q_i*q_j

            rtmp = rmin/rij
            rtmp2 = rtmp*rtmp
            rtmp6 = rtmp2*rtmp2*rtmp2

            enlj = eps*rtmp6*(rtmp6 - 2.0)
            delj = eps*rtmp6*(rtmp6-1.0)*(-12.0/rij)

            """
            if rij > 1.0:
                print("TOO OVERLAP ", ia, ja, enlj, rij)
                print("pos_i, pos_j ", pos[ia], pos[ja])
            """
            enqq = qij/rij
            deqq = -enqq/rij

            gij = (delj + deqq)*pij/rij

            grd[ia] += gij
            grd[ja] -= gij

            en += (enlj + enqq)

    return en, grd


def pot_grd_nonbond14(pos, nb14_parms, moving_atoms):

    cutoff = 10.0
    scee = 1.0/1.4
    scnb = 1.0/2.0

    grd = np.zeros(pos.shape)
    en = 0.0

    for (atom1, atom4, chg14, rmin14, eps14) in nb14_parms:
        if atom1 not in moving_atoms and atom4 not in moving_atoms:
            continue
        p14 = pos[atom1] - pos[atom4]
        r14 = np.sqrt(np.einsum('i,i', p14, p14))

        if r14 < cutoff:
            rtmp = rmin14/r14
            rtmp2 = rtmp*rtmp
            rtmp6 = rtmp2*rtmp2*rtmp2

            enlj = scnb*eps14*rtmp6*(rtmp6 - 2.0)
            delj = scnb*eps14*rtmp6*(rtmp6-1.0)*(-12.0/r14)

            enqq = scee*chg14/r14
            deqq = -enqq/r14

            g14 = (delj+deqq)*p14/r14
            grd[atom1] += g14
            grd[atom4] -= g14

            en += (enlj + enqq)

    return en, grd


def read_amber_prmtop(FN_prmtop):

    prm = amber_file_parser.PrmtopLoader(FN_prmtop)

    charges = prm.getCharges()
    masses = prm.getMasses()
    znums = np.array(prm._raw_data['ATOMIC_NUMBER'], dtype=int)

    # LJ
    nonBondTerms_omm = np.array(prm.getNonbondTerms())

    nonBondTerms = []
    nonBondExcls = {}
    for ia, (rmin2, eps) in enumerate(nonBondTerms_omm):
        rmin2 *= 10.0  # nm -->A
        eps /= 4.184  # kJ/mol --> kcal/mol
        chg = charges[ia]*18.2223  # CHARGECON --> kcal/mol
        mass = masses[ia]/418.4  # at.weight --> kcal/mol * (ps/A)
        nonBondTerms.append((float(rmin2),
                             float(eps),
                             float(chg),
                             float(mass)))
        nonBondExcls[ia] = []

    nonBondTerms = np.array(nonBondTerms)

    # Bond
    bonds_omm = prm.getBondsNoH() + prm.getBondsWithH()
    bonds = []
    forceConstConversionFactor = (unit.kilojoule_per_mole/(unit.nanometer*unit.nanometer)
                                  ).conversion_factor_to(unit.kilocalorie_per_mole/(unit.angstroms*unit.angstroms))

    for atom1, atom2, k, dist in bonds_omm:
        k *= forceConstConversionFactor
        dist *= 10.0
        bonds.append((int(atom1), int(atom2), float(k), float(dist)))

        nonBondExcls[atom1].append(atom2)
        nonBondExcls[atom2].append(atom1)

    # Angle
    angles_omm = prm.getAngles()
    angles = []
    forceConstConversionFactor = unit.kilojoule_per_mole.conversion_factor_to(
        unit.kilocalorie_per_mole)
    for (atom1, atom2, atom3, k, dist) in angles_omm:
        k *= forceConstConversionFactor
        angles.append((int(atom1),
                       int(atom2),
                       int(atom3),
                       float(k),
                       float(dist)))

        nonBondExcls[atom1].append(atom3)
        nonBondExcls[atom3].append(atom1)

    # Dihedrals
    torsions_omm = prm.getDihedrals()
    forceConstConversionFactor = unit.kilojoule_per_mole.conversion_factor_to(
        unit.kilocalorie_per_mole)

    torsions_cos_phase = {}
    torsions_kvals = {}
    for (atom1, atom2, atom3, atom4, k, phase, periodicity) in torsions_omm:
        k *= forceConstConversionFactor

        if (atom1, atom2, atom3, atom4) not in torsions_kvals:
            torsions_kvals[(atom1, atom2, atom3, atom4)
                           ] = np.zeros(6, dtype=np.float32)
            torsions_cos_phase[(atom1, atom2, atom3, atom4)
                               ] = np.zeros(6, dtype=np.float32)

        torsions_cos_phase[(atom1, atom2, atom3, atom4)
                           ][periodicity-1] = np.cos(phase)
        torsions_kvals[(atom1, atom2, atom3, atom4)
                       ][periodicity-1] = k

    # NonBond14
    nonBond14_omm = prm.get14Interactions()
    nonBond14 = []
    for (atom1, atom4, chg14, rmin14, eps14, scee, scnb) in nonBond14_omm:
        eps14 /= 4.184  # kJ/mol --> kcal/mol
        rmin14 *= 10.0  # nm -> A

        nonBond14.append((int(atom1),
                          int(atom4),
                          float(chg14),
                          float(rmin14),
                          float(eps14)))

        nonBondExcls[atom1].append(atom4)
        nonBondExcls[atom4].append(atom1)

    return bonds, angles, torsions_kvals, torsions_cos_phase, \
        nonBond14, nonBondTerms, znums, nonBondExcls


def nonbond_pairs(pos, rcut_half, neigh_list, nonBondExcls,
                  rcut_ovlp,
                  znums=None):
    """
    """

    unitcell_min = pos.min(axis=0)
    unitcell_max = pos.max(axis=0)
    box = unitcell_max - unitcell_min + 1.0

    # Cell Linked List
    ndim = np.array(box//rcut_half, dtype=np.int32)
    rcut_cell = box/ndim
    cell_list = _init_cell_list(ndim)

    for ia, ri in enumerate(pos):
        pi = ri - unitcell_min
        ic, jc, kc = np.array(pi//rcut_cell, dtype=np.int32)

        if (ic, jc, kc) in cell_list:
            cell_list[(ic, jc, kc)].append(ia)
        else:
            print('icell error ', ic, jc, kc, ndim)
            sys.exit(1)

    # Generate nonbond_pairs
    nb_pairs = []
    if znums is None:
        for ia, ri in enumerate(pos):
            pi = ri - unitcell_min
            icel = np.array(pi//rcut_cell, dtype=np.int32)
            ncel = icel + neigh_list

            for ic, jc, kc in ncel:
                if (ic, jc, kc) in cell_list:
                    for ja in cell_list[(ic, jc, kc)]:
                        if ja < ia:
                            if ja in nonBondExcls[ia]:
                                continue
                            pij = pos[ia] - pos[ja]
                            rij = np.sqrt(np.einsum('i,i', pij, pij))

                            if rij < rcut_ovlp:
                                nb_pairs.append((ia, ja))
    else:
        for ia, ri in enumerate(pos):
            if znums[ia] == 1:
                continue
            pi = ri - unitcell_min
            icel = np.array(pi//rcut_cell, dtype=np.int32)
            ncel = icel + neigh_list

            for ic, jc, kc in ncel:
                if (ic, jc, kc) in cell_list:
                    for ja in cell_list[(ic, jc, kc)]:
                        if ja < ia:
                            if znums[ja] == 1:
                                continue
                            if ja in nonBondExcls[ia]:
                                continue

                            pij = pos[ia] - pos[ja]
                            rij = np.sqrt(np.einsum('i,i', pij, pij))

                            if rij < rcut_ovlp:
                                nb_pairs.append((ia, ja))

    return nb_pairs


class PotGrad(object):

    def __init__(self, fname_prmtop):

        self.bonds, \
            self.angles, \
            self.torsions_kvals, \
            self.torsions_cos_phase, \
            self.nonBond14, \
            self.nonBondTerms, \
            self.znums, \
            self.nonBondExcls = read_amber_prmtop(fname_prmtop)

        self.rcut = 8.0
        self.rcut_half = self.rcut/2
        self.neigh_list = _init_neigh_list(1)

    def run(self, pos, l_with_H=False):
        """
        pos : np.array ( (natm, 3) )
        """

        # Generate NonBond Pair Lists
        if l_with_H:
            rcut_ovlp = 1.5
            nb_pairs = nonbond_pairs(pos, self.rcut_half,
                                     self.neigh_list,
                                     self.nonBondExcls,
                                     rcut_ovlp)
        else:
            rcut_ovlp = 2.0
            nb_pairs = nonbond_pairs(pos, self.rcut_half,
                                     self.neigh_list,
                                     self.nonBondExcls,
                                     rcut_ovlp,
                                     znums=self.znums)

        if len(nb_pairs) == 0:
            print("No overlapped Atoms")
            sys.exit(1)

        ovlp_atoms = []
        for ia, ja in nb_pairs:
            pij = pos[ia] - pos[ja]
            rij = np.sqrt(np.einsum('i,i', pij, pij))
            ovlp_atoms.append((ia, ja, rij))

        def sortFunc(ovlp_atom):
            ia, ja, rij = ovlp_atom
            return rij

        ovlp_atoms.sort(key=sortFunc)
        print('mimimum distance of ovlapped atoms ', ovlp_atoms[0])

        moving_atoms = []
        for ia, ja, rij in ovlp_atoms:
            if ia not in moving_atoms:
                moving_atoms.append(ia)
            if ja not in moving_atoms:
                moving_atoms.append(ja)

        print('moving atoms ', len(moving_atoms))

        moving_atoms.sort()

        en = 0.0
        grd = np.zeros(pos.shape)

        #
        en_tmp, grd_tmp = pot_grd_bond(pos, self.bonds)
        en += en_tmp
        grd += grd_tmp
        print('bonds ', en_tmp)
        #
        en_tmp, grd_tmp = pot_grd_angle(pos, self.angles)
        en += en_tmp
        grd += grd_tmp
        print('angles ', en_tmp)
        #
        en_tmp, grd_tmp = pot_grd_torsion(pos,
                                          self.torsions_cos_phase,
                                          self.torsions_kvals)

        en += en_tmp
        grd += grd_tmp
        print('torsions ', en_tmp)

        # pot_grd_nonbond14(pos, nb14_parms)
        en_tmp, grd_tmp = pot_grd_nonbond14(pos,
                                            self.nonBond14, moving_atoms)
        en += en_tmp
        grd += grd_tmp
        print('nbond14 ', en_tmp)
        #
        en_tmp, grd_tmp = pot_grd_nonbond(pos,
                                          self.nonBondTerms,
                                          nb_pairs)

        en += en_tmp
        grd += grd_tmp
        print('nonbonds ', en_tmp)

        return en, grd, moving_atoms


class MINI (object):

    def __init__(self, fname_prmtop, fname_inpcrd):
        """
        temp in K
        dt (time step) in ps 
        """
        self.pos = AmberInpcrdFile(
            fname_inpcrd).getPositions(asNumpy=True).value_in_unit(unit.angstroms)
        self.pot_grd = PotGrad(fname_prmtop)
        self.at_mass = self.pot_grd.nonBondTerms[:, 3]

    def run(self, l_with_H=False, nstep=1000, dt=2.0e-3, dmax=0.05):
        """
        dt = ps
        """

        en_pot, grd, moving_atoms = self.pot_grd.run(self.pos, l_with_H)

        # en_kin = self.kinetic_energy()

        nprint = 2

        for istep in range(nstep):

            for ia, gi in enumerate(grd):
                # if ia in moving_atoms:
                dpi = - dt*dt*gi/self.at_mass[ia]
                dri = np.sqrt(np.einsum('i,i', dpi, dpi))
                if dri > dmax:
                    scale = dmax/dri
                    self.pos[ia] += dpi*scale
                else:
                    self.pos[ia] += dpi

            if (istep + 1) % nprint == 0:
                print(istep+1, en_pot)
                print_inpcrd('updated.inpcrd', self.pos, vel=None, box=None)

            en_pot, grd, moving_atoms = self.pot_grd.run(self.pos, l_with_H)


if __name__ == '__main__':
    """
    1) remove the overlapped heavy atoms
    mini.run(l_with_H=False)
    2) remove the overlapped atoms
    mini.run(l_with_H=True)
    """
    #fname_inpcrd = 'prt_mem_wat_nacl.inpcrd'
    fname_inpcrd = 'restart.inpcrd'
    fname_prmtop = 'prt_mem_wat_nacl.prmtop'

    mini = MINI(fname_prmtop, fname_inpcrd)
    #mini.run(l_with_H=False)
    mini.run(l_with_H=True)
