from openmm.app import *
from openmm import *
from openmm.unit import * 

from sys import stdout


pdb = PDBFile('ligand.pdb')
prmtop = AmberPrmtopFile('vacuum.prmtop')
#inpcrd = AmberInpcrdFile('receptor.inpcrd')
system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds,implicitSolvent=None)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
state = simulation.context.getState(getEnergy=True)
E_lig = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
print ('E(kJ/mol)', E_lig)
