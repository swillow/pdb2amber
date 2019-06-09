from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

pdb = PDBFile('receptor.pdb')
forcefield = ForceField('amber14-all.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
state = simulation.context.getState(getEnergy=True)
E = state.getPotentialEnergy().value_in_unit(kilocalorie_per_mole)

print ("ENERGY ", E)

