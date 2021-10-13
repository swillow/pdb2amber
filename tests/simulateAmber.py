try:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *    
except:
    from openmm.app import *
    from openmm import *
    from openmm.unit import * 

from sys import stdout

prmtop = AmberPrmtopFile('receptor.prmtop')
inpcrd = AmberInpcrdFile('receptor.inpcrd')
system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
state = simulation.context.getState(getEnergy=True)
E = state.getPotentialEnergy().value_in_unit(kilocalorie_per_mole)

print ("ENERGY ", E)

