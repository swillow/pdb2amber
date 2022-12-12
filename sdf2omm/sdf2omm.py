from openff.toolkit import ForceField, Molecule, Topology 
from openmm import unit as openmm_unit
import openmm

fname_ligand = 'ligand.sdf'
ligand = Molecule.from_file (fname_ligand)
top = Topology.from_molecules (ligand)
sage_ff14sb = ForceField("openff-2.0.0.offxml")
interchange = sage_ff14sb.create_interchange(top)
omm_system = interchange.to_openmm()
omm_top = interchange.to_openmm_topology()

# Construct and configure a Langevin integrator at 300 K with an appropriate friction constant and time-step
integrator = openmm.LangevinIntegrator(
    300 * openmm_unit.kelvin,
    1 / openmm_unit.picosecond,
    0.002 * openmm_unit.picoseconds,
)


# Combine the topology, system, integrator and initial positions into a simulation
simulation = openmm.app.Simulation(omm_top, omm_system, integrator)
simulation.context.setPositions(top.get_positions().to_openmm())

state = simulation.context.getState(getEnergy=True, getForces=True)
print ("Energy", state.getPotentialEnergy())
