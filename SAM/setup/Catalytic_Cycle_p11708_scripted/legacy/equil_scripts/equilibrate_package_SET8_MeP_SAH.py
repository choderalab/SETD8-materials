from simtk.openmm import *
from simtk.openmm.app import *
import simtk.unit as u

### set structure name for every individual script
structure_name = 'SET8_MeP_SAH'

### parameters - these made to match Ensembler, except the barostat_frequency
# we are going with default OpenMM of 25, rather than 50 in Ensembler
verlet_timestep = 1.0 * u.femtoseconds
timestep = 2.0 * u.femtoseconds

cutoff = 0.9 * u.nanometers

friction = 1.0 / u.picoseconds
equil_friction = 20.0 / u.picoseconds

output_frequency = 1000
n_steps = 5000000

temperature = 300.0 * u.kelvin
pressure = 1.0 * u.atmospheres
barostat_frequency = 25
###

if not os.path.exists('../equil/'):
    os.mkdir('../equil')
if not os.path.exists('../package/'):
    os.mkdir('../package')
if not os.path.exists('../package/%s' % structure_name):
    os.mkdir('../package/%s' % structure_name)

pdb_filename = "../with_solvent/%s.pdb" % structure_name
out_pdb_filename = "../equil/%s.pdb" % structure_name
dcd_filename = "../equil/%s.dcd" % structure_name
log_filename = "../equil/%s.log" % structure_name

system_filename = "../package/%s/system.xml" % structure_name
integrator_filename = "../package/%s/integrator.xml" % structure_name
state_filename = "../package/%s/state.xml" % structure_name

print("Loading PDB...")
pdb = PDBFile(pdb_filename)

topology = pdb.topology
positions = pdb.positions

ff = ForceField('amber99sbildn.xml', 'tip3p.xml', '../parameters/gaff.xml', '../parameters/ffptm.xml', '../parameters/SAM.xml', '../parameters/SAH.xml')

platform = Platform.getPlatformByName('CUDA')
# for minimization do not specify properties - i.e. go with SINGLE precision - it does not minimize enough (i.e. energy decrease is only small - tested for SET8 system)
# but for equilibration and then production on FAH - MIXED precision
properties = {'CudaPrecision': 'mixed'}

# need to minimize with CutoffPeriodic (reaction field) to avoid blowing up
print("Preparing CutoffPeriodic minimization system...")
system = ff.createSystem(topology, nonbondedMethod=CutoffPeriodic, nonbondedCutoff=cutoff, constraints=HBonds)

print("Preparing simulation for minimization...")
integrator = VerletIntegrator(verlet_timestep)
simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)
print("Initial energy is %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy()
print("Energy after minimization is %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
positions = simulation.context.getState(getPositions=True).getPositions()
del simulation

print("Preparing new PME system for 10ns equilibration run...")
system = ff.createSystem(topology, nonbondedMethod=PME, nonbondedCutoff=cutoff, constraints=HBonds)
integrator = LangevinIntegrator(temperature, equil_friction, timestep)
system.addForce(MonteCarloBarostat(pressure, temperature, barostat_frequency))

print("Preparing simulation for equilibration 10ns run...")
simulation = Simulation(topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)

print("Initial energy is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy()
print("Energy after minimization is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))

simulation.context.setVelocitiesToTemperature(temperature)

print("Appending reporters...")
simulation.reporters.append(DCDReporter(dcd_filename, output_frequency))
simulation.reporters.append(PDBReporter(out_pdb_filename, n_steps))
simulation.reporters.append(StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True, potentialEnergy=True, temperature=True, density=True))

print("Running 10ns equilibration simulation steps...")
simulation.step(n_steps)
print('Done with the equilibration!')
print("Final energy is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
positions = simulation.context.getState(getPositions=True).getPositions()
vectors = simulation.context.getState().getPeriodicBoxVectors()
del simulation

# Do packaging
print('Initiating packaging...')
def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

print('Preparing new integrator and updating system with new PeriodicBoxVectors...')
system.setDefaultPeriodicBoxVectors(vectors[0], vectors[1], vectors[2])
integrator = LangevinIntegrator(temperature, friction, timestep)

print("Preparing new simulation...")
simulation = Simulation(topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)
print("Performing 1000 simulation steps...")
simulation.step(output_frequency)
print("Updating system with new PeriodicBoxVectors...")
vectors = simulation.context.getState().getPeriodicBoxVectors()
system.setDefaultPeriodicBoxVectors(vectors[0], vectors[1], vectors[2])

print("Writing system and integrator...")
write_file(system_filename, XmlSerializer.serialize(system))
write_file(integrator_filename, XmlSerializer.serialize(integrator))

print("Writing state...")
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
write_file(state_filename, XmlSerializer.serialize(state))
print('All done!')
