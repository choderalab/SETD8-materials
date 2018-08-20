from simtk import unit as u

timestep = 2.0 * u.femtoseconds
equil_timestep = 1.0 * u.femtoseconds

padding = 1.3 * u.nanometers
cutoff = 0.95 * u.nanometers

friction = 0.25 / u.picoseconds
equil_friction = 2.0 / u.picoseconds

output_frequency = 1000
n_steps = 5000000

temperature = 300. * u.kelvin
pressure = 1.0 * u.atmospheres
barostat_frequency = 25

discard_steps = 50000
