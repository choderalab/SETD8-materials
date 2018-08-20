# this is to repackage the FAH files generated with OpenMM 7.0, with OpenMM 6.3
# because we are using Core21 on FAH (OpenMM 6 core)

from simtk.openmm import XmlSerializer
import os

if not os.path.exists('package_OpenMM6/'):
    os.mkdir('package_OpenMM6')

structures = [
             'SET8', 'SET8_P', 'SET8_PH', 'SET8_MeP', 'SET8_SAM', 'SET8_SAH', 'SET8_P_SAM',
             'SET8_P_SAH', 'SET8_PH_SAM', 'SET8_PH_SAH', 'SET8_MeP_SAM', 'SET8_MeP_SAH',
             '4IJ8_SET8', '4IJ8_SET8_SAM', '4IJ8_SET8_SAH',
             '1ZKK_apo_SET8', '1ZKK_apo_SET8_P', '1ZKK_apo_SET8_PH', '1ZKK_apo_SET8_MeP',
             '1ZKK_apo_SET8_SAM', '1ZKK_apo_SET8_SAH', '1ZKK_apo_SET8_P_SAM', '1ZKK_apo_SET8_P_SAH',
             '1ZKK_apo_SET8_PH_SAM', '1ZKK_apo_SET8_PH_SAH', '1ZKK_apo_SET8_MeP_SAM', '1ZKK_apo_SET8_MeP_SAH'
             ]

def read_file(filename):
    with open(filename) as infile:
        return infile.read()

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

for structure in structures:

    if not os.path.exists('package/%s/system.xml' % structure):
        continue
    
    if not os.path.exists('package_OpenMM6/%s' % structure):
        os.mkdir('package_OpenMM6/%s' % structure)
    
    system = XmlSerializer.deserialize(read_file('package/%s/system.xml' % structure))
    write_file('package_OpenMM6/%s/system.xml' % structure, XmlSerializer.serialize(system))

    integrator = XmlSerializer.deserialize(read_file('package/%s/integrator.xml' % structure))
    write_file('package_OpenMM6/%s/integrator.xml' % structure, XmlSerializer.serialize(integrator))

    state = XmlSerializer.deserialize(read_file('package/%s/state.xml' % structure))
    write_file('package_OpenMM6/%s/state.xml' % structure, XmlSerializer.serialize(state))
