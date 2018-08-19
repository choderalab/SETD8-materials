# this is to repackage the FAH files generated with OpenMM 7.0, with OpenMM 6.3
# because we are using Core21 on FAH (OpenMM 6 core)

from simtk.openmm import XmlSerializer
import os

if not os.path.exists('package_OpenMM6/'):
    os.mkdir('package_OpenMM6')

structures = ['3F9W_B', '3F9W_C', '3F9X_A', 'TDIX_A', 'TDIY_B']

def read_file(filename):
    with open(filename) as infile:
        return infile.read()

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

for structure in structures:
    
    if not os.path.exists('package_OpenMM6/%s' % structure):
        os.mkdir('package_OpenMM6/%s' % structure)

    state = XmlSerializer.deserialize(read_file('%s/state0.xml' % structure))
    write_file('package_OpenMM6/%s/state0.xml' % structure, XmlSerializer.serialize(state))
