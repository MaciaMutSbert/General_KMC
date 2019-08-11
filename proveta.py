from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system, check_finish
from analysis_functions import get_trajectory
import matplotlib.pyplot as plt
import numpy as np
import json
from systems.molecules import Molecule


"Physical conditions"
conditions = {'temperature': 273.15,    # Kelvin
              'refractive_index': 1,    # adimensional
              'neighbourhood_radius': 1.1}  # nm. Maximum interaction distance

"""
Generic molecule initialization
Possible states: 
    'g_s': ground state 
    's_1': first singlet state 
    't_1': first triplet state 
All energies must be given in eV. By default initialized at g_s.
"""

state_energies = {'g_s': 0, 's_1': 0}          # eV Tetracene

relaxation_energies = {'g_s': 0, 's_1': 0.7}     # eV Tetracene

transition_moment = np.array([1.2, 0, 0])        # a.u.  Tetracene value

generic_molecule = Molecule(state_energies, relaxation_energies, transition_moment)


"Morphology parameters"
order = 'ordered'                 # can take 'ordered' or 'disordered'
dimensions = [140, 140]                 # molecules per side. dimensionality = len(dimensions)
lattice_parameter = 1.0           # nm. Molecular site in an ordered system. Not used for disordered systems
# num_molecules = 0               # int. Number of molecules in the system. Not used in ordered systems

orientation = 'parallel'                              # orientation between molecules
reference_orientation = np.array([1.0, 0, 0])         # necessary only when (anti)parallelism is required

"""
Excitons.
Possible positions commands: 'random', 'first', 'last', 'centre', 'furthest'.
The dictionary takes the state as key and a list with the position of every exciton.
"""
excitons = {'s_1': ['centre']}


trajectories = []                   # list with the trajectories of all excitons
num_trajectories = 1
num_steps = 1
for j in range(num_trajectories):
    key = str(len(dimensions))
    system = get_homogeneous_system['ordered'][key](conditions, generic_molecule, dimensions, lattice_parameter,
                                                    orientation, reference_orientation, excitons)

    """
    system is a dictionary with three keys:
        molecules: List of objects class Molecule
        conditions: dictionary with the physical conditions of the system such as temperature, refractive index...
        centre_indexes: list with the indexes of the excited molecules.
    """

    total_time = [0.0]
    path_list = []
    """
    Veure com definim el màxim d'iteracions.
    """
    finished = False
    for i in range(num_steps):
        " OJO. A LA FUNCIÓ NEIGHBOURHOOD HI HA UN TRUQUILLO QUE HEM DE CANVIAR EN CANVIAR LA DIMENSIONALITAT"
        path, time = update_system(system)
        total_time.append(total_time[-1]+time)
        path_list.append(path)
        print(i)
        if check_finish(path_list) is True:
            break
        if total_time[-1] == 480:
            break

    print('j=', j)
    trajectories.append(get_trajectory(path_list, total_time, system))


###########################################################################################################
# We list all the outputs in a dictionary system_information and write it in the output_file

system_information = {'conditions': conditions, 'state_energies': state_energies,
                      'relaxation_energies': relaxation_energies, 'transition_moment': list(transition_moment),
                      'order': order, 'dimensions': dimensions, 'lattice_parameter':  lattice_parameter,
                      'orientation':  orientation,  'reference_orientation': list(reference_orientation),
                      'excitons': excitons}

output = {'system_information': system_information, 'trajectories': trajectories, 'steps': num_steps}

with open("results_file_test_neighbour.json", 'w') as write_file:
    json.dump(output, write_file)




