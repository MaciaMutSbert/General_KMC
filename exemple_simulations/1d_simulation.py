from initialize_systems.__init__ import get_homogeneous_system
from core.__init__ import update_system, check_finish
from analysis_functions import get_trajectory
import numpy as np
import json
from molecules import Molecule

"Output file name"
output_file_name = '1d_simulation_trajectories_random.json'

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

state_energies = {'g_s': 0, 's_1': 2.5}          # eV Tetracene

relaxation_energies = {'g_s': 0, 's_1': 0.7}

transition_moment = np.array([1.0, 0, 0])        # a.u.  Tetracene value

generic_molecule = Molecule(state_energies, relaxation_energies, transition_moment)


"Morphology parameters"
order = 'ordered'                 # can take 'ordered' or 'disordered'
dimensions = [1000]                 # molecules per side. dimensionality = len(dimensions)
lattice_parameter = 1.0           # nm. Molecular site in an ordered system. Not used for disordered initialize_systems
# num_molecules = 0               # int. Number of molecules in the system. Not used in ordered initialize_systems

orientation = 'random'                              # orientation between molecules
reference_orientation = np.array([1.0, 0, 0])         # necessary only when (anti)parallelism is required

"""
Excitons.
Possible positions commands: 'random', 'first', 'last', 'centre', 'furthest'.
The dictionary takes the state as key and a list with the position of every exciton.
"""
excitons = {'s_1': ['centre']}


trajectories = []                   # list with the trajectories of all excitons
trajectories_steps = []
num_trajectories = 2000
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

    finished = False
    it = 0
    while finished is False:
        path, time = update_system(system)
        total_time.append(total_time[-1]+time)
        path_list.append(path)

        if check_finish(path_list) is True:
            break

        if total_time[-1] == 240:               # tallam a 4 cops la vida mitjana de l'exctió.
            break
        it = it+1
    print(it)
    print( 'j=', j)
    trajectories.append(get_trajectory(path_list, total_time, system))
    trajectories_steps.append(len(trajectories[-1]['positions']))


###########################################################################################################
# We list all the outputs in a dictionary system_information and write it in the output_file

system_information = {'conditions': conditions, 'state_energies': state_energies,
                      'relaxation_energies': relaxation_energies, 'transition_moment': list(transition_moment),
                      'order': order, 'dimensions': dimensions, 'lattice_parameter':  lattice_parameter,
                      'orientation':  orientation,  'reference_orientation': list(reference_orientation),
                      'excitons': excitons, 'type': system['type']}

"""
For an ordered system with parallel orientated dipoles a theorethical value for the diffusion constant and
diffusion length can be computed.
"""



output = {'system_information': system_information, 'trajectories': trajectories,
          'trajectory_steps': trajectories_steps}

with open(output_file_name, 'w') as write_file:
    json.dump(output, write_file)




