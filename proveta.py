from initialize_systems import get_system
from core.__init__ import update_system, check_finish
from analysis_functions import get_trajectory
import numpy as np
import json
from molecules import Molecule


#######################################################################################################################

output_file_name = 'example_1d_simulation.json'         # name of the output file where the trajectories will be saved
                                                        # .json format

#######################################################################################################################

"""
Generic molecule initialization
Possible states: 
    'gs': ground state 
    's1': first singlet state 
All energies must be given in eV. By default initialized at gs.
"""

molecule = Molecule(state_energies={'gs': 0, 's1': 2.5},                # excitation energies of the electronic states (eV)
                    reorganization_energies={'gs': 0, 's1': 0.7},       # reorganization energies of the states (eV)
                    transition_moment=np.array([1.2, 0, 0]))            # transition dipole moment of the molecule (a.u)

#######################################################################################################################

# physical conditions of the system (as a dictionary)
conditions = {'temperature': 273.15,            # temperature of the system (K)
              'refractive_index': 1,            # refractive index of the material (adimensional)
              'neighbourhood_radius': 1.1}      # maximum interaction distance (nm)

#######################################################################################################################

trajectories = []                               # list with the trajectories of all excitons
num_trajectories = 1000                         # number of trajectories that will be simulated
max_steps = 10000                               # maximum number of steps for trajectory allowed

for j in range(num_trajectories):

    system = get_system(conditions=conditions,                  # dictionary with the physical conditions of the simulation
                        molecule=molecule,
                        lattice={'dimensions': [1000],          # dictionary lattice. dimensions: supercell of the crystall
                                 'lattice_parameter': 1.0},     # lattice parameter (nm)
                        amorphous=None,                         # dictionary with a set of parameters for a disordered system
                        orientation='parallel',                 # orientation of the molecules: 'parallel', 'antiparallel', 'random'
                        initial_excitation={'s1': ['centre']})  # intial excitation of the system (excited states and the positions of the excitons in the lattice)


#    system is a dictionary with three keys:
#       molecules: List of objects class Molecule
#        conditions: dictionary with the physical conditions of the system such as temperature, refractive index...
#        centre_indexes: list with the indexes of the excited molecules.
#######################################################################################################################

    time = [0.0]
    path = []

    for i in range(max_steps):

        change_step, step_time = update_system(system)

        time.append(time[-1] + step_time)
        path.append(change_step)
        print(i)

        if check_finish(system) is True:
            break

        if i == max_steps-1:
            print('Maximum number of steps reached!!')

    print('j= ', j)
    trajectories.append(get_trajectory(path, time, system))


###########################################################################################################
# We collect all the outputs in a dictionary system_information and write it in the output_file

system_information = {'conditions': conditions, 'state_energies': state_energies,
                      'relaxation_energies': relaxation_energies, 'transition_moment': list(transition_moment),
                      'order': order, 'dimensions': dimensions, 'lattice_parameter':  lattice_parameter,
                      'orientation':  orientation,  'reference_orientation': list(reference_orientation),
                      'excitons': excitons}

output = {'system_information': system_information, 'trajectories': trajectories, 'steps': max_steps}

with open("results_file_test_neighbour.json", 'w') as write_file:
    json.dump(output, write_file)




