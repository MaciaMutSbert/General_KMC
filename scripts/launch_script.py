from initialize_systems import get_system
from core import update_system, check_finish
from analysis import update_trajectory
import numpy as np
import json
from molecules import Molecule


#######################################################################################################################

output_file_name = 'example_1d_simulation_mu_0-6.json'       # name of the output file where the trajectories will be saved
                                                        # .json format

#######################################################################################################################

"""
Generic molecule initialization
Possible states: 
    'gs': ground state 
    's1': first singlet state 
All energies must be given in eV. By default initialized at gs.
"""
state_energies = {'gs': 0, 's1': 2.5}                           # excitation energies of the electronic states (eV)
reorganization_energies = {'gs': 0, 's1': 0.7}                  # reorganization energies of the states (eV)
transition_dipole_moment = np.array([0.6, 0, 0])                # transition dipole moment of the molecule (a.u)

molecule = Molecule(state_energies=state_energies,
                    reorganization_energies=reorganization_energies,
                    transition_moment=transition_dipole_moment)

#######################################################################################################################

# physical conditions of the system (as a dictionary)
conditions = {'temperature': 273.15,            # temperature of the system (K)
              'refractive_index': 1,            # refractive index of the material (adimensional)
              'neighbourhood_radius': 1.1}      # maximum interaction distance (nm)

#######################################################################################################################

trajectories = []                               # list with the trajectories of all excitons
num_trajectories = 300                         # number of trajectories that will be simulated
max_steps = 10000                               # maximum number of steps for trajectory allowed

for j in range(num_trajectories):

    system = get_system(conditions=conditions,                  # dictionary with the physical conditions of the simulation
                        molecule=molecule,
                        lattice={'dimensions': [1000],          # dictionary lattice. dimensions: supercell of the crystall
                                 'lattice_parameter': [1.0]},     # lattice parameter (nm)
                        amorphous=None,                         # dictionary with a set of parameters for a disordered system
                        orientation='parallel',                 # orientation of the molecules: 'parallel', 'antiparallel', 'random'
                        initial_excitation={'s1': ['centre']})  # intial excitation of the system (excited states and the positions of the excitons in the lattice)

#    system is a dictionary with three keys:
#       molecules: List of objects class Molecule
#       conditions: dictionary with the physical conditions of the system such as temperature, refractive index...
#       lattice/amorphous: dictionary with the morphology information
#   Tricky entrances:
#       centres: list with the indexes of the excited molecules.
#       type (label for the system)
#######################################################################################################################

    path = []                          # list with the path followed by the system

    trajectory = {'time': [0.0], 'n': [], 'positions': [], 'process': []}
    # trajectory of the system. Gives the number of excited states, its positions and the process occurred at each time
    # first output of the programm

    for i in range(max_steps):

        change_step, step_time = update_system(system)
        # returns the process occurred {donor, process, acceptor} and the duration of this process
        # also modifies system updating the dictionary with the information of change_step.

        path.append(change_step)
        # print(i)

        update_trajectory(trajectory, change_step, step_time, system)
        # the dictionary trajectory is updated with the information of the occurred process

        if check_finish(system) is True:
            # checks if all excitons have decayed
            break

        if i == max_steps-1:
            # warns the user if the maximum number of steps has been reached.
            print('Maximum number of steps reached!!')

    print('j= ', j)
    trajectories.append(trajectory)


###########################################################################################################
# We collect all the outputs in a dictionary system_information and write it in the output_file

# s'ha de veure com s'enllaça la informació que donam aquí amb la que es defineix a l'input #############

system_information = {'conditions': system['conditions'], 'lattice': system['lattice'],
                      'orientation': 'parallel', 'type': system['type'], 'excitons': {'s1': ['centre']},
                      'state_energies': state_energies, 'reorganization_energies': reorganization_energies,
                      'transition_moment': list(transition_dipole_moment)}


output = {'system_information': system_information,
          'trajectories': trajectories, 'max_steps': max_steps}

with open(output_file_name, 'w') as write_file:
    json.dump(output, write_file)




