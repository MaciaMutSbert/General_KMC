from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system, check_finish
from result_analysis import path_distance, final_position, x_y_splitter
import matplotlib.pyplot as plt
import numpy as np


"""
In conditions the physical conditions of the system and molecule are defined.
Global condition: temperature
Material parameters: refractive index (n), orientational factor, k
Molecule parameters:    singlet (excitation) energy, must be given in eV
                        transition dipole, must be given in atomic units
                        characteristic_length, in nm
Absorption and emission spectra deviation  (considered as gaussian)       (unclassified parameter)
"""
conditions = {'temperature': 273.15, 'refractive_index': 2, 'orientational_factor_1': 2, 'neighbourhood_radius': 0.5,
              'a_e_spectra_deviation': 0.3,
              'molecule_type': 1, 'singlet_energy': 2.5, 'transition_dipole': 1, 'characteristic_length': 10**-8}

lattice_parameter = 0.45       # nm
dimensions = [90, 90]          # nm


final_positions = []
time_list = []
exciton_distance = []

for i in range(1):
    system = get_homogeneous_system(conditions, dimensions=dimensions, lattice_parameter=lattice_parameter)
    total_time = 0
    distance = 0

    path_list = []

    finished = False
    it = 0
    while finished is False:
        path, time = update_system(system)
        total_time += time
        print(total_time)
        distance += path_distance(path, system)
        print(distance)
        path_list.append(path)
        it += 1

        finished = check_finish(path_list)
    final_positions.append(final_position(path_list, system))
    time_list.append(total_time)
    exciton_distance.append(distance)

final_x_list, final_y_list = x_y_splitter(final_positions)


diffusion_length = np.average(np.array(exciton_distance))
life_time = np.average(np.array(time_list))

print(exciton_distance)
print(diffusion_length)
print(life_time)
