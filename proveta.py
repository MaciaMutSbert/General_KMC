from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system, check_finish
from result_analysis import path_distance, final_position, x_y_splitter
import matplotlib.pyplot as plt
import numpy as np
from systems.molecules import Molecule


conditions = {'temperature': 273.15, 'refractive_index': 2}
lattice_parameter = 1   # nm
dimensions = [100, 100]    # nm


final_positions = []
time_list = []
exciton_distance = []

for i in range(1):
    system = get_homogeneous_system(conditions, dimensions=dimensions, lattice_parameter=lattice_parameter)
    total_time = 0
    distance = 0

    path_list = []

    general_molecule = Molecule([0, 0, 0], state=1)
    decay_rates = general_molecule.decay_rate()
    rate = decay_rates['Singlet_radiative_decay_rate']
    time_max = 1 / rate
    while total_time < time_max:
        path, time = update_system(system)
        total_time += time
        distance += path_distance(path, system)
        path_list.append(path)

        if check_finish(path_list) is True:
            break

    final_positions.append(final_position(path_list, system))
    time_list.append(total_time)
    exciton_distance.append(distance)

final_x_list, final_y_list = x_y_splitter(final_positions)

print(final_positions)
print(exciton_distance)
print(time_list)
"""
                    Results

plt.plot(final_x_list, final_y_list, 'ro', label='Position')
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Final positions')
plt.show()
"""