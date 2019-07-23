from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system, check_finish, get_centers
from result_analysis import initial_position, final_position, x_y_splitter, moved_length
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


lattice_parameter = 0.45           # nm
dimensions = [26.1, 26.1]          # nm
num_dimensions = 2
excitons = {'number': 1, 'positions': 1653}


conditions = {'temperature': 273.15, 'refractive_index': 2, 'neighbourhood_radius': 0.46,
              'a_e_spectra_deviation': 0.3, 'a_e_centre_shift': 1.65,
              'singlet_energy': 2.5, 'transition_dipole': 1, 'characteristic_length': 10 ** -8}

# Lists for further analysis
initial_positions = []
final_positions = []
time_list = []
life_time = []
diffusion_length = []
exciton_lengths = []

# Memory of the system
memory = {}

for j in range(100):
    system = get_homogeneous_system(conditions, dimensions=dimensions, lattice_parameter=lattice_parameter,
                                    excitons=excitons)

    total_time = 0
    path_list = []

    finished = False
    it = 0
    path = None
    while finished is False:
        get_centers(system, path)
        path, time = update_system(system, memory)

        total_time += time
        path_list.append(path)

        it += 1
        finished = check_finish(path_list)

    final_positions.append(final_position(path_list, system))
    exciton_lengths.append(moved_length(initial_position(path_list, system), final_position(path_list, system)))
    diffusion_length.append(np.average(np.array(exciton_lengths)))

    time_list.append(total_time)
    life_time.append(np.average(np.array(time_list)))


final_x_list, final_y_list = x_y_splitter(final_positions)

df_deviation = np.float(np.std(np.array(diffusion_length)))

lt_deviation = np.float(np.std(np.array(life_time)))


iterations = np.arange(0, 100, 1)

plt.plot(iterations, diffusion_length, 'ro')
plt.xlabel('Number of iterations')
plt.ylabel('Diffusion length')
plt.show()

plt.plot(iterations, life_time, 'ro')
plt.xlabel('Number of iterations')
plt.ylabel('Life_time')
plt.show()

print('Average distance: %.5f +- %.5f' % (diffusion_length[-1], df_deviation))
print('Average time %.5f +- %.5f' % (life_time[-1], lt_deviation))

