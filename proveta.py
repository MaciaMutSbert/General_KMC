from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system, check_finish
from result_analysis import initial_position, final_position, x_y_splitter, moved_length
from diffusion_lengths import final_position_distance, diffusion_length_time_residence
import matplotlib.pyplot as plt
import numpy as np
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

state_energies = {'g_s': 0, 's_1': 2.5}       # eV  Tetracene value

relaxation_energies = {'g_s': 0, 's_1': 0.7}     # eV  Tetracene value

transition_moment = np.array([1, 0, 0])     # a.u.  Tetracene value

generic_molecule = Molecule(state_energies, relaxation_energies, transition_moment)


"Morphology parameters"
order = 'ordered'               # can take 'ordered' or 'disordered'
dimensions = [30, 30]          # molecules per side. dimensionality = len(dimensions)
lattice_parameter = 1           # nm. Molecular site in an ordered system. Not used for disordered systems
# num_molecules               # int. Number of molecules in the system. Not used in ordered systems

orientation = 'parallel'        # orientation between molecules
reference_orientation = np.array([1, 0, 0])         # necessary only when (anti)parallelism is required

"""
Excitons.
Possible positions commands: 'random', 'first', 'last', 'centre', 'furthest'.
The dictionary takes the state as key and a list with the position of every exciton.
"""
excitons = {'s_1': ['centre']}


"""
    ----------------------------    PROVA 1D    ------------------------------

    Estudi de la longitud de difusió.
    La calculam com la distància entre la posició final i inicial
                a partir del temps de residència
"""

final_positions = []
exciton_shift = []
diffusion_length = []

ld_residence_time = []
ld_residence_time_average = []

time_list = []
life_time = []


for j in range(100):
    key = str(len(dimensions))
    system = get_homogeneous_system['ordered'][key](conditions, generic_molecule, dimensions, lattice_parameter,
                                                    orientation, reference_orientation, excitons)


    """
    system is a dictionary with three keys:
        molecules: List of objects class Molecule
        conditions: dictionary with the physical conditions of the system such as temperature, refractive index...
        centre_indexes: list with the indexes of the excited molecules.

    !!!!!!!!!!!!!  Veure com definim el màxim d'iteracions.  !!!!!!!!!!!!!!!!
    """

    total_time = 0
    path_list = []
    residence_time = []
    finished = False
    while finished is False:
        path, time = update_system(system)
        total_time += time
        path_list.append(path)
        residence_time.append(time)

        if check_finish(path_list) is True:
            break

        if total_time == 600:
            break

    time_list.append(total_time)
    life_time.append(np.average(np.array(time_list)))

    final_positions.append(final_position(path_list, system))
    exciton_shift.append(moved_length(initial_position(path_list, system), final_position(path_list, system)))
    diffusion_length.append(np.average(np.array(exciton_shift)))

    ld_residence_time.append(diffusion_length_time_residence(residence_time, lattice_parameter))
    ld_residence_time_average.append(np.average(np.array(ld_residence_time)))

    print(j)

ld_teo_1 = 13.997           # nm
ld_teo_2 = ld_teo_1 * np.sqrt(4)           # nm
ld_teo_3 = ld_teo_1 / np.sqrt(4)         # nm

df_deviation = np.float(np.std(diffusion_length))
df_residence_time_deviation = np.float(np.std(ld_residence_time_average))
lt_deviation = np.float(np.std(np.array(life_time)))

print('Valors teòrics: %.5f, %.5f, %.5f' % ld_teo_1, ld_teo_2, ld_teo_3)


print('Average distance: %.5f +- %.5f' % (diffusion_length[-1], df_deviation))
print('Average distance (computed with the residence time): %.5f +- %.5f' % (ld_residence_time_average[-1],
                                                                             df_residence_time_deviation))

final_x_list, final_y_list = x_y_splitter(final_positions)
print('Average time %.5f +- %.5f' % (life_time[-1], lt_deviation))
plt.hist2d(final_x_list, final_y_list, bins=20)
plt.show()


plt.plot(final_x_list, final_y_list, 'bo')
plt.xlabel('final x positions')
plt.ylabel('final y positions')
plt.title('Final positions')
plt.show()

iterations = np.arange(0, 100, 1)

plt.plot(iterations, diffusion_length, 'ro')
plt.xlabel('Number of iterations')
plt.ylabel('Diffusion length')
plt.show()


plt.hist(final_positions, bins=20)
plt.title('histograma amb les distàncies finals')
plt.show()

plt.plot(iterations, ld_residence_time_average, 'ro')
plt.xlabel('Number of iterations')
plt.ylabel('Diffusion length. Ld_teo = 24,7274')
plt.show()

plt.plot(iterations, life_time, 'ro')
plt.xlabel('Number of iterations, Ld_teo = 24,7274')
plt.ylabel('Life_time')
plt.show()




