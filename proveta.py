from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system, check_finish
from result_analysis import initial_position, final_position, x_y_splitter, moved_length
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

state_energies = {'g_s': 0, 's_1': 2.5}       # eV

relaxation_energies = {'g_s': 0, 's_1': 0.33}     # eV

transition_moment = np.array([1, 0, 0])     # a.u.  Tetracene value

generic_molecule = Molecule(state_energies, relaxation_energies, transition_moment)


"Morphology parameters"
order = 'ordered'               # can take 'ordered' or 'disordered'
dimensions = [10, 10]           # molecules per side. dimensionality = len(dimensions)
lattice_parameter = 1           # nm. Molecular site in an ordered system. Not used for disordered systems
num_molecules = 0               # int. Number of molecules in the system. Not used in ordered systems

orientation = 'parallel'        # orientation between molecules
reference_orientation = np.array([1, 0, 0])         # necessary only when (anti)parallelism is required

"""
Excitons.
Possible positions commands: 'random', 'first', 'last', 'centre', 'furthest'.
The dictionary takes the state as key and a list with the position of every exciton.
"""
excitons = {'s_1': ['centre']}


# Lists for further analysis


for j in range(1):
    key = str(len(dimensions))
    system = get_homogeneous_system['ordered'][key](conditions, generic_molecule, dimensions, lattice_parameter,
                                                    orientation, reference_orientation, excitons)

    """
    system is a dictionary with three keys:
        molecules: List of objects class Molecule
        conditions: dictionary with the physical conditions of the system such as temperature, refractive index...
        centre_indexes: list with the indexes of the excited molecules.
    """

    total_time = 0
    path_list = []
    """
    Veure com definim el màxim d'iteracions.
    """
    for it in range(1):
        path, time = update_system(system)
        total_time += time
        path_list.append(path)

        print(path)

        if check_finish(path_list) is True:
            break


"""
 final_positions.append(final_position(path_list, system))
    exciton_lengths.append(moved_length(initial_position(path_list, system), final_position(path_list, system)))
    diffusion_length.append(np.average(np.array(exciton_lengths)))
    print(j)

    time_list.append(total_time)
    life_time.append(np.average(np.array(time_list)))


final_x_list, final_y_list = x_y_splitter(final_positions)

df_deviation = np.float(np.std(np.array(diffusion_length)))

lt_deviation = np.float(np.std(np.array(life_time)))


iterations = np.arange(0, 1000, 1)

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
"""