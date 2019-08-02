from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system, check_finish, get_centres
from result_analysis import initial_position, final_position, x_y_splitter, moved_length
import matplotlib.pyplot as plt
import numpy as np
from systems.molecules import Molecule

"""
ORGANITZAR-HO EN FUNCIONS ESPECÍFIQUES
"""


"Physical conditions"
conditions = {'temperature': input('Temperature in K:'),
              'refractive_index': input('Refractive index of the material:'),
              'neighbour_radius': input('Maximum distance considered for the interaction:')}


"Generic molecule initialization"
print("""\n Possible states: 
         \t 'g_s': ground state 
         \t 's_1': first singlet state 
         \t 't_1': first triplet state 
         \n All energies must be given in eV.""")

state_energies = {'g_s': float(input('Ground state energy (eV): ')),
                  's_1': float(input('First singlet state energy (eV): ')),
                  't_1': float(input('First triplet state energy (eV): '))}
print('\n')
relaxation_energies = {'g_s': float(input('Ground state relaxation energy (eV): ')),
                       's_1': float(input('First singlet state relaxation energy (eV): ')),
                       't_1': float(input('First triplet state relaxation energy (eV): '))}

print("\n Dipole transition moment vector (a.u)")
transition_moment = np.array([float(input('x_component:')),
                              float(input('y_component: ')),
                              float(input('z_component: '))])

generic_molecule = Molecule(state_energies, relaxation_energies, transition_moment)
print('\n')

"Morphology parameters"
order = input("'ordered' / 'disordered': ")
print(order)
print(type(order))
print("\n Number of molecules per side: ")
dimensions = [int(input('Molecules in x: ')),
              int(input('Molecules in y: ')),
              int(input('Molecules in z: '))]
for dimension in dimensions:
    if dimension == 0:
        dimensions.remove(dimension)
print(dimensions)

if order is 'ordered':
    lattice_parameter = float(input("Lattice parameter (nm), only needed for 'ordered': "))
    num_molecules = None
if order is 'disordered':
    num_molecules = int(input("Number of molecules (only for 'disordered'): "))
    lattice_parameter = None

orientation = input("'random'/'parallel'/'antiparallel': ")
reference_orientation = np.array([1, 0, 0])
print('\n')

"Excitons"
print("""Possible positions: random, first, last, centre and furthest \n""")
excitons = {}
num_s1 = int(input('Number of s_1 excitons:'))
s_1_positions = []
for i in range(num_s1):
    s_1_positions.append(input('s_1 position:'))
excitons['s_1'] = s_1_positions

num_t1 = int(input('Number of t_1 excitons'))
t_1_positions = []
for i in range(num_t1):
    t_1_positions.append(input('t_1 position: '))
excitons['t_1'] = t_1_positions
print(excitons)




# Lists for further analysis



for j in range(1000):
    system = get_homogeneous_system(conditions, dimensions=dimensions, lattice_parameter=lattice_parameter,
                                    excitons=excitons)

    """
    system is a dictionary with three keys:
        molecules: List of objects class Molecule
        conditions: dictionary with the physical conditions of the system such as temperature, refractive index...
        centre_indexes: list with the indexes of the excited molecules.
    """
    total_time = 0
    path_list = []

    path = None
    """
    Veure com definim el màxim d'iteracions.
    """
    for it in range(max):
        get_centres(system, path)
        path, time = update_system(system)

        total_time += time
        path_list.append(path)

        if check_finish(path_list) is True:
            break

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
