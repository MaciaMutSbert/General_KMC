import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt


def statistical_diffusivity(trajectories, trajenctories_steps, dimensionality):
    """
    :param trajectories: list with all the trajectories, each is a dictionary with a list with the positions,
    time to reach each position and the exciton transfer process.
    :param hops: Number of steps of the trajectory (len of the lists in trajectory), the same for all of them
    :param dimensionality
    :return: the diffusivity constant computed using all the points of all the trajectories.
        For each step we computed the mean square of the distance and the mean time.
        Doing the lineal regression of all these couple of values, the diffusivity constant is defined as the
        slope of the line (except for a dimensional parameter depending of the dimensionality.
    """

    mean_square_distance_list = []
    mean_lifetime_list = []
    set_of_points = []

    steps = 250               # prenem un número fixe de steps per fer l'estudi estadístic (triat per tenir el resultat òptim).
    print(steps)
    for i in range(1, steps):
        squared_distances = []
        time = []
        for trajectory in trajectories:
            if i >= len(trajectory['positions']):
                count = 0
            else:
                distance_i_squared = np.linalg.norm(np.array(trajectory['positions'][i]))**2
                time_i = trajectory['time_advance'][i]

                squared_distances.append(distance_i_squared)
                time.append(time_i)

        mean_square_distance_list.append(np.average(np.array(squared_distances)))
        mean_lifetime_list.append(np.average(np.array(time)))

        set_of_points.append([np.average(squared_distances), np.average(time)])

    diffusion_constant_list = np.array(mean_square_distance_list) / (np.array(mean_lifetime_list) *2*dimensionality)

    return {'diffusion_constants': diffusion_constant_list, 'mean_square_distances':mean_square_distance_list,
            'life_times': mean_lifetime_list}


def diffusion_parameters(trajectories, dimensionality):
    """
    :param trajectories: list with all the trajectories, each is a dictionary with a list with the positions,
    time to reach each position and the exciton transfer process.
    :param dimensionality
    :return: the diffusion length computed as the root mean square of the final positions of the trajectory
             the exciton life time computed as the average of the total time advance
             the diffusion constant computed with the two parameters above
    """
    squared_distances = []
    lifetimes = []
    for trajectory in trajectories:
        final_distance_squared = np.linalg.norm(np.array(trajectory['positions'][-1])) ** 2
        final_time = trajectory['time_advance'][-1]

        squared_distances.append(final_distance_squared)
        lifetimes.append(final_time)

    "Histograma amb les posicions finals"
    plt.hist(np.sqrt(squared_distances), bins=15)
    plt.title('Final positions histogram')
    plt.show()

    return {'diffusion_length': np.sqrt(np.average(squared_distances)),
            'exciton_lifetime': np.average(lifetimes),
            'diffusion_constant': np.average(squared_distances) / (2 * dimensionality * np.average(lifetimes))}


def statistics(sample_set):
    average = np.average(np.array(sample_set))
    deviation = np.std(np.array(sample_set))

    return average, deviation


###############################################################################################################

#               NOT USED FUNCTIONS


def moved_length(initial_position, final_position):
    return distance.euclidean(np.array(initial_position), np.array(final_position))


def final_position(path_list, system):
    """
    :param path_list: List indicating the path of the exciton with dict(donor, process, acceptor)
    :param system: includes the list of molecules
    :return: the coordinates of the last acceptor
    """
    final_index = path_list[-1]['acceptor']

    return system['molecules'][final_index].coordinates


def initial_position(path_list, system):
    """
    :param path_list: List indicating the path of the exciton with dict(donor, process, acceptor)
    :param system: includes the list of molecules
    :return: the coordinates of the first donor
    """
    initial_index = path_list[0]['donor']

    return system['molecules'][initial_index].coordinates


def exciton_shift(initial_positions, final_positions):
    """
    :param initial_positions:
    :param final_positions:
    :return: List of the exciton covered distances.
    """
    shift_list = []

    for i, initial in enumerate(initial_positions):
        shift = distance.euclidean(np.array(initial), np.array(final_positions[i]))
        shift_list.append(shift)

    return shift_list


def x_y_splitter(final_positions):
    """
    :param final_positions:
    :return: List with final x positions and analogously for y.
    """
    x_list = []
    y_list = []

    for final in final_positions:
        x_list.append(final[0])
        y_list.append(final[1])

    return x_list, y_list
