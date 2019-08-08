import numpy as np
from scipy.spatial import distance


def get_trajectory(path_collector, time_advance, system_):
    molecules = system_['molecules']
    positions = []
    states = []

    positions.append(list(molecules[path_collector[0]['donor']].molecular_coordinates()))

    for path in path_collector:
        positions.append(list(molecules[path['acceptor']].molecular_coordinates()))
        states.append(path['process'])

    return {'positions': positions, 'time_advance': time_advance, 'state': states}


def stadistical_diffusivity(trajectories, hops, dimensionality):
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
    for i in range(1, hops):
        squared_distances = []
        time = []
        for trajectory in trajectories:
            distance_i_squared = np.linalg.norm(np.array(trajectory['positions'][i]))**2
            time_i = trajectory['time_advance'][i]

            squared_distances.append(distance_i_squared)
            time.append(time_i)

        mean_square_distance_list.append(np.average(np.array(squared_distances)))
        mean_lifetime_list.append(np.average(np.array(time)))

        set_of_points.append([np.average(squared_distances), np.average(time)])

    diffusion_constant_list = np.array(mean_square_distance_list) / (np.array(mean_lifetime_list) *2*dimensionality)
#    statistical_diffusion_constant = LinearRegression().fit(set_of_points, y)

    return diffusion_constant_list


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

    return {'diffusion_length': np.sqrt(np.average(squared_distances)),
            'exciton_lifetime': np.average(lifetimes),
            'diffusion_constant': np.average(squared_distances) / (2 * dimensionality * np.average(lifetimes))}


def statistics(sample_set):
    average = np.average(np.array(sample_set))
    deviation = np.std(np.array(sample_set))

    return average, deviation


"""
----------------------------------------------------------------------------------
                            FUNCIONS QUE JA NO USAM
----------------------------------------------------------------------------------
"""


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
