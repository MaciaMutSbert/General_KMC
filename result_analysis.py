import numpy as np
from scipy.spatial import distance


def get_trajectory(path_collector, time_advance, system_):
    molecules = system_['molecules']
    positions = []
    states = []

    positions.append(molecules[path_collector[0]['donor']].molecular_coordinates())

    for path in path_collector:
        positions.append(molecules[path['acceptor']].molecular_coordinates())
        states.append(molecules[path['acceptor']].state)

    return {'positions': positions, 'time_advance': time_advance, 'state': states}


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
