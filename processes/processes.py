import numpy as np
from scipy.spatial import distance


def get_transfer_rates(center, neighbour_index, system):
    """
    :param center: Index of the studies excited molecule
    :param neighbour_index: Index of a nearby molecule (possible acceptor of the exciton)
    :param system: Dictionary with the list of molecules and additional physical information
    :return: Dictionary with the possible transfer processes between the two involved molecules
    as keys and its rates as arguments.
    In this simple example we will consider that the rate is given by:
    rate = asb(directional_rates * u_direction)
    Where u_direction is the unit directional vector between both molecules
    """
    molecules = system['molecules']
    directional_rates = np.array([1/5, 1/7])

    process = 'Singlet transfering'

    center_vector = np.array(molecules[center].coordinates)
    neighbour_vector = np.array(molecules[neighbour_index].coordinates)
    u_direction = (center_vector - neighbour_vector) / distance.euclidean(center_vector, neighbour_vector)

    rate = np.abs(np.inner(u_direction, directional_rates))

    return {process: rate}


def update_step(choosen_process, time, system):
    """
    :param choosen_process: dictionary like dict(center, process, neighbour)
    :param time: Duration of the process
    :param system: dictionary with the molecule list and the additional information of the system
    Modifies system.
    """
    molecules = system['molecules']

    if choosen_process['process'] is 'Singlet transfering':
        molecules[choosen_process['donor']].state = 0
        molecules[choosen_process['acceptor']].stater = 1

