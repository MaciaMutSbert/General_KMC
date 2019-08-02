import numpy as np
from scipy.spatial import distance
from kmc_implementation.kmc_file import kmc_algorithm
from processes.processes import get_transfer_rates, update_step, get_decay_rates


def update_system(system):
    """
    :param system: Dictionary with all the information of the system
    Dictionary system already has the indexes of the excited molecules
    1. Looks for the neighbourhood of every centre.
    2. Chooses a process for every exciton (KMC). This path variables include: dict(centre, process, new molecule)
    3. Considering all the calculated rates computes the time interval for each process.
    4. Updates the system according to the chosen path and the time passed.
    :return: the chosen process and the advanced time
    """
    molecules = system['molecules']
    center_indexes = system['centres']

    rate_collector = []
    process_collector = []
    for center in center_indexes:

        neighbour_indexes = neighbourhood(center, molecules, radius=system['conditions']['neighbourhood_radius'])
        path_list, rate_list = get_rates_process(center, neighbour_indexes, system)

        rate_collector += rate_list
        process_collector += path_list

    chosen_process, time = kmc_algorithm(rate_collector, process_collector)
    update_step(chosen_process, molecules)

    return chosen_process, time


def check_finish(path_list):

    if path_list[-1]['donor'] == path_list[-1]['acceptor']:
        return True

    else:
        return False


def get_centres(system, path):
    """
    :param system: Dictionary with the information of the system. Only the list of the excited molecules
    indexes is used.
    :param path: dict(donor, process, acceptor) the acceptor is the new center
    Changes the list of centers when the exciton is transferred from a 'donor' to an 'acceptor'.
    """

    """
        molecules = system['molecules']
        for i, molecule in enumerate(molecules):
            if molecule.state != 0:
                center_indexes.append(i)
        """

    if path is not None:
        system['centres'].remove(path['donor'])
        system['centres'].append(path['acceptor'])


#############################################################################################
def neighbourhood(center, molecules, radius=1.05):
    """
    :param center: Index of an excited Molecule object
    :param molecules: List of objects Molecule
    :param radius: Effective distance where interaction may be considerable. Default 0.11
    :return: List of indexes of molecules in a neighbourhood of center
    If there is not any neighbours in the defined neighbourhood an alert is printed.
    """
    center_position = np.array(molecules[center].coordinates)

    neighbours = []
    for i, molecule in enumerate(molecules):
        coordinates = np.array(molecule.coordinates)

        if 0 < distance.euclidean(center_position, coordinates) < radius:
            neighbours.append(i)

    if len(neighbours) == 0:
        print('No neighbours found')

    return neighbours


def get_rates_process(centre, neighbour_indexes, system):
    """
    :param centre: Index of the studied excited molecule. Donor
    :param neighbour_indexes: Indexes of the neighbours (candidates to acceptors)
    :param system: Dictionary with the information of the system.
    Computes the transfer and decay rates and builds two dictionaries:
            One with the decay process as key and its rate as argument
            One with the transferred molecule index as key and {'process': rate} as argument

    :return:    process_list: List of elements like dict(center, process, new molecule)
                rate_list: List with the respective rates
                The list indexes coincide.
    """
    transfer_rates = {}

    for i in neighbour_indexes:
        i_rates = get_transfer_rates(centre, i, system)
        transfer_rates[str(i)] = i_rates

    decay_rates = get_decay_rates(system, centre)
    return process_rate_splitter(transfer_rates, decay_rates, centre)


def process_rate_splitter(transfer_rates, decay_rates, centre_index):
    """
    :param transfer_rates: Dictionary with the transferred molecule index as key and
    another dictionary {'process': rate} as argument
    :param decay_rates: Dictionary with the decay process as key and its rate as argument
    :param centre_index: Index of the excited molecule
    :return: Two lists:
            Process list. List with elements dict(donor(index), process, acceptor(index))
            Rate_list. List with the respective rates
            (The indexes coincide)
    """
    process_list = []
    rates_list = []

    for decay in decay_rates:
        process_list.append({'donor': centre_index, 'process': decay, 'acceptor': centre_index})
        rates_list.append(decay_rates[decay])

    for neighbour_index in transfer_rates:
        for process in transfer_rates[neighbour_index]:
            process_list.append({'donor': centre_index, 'process': process, 'acceptor': int(neighbour_index)})
            rates_list.append(transfer_rates[neighbour_index][process])

    return process_list, rates_list





