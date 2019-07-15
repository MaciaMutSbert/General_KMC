import numpy as np
from systems.molecules import Molecule
from scipy.spatial import distance
from kmc_implementation.kmc_file import kmc_algorithm
from processes.processes import get_transfer_rates, update_step


def update_system(system):
    """
    :param system: Dictionary with all the information of the system
    1. Looks for the excited molecules (center_indexes).
    2. Looks for the neighbourhood of every centre.
    2(bis). Chooses a process for every exciton. This paths includes: dict(centre, process, new molecule)
    3. Considering all the calculated rates computes the time interval for each process.
    4. Updates the system according to the chosen path and the time passed.
    :return:
    """
    molecules = system['molecules']

    center_indexes = get_centers(molecules)

    rate_collector = []
    process_collector = []
    for center in center_indexes:
        neighbours_index = neighbourhood(center, molecules, radius=1.1)

        path_list, rate_list = get_rates_process(center, neighbours_index, system)     # Include a time
        rate_collector += rate_list
        process_collector += path_list

    chosen_process, time = kmc_algorithm(rate_collector, process_collector)
    update_step(chosen_process, time, system)
    return chosen_process, time


def check_finish(path_list):

    if path_list[-1]['donor'] == path_list[-1]['acceptor']:
        return True


def get_centers(molecules):
    """
    :param molecules: List of objects Molecule
    :return: The indexes of the excited elements
    """
    centers = []
    for i, molecule in enumerate(molecules):
        if molecule.state != 0:
            centers.append(i)

    return centers


def neighbourhood(center, molecules, radius=0.11):
    """
    :param center: Index of an excited Molecule object
    :param molecules: List of objects Molecule
    :param radius: Effective distance where interaction may be considerable. Default 0.11
    :return: List of indexs of molecules in a neighbourhood of center
    """
    center_position = np.array(molecules[center].coordinates)

    neighbours = []
    for i, molecule in enumerate(molecules):
        coordinates = np.array(molecule.coordinates)

        if 0 < distance.euclidean(center_position, coordinates) <\
                radius:
            neighbours.append(i)

    return neighbours


def get_rates_process(centre, neighbour_index, system):
    """
    :param centre: Index of the studied excited molecule
    :param neighbour_index: Indexes of the neighbours (candidates to accept the exciton)
    :param system: Dictionary with the information of the system.
    Computes the transfer and decay rates and builds two dictionaries:
            One with the decay process as key and its rate as argument
            One with the transferred molecule index as key and {'process': rate} as argument

    :return:    process_list: List of elements like dict(center, process, new molecule)
                rate_list: List with the respective rates
                The list indexes coincide.
    """
    transfer_rates = {}
    for i in neighbour_index:
        i_rates = get_transfer_rates(centre, i, system)
        transfer_rates[str(i)] = i_rates

    decay_rates = system['molecules'][centre].decay_rate()
    process_list, rate_list = process_rate_splitter(transfer_rates, decay_rates, centre)

    return process_list, rate_list


def process_rate_splitter(transfer_rates, decay_rates, centre_index):
    """
    :param transfer_rates: Dictionary with the transferred molecule index as key and
    a new dictionary {'process': rate} as argument
    :param decay_rates: Dictionary with the decay process as key and its rate as argument
    :param centre_index: Index of the studied excited molecule
    :return: Two lists:
            Process list. List with elements dict(center(index), process, new molecule(index))
            Rate_list. List with the respective rates
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





