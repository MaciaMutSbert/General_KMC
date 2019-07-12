import numpy as np
from scipy.spatial import distance
from kmc_implementation.kmc_file import kmc_algorithm


def update_system(system):
    """
    :param system: Dictionary with all the information of the system
    1. Looks for the excited molecules (center_indexes).
    2. Looks for the neighbourhood of every molecule.
    2(bis). Chooses a path for every exciton
    3. Considering all the calculated rates computes the time interval for each process.
    4. (Merge plans)
    :return:
    """
    molecules = system['molecules']

    center_indexes = get_centers(molecules)

    rate_collector = []
    path_collector = []
    for center in center_indexes:
        neighbours_index = neighbourhood(center, molecules)

        path_list, rate_list = get_rates_paths(center, neighbours_index, system)
        rate_collector += rate_list
        path_collector += path_list

    chosen_path, time = kmc_algorithm(rate_collector, path_collector)
    plan = update_step(chosen_path)

    return plan, time



def get_centers(molecules):
    """
    :param molecules: List of objects Molecule
    :return: The indexs of the excited elements
    """
    centers = []
    for i, molecule in enumerate(molecules):
        if molecule.state != 0:
            centers.append(i)

    return centers


def neighbourhood(center, molecules, radius = 0.2):
    """
    :param center: Index of an excited Molecule object
    :param molecules: List of objects Molecule
    :param radius: Effective distance where interaction may be considerable. Default 0.2
    :return: List of indexs of molecules in a neighbourhood of center
    """
    center_position = np.array(molecules[center].coordinates)

    neighbours = []
    for i, molecule in enumerate(molecules):
        coordinates = np.array(molecule.coordinates)

        if distance.euclidean(center_position, coordinates) <\
                radius:
            neighbours.append(i)

    return neighbours


def get_rates_paths(center, neighbour_index, system):
    """
    :param center: Index of the studied excited molecule
    :param neighbour_index: Indexes of the neighbours (candidates to accept the exciton)
    :param system: Dictionary with the information of the system.
    Computes the transfer and decay rates and builts two dictionaries:
            One with the decay process as key and its rate as argument
            One with the transferred molecule index as key and the list(process, rate) as argument

    :return:    path_list: List of elements like list(index of the new excited molecule, process)
                rate_list: List with the respective rates
                The list indexes coincide.
    """
    physical_conditions = system['conditions']

    transfer_rates = {}
    for i in neighbour_index:
        i_rates = get_transfer_rate(center, i, physical_conditions)
        transfer_rates[str(i)] = i_rates

    decay_rates = system['molecules'][center].decay_rate

    path_list, rate_list = path_rate_splitter(transfer_rates, decay_rates, center)

    return path_list, rate_list


def path_rate_splitter(transfer_rates, decay_rates, center_index):
    """
    :param transfer_rates: Dictionary with the transferred molecule index as key and the list(process, rate) as argument
    :param decay_rates: Dictionary with the decay process as key and its rate as argument
    :return: Two lists:
            Process list. List with elements (affected molecule, process)
            Rate_list. List with the respective rates
    """
    process_list = []
    rates_list = []

    for decay in decay_rates:
        process_list.append([center_index, decay])
        rates_list.append(decay_rates[decay])

    for neighbour_index in transfer_rates:
        for process in transfer_rates[neighbour_index]:
            process_list.append([int(neighbour_index), process])
            rates_list.append(transfer_rates[neighbour_index][process])

    return process_list, rates_list





