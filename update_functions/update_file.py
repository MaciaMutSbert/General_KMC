import numpy as np
from scipy.spatial import distance
from kmc_implementation.kmc_file import kmc_algorithm


def update_system(system):
    molecules = system['molecules']

    center_indexs = get_centers(molecules)

    for center in center_indexs:
        neigbours_index = neighbourhood(center, molecules)

        something = evolution(center, neigbours_index, system)

    return something


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


def evolution(center, neighbour_index, system):
    physical_conditions = system['conditions']

    transfer_rates = {}
    for i in neighbour_index:
        i_rates = get_transfer_rate(center, i, physical_conditions)
        transfer_rates[str(i)] = i_rates

    decay_rates = system['molecules'][center].state

    path, time = kmc_algorithm(decay_rates, transfer_rates)

    return path, time




