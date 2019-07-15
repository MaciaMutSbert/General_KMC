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

    process = 'Singlet_transfer'

    k = 2      # orientational factor. Taken as a constant in an ideal system

    alfa = 1.15      # correcting factor for short distances
    hopping_distance = distance.euclidean(np.array(molecules[center].coordinates), np.array(molecules[neighbour_index].coordinates))
    corrected_distance = alfa*molecules[int(center)].transition_dipole + hopping_distance

    n = system['conditions']['refractive_index']
    sigma = 0.3     # deviation of the absorption and emission spectra considered gaussian

    factor_1 = k**2 * np.pi ** 2 * molecules[int(center)].transition_dipole**4
    factor_2 = n**4 * corrected_distance**6 * sigma

    rate = factor_1 / factor_2
    return {process: rate}


def update_step(chosen_process, time, system):
    """
    :param chosen_process: dictionary like dict(center, process, neighbour)
    :param time: Duration of the process
    :param system: dictionary with the molecule list and the additional information of the system
    Modifies system.
    """

    if chosen_process['process'] is 'Singlet_transfer':
        system['molecules'][chosen_process['donor']].state = 0
        system['molecules'][chosen_process['acceptor']].state = 1

    if chosen_process['process'] is 'Singlet_radiative_decay_rate':
        system['molecules'][chosen_process['donor']].state = 0






