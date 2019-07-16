import numpy as np
from scipy.spatial import distance
from numpy import pi


def get_transfer_rates(center, neighbour_index, system):
    """
    :param center: Index of the studies excited molecule
    :param neighbour_index: Index of a nearby molecule (possible acceptor of the exciton)
    :param system: Dictionary with the list of molecules and additional physical information
    :return: Dictionary with the possible transfer processes between the two involved molecules
    as keys and its rates as arguments.
    """
    molecules = system['molecules']
    transfer_rates = {}

    if molecules[center].type == 1:
        if molecules[neighbour_index].type == 1:
            k = system['conditions']['orientational_factor_1']

            if molecules[center].state == 1:
                if molecules[neighbour_index].state == 0:
                    process = 'Singlet_transfer'

                    n = system['conditions']['refractive_index']
                    alfa = 1.15                                                          # short distances correction
                    sigma = system['conditions']['a_e_spectra_deviation'] / 27.211       # atomic units

                    u = system['conditions']['transition_dipole']
                    center_position = np.array(molecules[center].coordinates)
                    neighbour_position = np.array(molecules[neighbour_index].coordinates)
                    inter_distance = distance.euclidean(center_position, neighbour_position)/0.053    # atomic units

                    factor_1 = k**2 * pi**(3/2) * u**4
                    factor_2 = n**4 * (alfa*u + inter_distance)**6 * sigma

                    rate = factor_1 / factor_2                              # atomic units
                    transfer_rates[process] = rate / (2.4189*10**-8)      # in ns⁻¹

    return transfer_rates


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

    if chosen_process['process'] is 'Singlet_radiative_decay':
        system['molecules'][chosen_process['donor']].state = 0






