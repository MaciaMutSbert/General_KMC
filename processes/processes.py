import numpy as np
from scipy.spatial import distance
from numpy import pi


def get_transfer_rates(center, neighbour_index, system, memory):
    """
    :param center: Index of the studies excited molecule
    :param neighbour_index: Index of a nearby molecule (possible acceptor of the exciton)
    :param system: Dictionary with the list of molecules and additional physical information
    :param memory: dictionary with the calculated rates as arguments, the hash of the characteristic
    parameters is used as key.
    :return: Dictionary with the possible transfer processes between the two involved molecules
    as keys and its rates as arguments.
    """
    molecules = system['molecules']
    transfer_rates = {}

    fcwd = compute_fcwd_gaussian(system, memory)

    if molecules[center].state == 1:
        if molecules[center].state == 0:

            process = 'Singlet_transfer'
            forster = compute_forster_coupling(center, neighbour_index, system, memory)
            transfer_rates[process] = fcwd * forster**2

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


def get_decay_rates(system, centre, memory):
    """
    :param system: Dictionary with all the information of the system
    :param centre: index of the excited molecule
    :param memory: dictionary with the calculated rates as arguments, the hash of the characteristic
    parameters is used as key.
    :return: A dictionary with the possible decay rates
    """

    molecule = system['molecules'][centre]

    molecular_state = molecule.state

    particular_info = str(hash(molecular_state))

    if particular_info in memory:
        decay_rates = memory[particular_info]

    else:
        decay_rates = molecule.decay_rate()
        memory[particular_info] = decay_rates

    return decay_rates


def compute_fcwd_gaussian(system, memory):
    """
    :param system: dictionary with all the physical information of the system
    :param memory: dictionary with the calculated rates as arguments, the hash of the characteristic
    parameters is used as key.
    :return: Franck-Condon-weighted density of states in gaussian aproximation
    """
    delta = system['conditions']['a_e_centre_shift'] / 27.211       # atomic units
    sigma = system['conditions']['a_e_deviation'] / 27.211          # atomic units

    info = str(hash((delta, sigma)))

    if info in memory:
        fcwd = memory[info]

    else:
        fcwd = np.sqrt(pi) / (2 * sigma) * np.exp(- delta**2 / (2*sigma)**2)
        memory[info] = fcwd

    return fcwd


def compute_forster_coupling(center, neighbour_index, system, memory):
    """
    :param center: Index of the excited molecule
    :param neighbour_index: Index of a possible acceptor
    :param system: Dictionary with all the physical information
    :param memory: dictionary with the calculated rates as arguments, the hash of the characteristic
    parameters is used as key.
    :return: Forster coupling between both molecules. We do not implement
    any correction for short distances.
    """
    molecules = system['molecules']

    u_a = molecules[center].transition_dipole
    u_b = molecules[neighbour_index].transition_dipole
    k = 2     # we shall define a function with the right expression for k when it is not constant

    r_a = np.array(molecules[center].coordinates)
    r_b = np.array(molecules[center].coordinates)
    inter_distance = distance.euclidean(r_a, r_b)

    particular_info = str(hash((u_a, u_b, inter_distance)))

    if particular_info in memory:
        forster_coupling = memory[particular_info]

    else:
        n = system['conditions']['refractive_index']
        forster_coupling = k * np.dot(u_a, u_b) / (n**2 * inter_distance**3)
        memory[particular_info] = forster_coupling

    return forster_coupling
