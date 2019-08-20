import numpy as np
from numpy import pi
from core.processes.coupling_functions import couplings
from conversion_functions import from_ns_to_au, from_ev_to_au


# Memory for the calculated decay rates and spectral overlaps is introduced.
decay_memory = {}
overlap_memory = {}


###########################################################################################################
#                                   FUNCTION 1: TRANSFER RATES
###########################################################################################################


def get_transfer_rates(centre, neighbour_indexes, system):
    """
    :param centre: Index of the studies excited molecule
    :param neighbour_indexes: neighbour indexes list
    :param system: Dictionary with the list of molecules and additional physical information
    :return: Dictionary with the possible transfer processes between the two involved molecules
    as keys and its rates as arguments.
    For each possible acceptor in neighbour_indexes computes the transfer rate using the Fermi's Golden Rule:
        - For the spectral overlap the Marcus Formula is used for all cases.
        - For the electronic coupling an external dictionary is defined. It contains the possible couplings between
            two states (more than one allowed). The keys of this dictionary are like:
                'state1_state2' + str(additional information)
            If the key 'state1_state2' is not in the dictionary the electronic coupling shall be taken as 0.
    """

    conditions = system['conditions']           # physical conditions of the system

    donor = system['molecules'][centre]         # excited molecule

    transfer_rates = []                         # list that collects the transfer rates (only the numerical values)
    transfer_processes = []                     # list that collects the transfer processes dict(donor,process,acceptor)

    for neighbour in neighbour_indexes:
        acceptor = system['molecules'][neighbour]
        # acceptor = molecule instance for each neighbour index.

        spectral_overlap = marcus_overlap_formula(donor, acceptor, conditions)
        # compute the spectral overlap with the Marcus Formula (for all possible couplings)

        prefix_key = '{}_{}'.format(donor.electronic_state(), acceptor.electronic_state())
        # first key. Only with the information of both implicated states: 'state1_state2'

        possible_couplings = slice_dict(couplings, prefix_key)
        # new dictionary with all the keys of couplings that start with prefix_key

        if len(possible_couplings) == 0:
            # when no couplings are found, the rate is taken as 0
            rate = 0
            transfer_rates.append(rate)

            process = {'donor': centre, 'process': prefix_key, 'acceptor': neighbour}
            transfer_processes.append(process)

        else:
            # for the possible couplings computes the transfer rate by the Fermi's Golden Rule.
            for key in possible_couplings:
                e_coupling = possible_couplings[key](donor, acceptor, conditions)

                rate = 2*pi * e_coupling**2 * spectral_overlap          # rate in a.u -- Fermi's Golden Rule
                transfer_rates.append(from_ns_to_au(rate, 'direct'))    # rate in ns⁻¹

                transfer_processes.append({'donor': int(centre), 'process': key, 'acceptor': int(neighbour)})

    return transfer_processes, transfer_rates


###########################################################################################################
#                               FUNCTION 2: DECAY RATES
###########################################################################################################


def get_decay_rates(centre, system):
    """
    :param centre: index of the excited molecule
    :param system: Dictionary with all the information of the system
    :return: A dictionary with the possible decay rates
    For computing them the method get_decay_rates of class molecule is call.
    """
    donor = system['molecules'][centre]

    info = str(hash(donor.state))
    # we define a compact string with the characteristic information of the decays: electronic state

    decay_processes = []
    decay_rates = []

    if info in decay_memory:
        decay_complete = decay_memory[info]
        # the decay memory defined is used if the decay have been already computed

    else:
        decay_complete = donor.decay_rates()
        decay_memory[info] = decay_complete

    for key in decay_complete:
        decay_processes.append({'donor': centre, 'process': key, 'acceptor': centre})
        decay_rates.append(decay_complete[key])

    return decay_processes, decay_rates


###########################################################################################################
#                           FUNCTION 3: UPDATE OF THE SYSTEM AND CENTRE INDEXES
###########################################################################################################


def update_step(chosen_process, molecules, centre_indexes):
    """
    :param chosen_process: dictionary like dict(center, process, neighbour)
    :param molecules: list of instances of molecule
    :param centre_indexes: list of the indexes of the excited molecules
    Modifies the state of the donor and the acceptor. Removes the donor from the centre_indexes list
    and includes the acceptor. If its a decay only changes and removes the donor

    New if(s) entrances shall be defined for more processes.
    """
    if chosen_process['process'] is 's1_gs':
        molecules[chosen_process['donor']].set_state('gs')          # des excitation of the donor
        molecules[chosen_process['acceptor']].set_state('s1')       # excitation of the acceptor

        centre_indexes.remove(chosen_process['donor'])
        centre_indexes.append(chosen_process['acceptor'])
        # modification of the excited molecules indexes list

    if chosen_process['process'] == 'Singlet_radiative_decay':
        molecules[chosen_process['donor']].set_state('gs')          # des excitation of the donor

        centre_indexes.remove(chosen_process['acceptor'])
        # modification of the excited molecules indexes list


    # No return function. Updates molecules and centre_indexes.


###########################################################################################################
#                           ASSISTANT FUNCTIONS: spectral overlaps
###########################################################################################################


def marcus_overlap_formula(donor, acceptor, conditions):
    """
    :param donor:
    :param acceptor:
    :param conditions:
    :return: The spectral overlap between the donor and the acceptor according to Marcus formula.
    """
    kb = 8.617333 * 10**(-5)            # Boltzmann constant in eV * K^(-1)
    T = conditions['temperature']       # temperature (K)

    excited_state = donor.electronic_state()
    gibbs_energy = donor.state_energies[excited_state] - acceptor.state_energies[excited_state]
    # Gibbs energy: energy difference between the equilibrium points of the excited states

    reorganization = donor.get_reorganization_state_energy()
    # donor reorganization energy of the excited state

    info = str(hash((T, gibbs_energy, reorganization, 'marcus')))
    # we define a compact string with the characteristic information of the spectral overlap

    if info in overlap_memory:
        # the memory is used if the overlap has been already computed
        overlap = overlap_memory[info]

    else:
        overlap = 1 / (2 * np.sqrt(pi*kb*T*reorganization)) * \
                  np.exp(-(gibbs_energy+reorganization)**2 / (4*kb*T*reorganization))

        overlap_memory[info] = overlap
        # new values are added to the memory

    return from_ev_to_au(overlap, 'inverse')
    # Since we have a quantity in 1/eV, we use the converse function from_ev_to_au in inverse mode
    # to have a 1/au quantity.


def compute_fcwd_gaussian(system):
    """
    :param system: dictionary with all the physical information of the system
    parameters is used as key.
    :return: Franck-Condon-weighted density of states in gaussian aproximation
    """
    delta = system['conditions']['a_e_spectra_centre_shift'] / 27.211       # atomic units
    sigma = system['conditions']['a_e_spectra_deviation'] / 27.211          # atomic units

    info = str(hash((delta, sigma)))

    if info in overlap_memory:
        fcwd = overlap_memory[info]

    else:
        fcwd = np.exp(- delta**2 / (2*sigma)**2) / (2 * np.sqrt(pi) * sigma)
        overlap_memory[info] = fcwd

    print(fcwd)
    return fcwd


###########################################################################################################
#          Take all the keys of a dictionary that start with the same string
###########################################################################################################


def slice_dict(original_dict, key_string):
    """
    :param original_dict: Original (and larger) dictionary with all entrances
    :param key_string: key string (desired starting string)
    :return: A new dictionary with all the entrances of the original that start with key_string
    """
    newdict = {}

    for key in original_dict:
        if key.startswith(key_string):
            newdict[key] = original_dict[key]

    return newdict


