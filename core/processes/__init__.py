import numpy as np
from numpy import pi
from core.processes.coupling_functions import couplings
from conversion_functions import from_ns_to_au, from_ev_to_au

decay_memory = {}
overlap_memory = {}

###########################################################################
#           FUNCIÓ 1: rates de transfrència
###########################################################################


def get_transfer_rates(centre, neighbour_index, system):
    """
    :param centre: Index of the studies excited molecule
    :param neighbour_index: Index of a nearby molecule (possible acceptor of the exciton)
    :param system: Dictionary with the list of molecules and additional physical information
    :return: Dictionary with the possible transfer processes between the two involved molecules
    as keys and its rates as arguments.
    """
    conditions = system['conditions']
    donor = system['molecules'][centre]
    acceptor = system['molecules'][neighbour_index]

    transfer_rates = {}

    spectral_overlap = marcus_overlap_formula(donor, acceptor, conditions)
    """
    Les referències dels possibles processos de transfarència seran les que s'usaran com a claus dels diccionaris:
    transfer_rates i couplings. Així com les que s'usaran per identificar com s'ha d'actualitzar el sistema. És per això
    que s'ha de seguir un conveni. La clau es construirà com: 'state1_state2' i farà referència a una excitació en state1
    que es propaga a una molècula en state2. E.g, si tenim la propagació d'un singlet_1 a una molècula en l'estat 
    fonamental la clau del procés serà: 's_1_g_s'.
    """

    key = '{}_{}'.format(donor.electronic_state(), acceptor.electronic_state())
    electronic_coupling = couplings[key](donor, acceptor, conditions)

    rate = (2*pi) * spectral_overlap * electronic_coupling**2
    transfer_rates[key] = from_ns_to_au(rate, 'direct')

    return transfer_rates


###########################################################################
#           FUNCIÓ 2: rates de decaïment
###########################################################################


def get_decay_rates(system, centre):
    """
    :param system: Dictionary with all the information of the system
    :param centre: index of the excited molecule
    :return: A dictionary with the possible decay rates
    """
    donor = system['molecules'][centre]
    info = str(hash(donor.state))

    if info in decay_memory:
        decay_rates = decay_memory[info]

    else:
        decay_rates = donor.decay_rates()
        decay_memory[info] = decay_rates

    return decay_rates


###########################################################################
#           FUNCIÓ 3: actualització del sistema
###########################################################################


def update_step(chosen_process, molecules, centre_indexes):
    """
    :param chosen_process: dictionary like dict(center, process, neighbour)
    :param molecules: list of instances of molecule
    :param centre_indexes: list of the indexes of the excited molecules
    Modifies the state of the donor and the acceptor. Removes the donor from the centre_indexes list
    and includes the acceptor. If its a decay only changes and removes tha acceptor
    """
    if chosen_process['process'] == 's_1_g_s':
        molecules[chosen_process['donor']].set_state('g_s')
        molecules[chosen_process['acceptor']].set_state('s_1')
        centre_indexes.remove(chosen_process['donor'])
        centre_indexes.append(chosen_process['acceptor'])

    if chosen_process['process'] == 'Singlet_radiative_decay':
        molecules[chosen_process['donor']].set_state('g_s')
        centre_indexes.remove(chosen_process['donor'])


###########################################################################
#           FUNCIONS AUXILIARS: solapament espectral
###########################################################################


def marcus_overlap_formula(donor, acceptor, conditions):
    """
    :param donor:
    :param acceptor:
    :param conditions:
    :return: The spectral overlap between the donor and the acceptor according to Marcus formula.
    """
    kb = 8.617333 * 10**(-5)            # Boltzmann constant in eV * K^(-1)
    T = conditions['temperature']       # K

    excited_state = donor.electronic_state()
    gibbs_energy = donor.state_energies[excited_state] - acceptor.state_energies[excited_state]
    relax = donor.get_relaxation_state_energy()

    info = str(hash((T, gibbs_energy, relax, 'marcus')))

    if info in overlap_memory:
        overlap = overlap_memory[info]

    else:
        overlap = 1 / (2 * np.sqrt(pi*kb*T*relax)) * np.exp(-(gibbs_energy+relax)**2 / (4*kb*T*relax))
        overlap_memory[info] = overlap

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



