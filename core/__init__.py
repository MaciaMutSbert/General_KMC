from core.neighbourhood_and_connectivity import neighbourhood
from core.kmc_implementation import kmc_algorithm
from core.processes import get_transfer_rates, update_step, get_decay_rates


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

    molecules = system['molecules']                 # list of all instances of Molecule
    centre_indexes = system['centres']              # tricky list with the indexes of all excited molecules

    rate_collector = []                             # list with all the rates (for all centers)
    process_collector = []                          # list with the respective processes (for all centers)
    # the indexes of both lists coincide.

    for centre in centre_indexes:
        neighbour_indexes = neighbourhood(centre, system, radius=system['conditions']['neighbourhood_radius'])
        # looks for the all molecules in a circle of radius centered at the position of the excited molecule

        process_list, rate_list = get_processes_and_rates(centre, neighbour_indexes, system)
        # for each center computes all the decay rates and all the transfer rates for all neighbours
        # return them as a list

        rate_collector += rate_list
        process_collector += process_list
        # merging of the new list with the rates and processes previously computed

    chosen_process, time = kmc_algorithm(rate_collector, process_collector)
    # chooses one of the processes using the Kinetic Monte-Carlo algorithm

    update_step(chosen_process, molecules, centre_indexes)        # updates both lists according to the chosen process

    # finally the chosen process and the advanced time are returned
    return chosen_process, time


#############################################################################################


def check_finish(system):
    """
    :param system:
    :return: Checks if the list of excited molecules is empty. When empty the excitations have all decayed and
    the simulation will finish.
    """

    if len(system['centres']) == 0:
        return True

    else:
        return False


##############################################################################################################
#                                   ASSISTANT FUNCTIONS
##############################################################################################################


def get_processes_and_rates(centre, neighbour_indexes, system):
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

    transfer_processes, transfer_rates = get_transfer_rates(centre, neighbour_indexes, system)
    # calls an external function that computes the transfer rates for all possible transfer processes between
    # the centre and all its neighbours

    decay_processes, decay_rates = get_decay_rates(centre, system)
    # calls an external function that computes the decay rates for all possible decay processes of the centre.
    # Uses a method of the class Molecule

    # merges all processes in a list and equally for the rates
    # the indexes of the list must coincide (same length. rate 'i' is related to process 'i')
    process_list = decay_processes + transfer_processes
    rate_list = decay_rates + transfer_rates

    return process_list, rate_list

"""
Things to be changed:
    process and rates must be splitted when calling get_transfer/decay_rates
    in get_decay_rates the order of the arguments has to be changed (centre, system)
    process_rate_splitter won't be needed anymore.
"""


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





