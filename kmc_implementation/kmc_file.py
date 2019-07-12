import numpy as np
import random as rd


def kmc_algorithm(decay_rates, transfer_rates, center_index):
    """
    :param decay_rates: Dictionary with decay process as keys and decay rates as arguments
    :param transfer_rates: Dictionary with the transferred molecule index as key
    and a list [process between molecules, transfer_rate] as argument.
    :param center_index: Index of the studied exciton
    Two list are built:
        process_list collects all the possible processes
        rates_list collects the respective rates
        The indexes for each process coincides with the index in rates
    We implement the KMC algorithm to rates_list.
    :return: path: chosen process:
                    if it is a transfer process: list(transferred molecule index, process)
                    if it is a decay process:  decay_process
             rates_list
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

    process_index = select_process(rates_list)
    path = process_list[process_index]

    return path, rates_list


def select_process(constant_list):
    """
    :param constant_list: List with the constant rates
    :return: Chooses a position of the list proportionally to its value.
    """
    r = np.sum(constant_list) * rd.random()
    list_ = np.where(r > np.cumsum(constant_list))
    return len(list_[0])


def time_advance(rate_list, centers_index):
    """
    :param rate_list: List with all the rates. Considering all the processes for all exciton
    :param centers_index: List with the indexes of every excited molecule.
    :return: Process duration for each excitation.
    """
    center_time = {}
    for index in centers_index:
        r = rd.random() + 0.0001        # In order to exclude r = 0 we add a quite small quantity.
        center_time[str(index)] = (-np.log(r))/np.sum(rate_list)

    return center_time
