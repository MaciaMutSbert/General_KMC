import numpy as np
import random as rd


def kmc_algorithm(decay_rates, transfer_rates):
    """
    :param decay_rates: Dictionary with decay process as keys and decay rates as arguments
    :param transfer_rates: Dictionary with the transferred molecule index as key
    and a list [process between molecules, transfer_rate] as argument.
    Two list are built:
        process_list collects all the possible processes
        rates_list collects the respective rates
        The indexes for each process coincides with the index in rates
    We implement the KMC algorithm to rates_list.
    :return: path: chosen process:
                    if it is a transfer process: list(process, transferred molecule index)
                    if it is a decay process:  decay_process
             time: advanced time
    """
    process_list = []
    rates_list = []

    for decay in decay_rates:
        process_list.append(decay)
        rates_list.append(decay_rates[decay])

    for neighbour_index in transfer_rates:
        for process in transfer_rates[neighbour_index]:
            process_list.append([int(neighbour_index), process])
            rates_list.append(transfer_rates[neighbour_index][process])

    process_index = select_process(rates_list)
    path = process_list[process_index]

    time = time_advance(rates_list)

    return path, time


def select_process(constant_list):
    """
    CHOOSE A PATHWAY
    constant_list: rate constants list
    return: chosen path index
    """
    r = np.sum(constant_list) * rd.random()
    list_ = np.where(r > np.cumsum(constant_list))
    return len(list_[0])


def time_advance(constant_list):
    """
    TIME ADVANCE
    constant_list: list of rate constants
    return: advanced time
    """
    r = rd.random()
    return -np.log(r)/np.sum(constant_list)