import numpy as np
import random as rd


def kmc_algorithm(rate_list, process_list):
    """
    :param rate_list: List with all the computed rates for all the neighbours for all the centers
    :param process_list: The associated process to each rate and the new molecule affected
    :return:    plan: The chosen proces and the new molecule affected
                time: the duration of the process
    """

    process_index = select_process(rate_list)
    chosen_process = process_list[process_index]

    time = time_advance(rate_list)

    return chosen_process, time


def select_process(constant_list):
    """
    :param constant_list: List with the constant rates
    :return: Chooses a position of the list proportionally to its value.
    """
    r = np.sum(constant_list) * rd.random()
    list_ = np.where(r > np.cumsum(constant_list))
    return len(list_[0])


def time_advance(rate_list, ):
    """
    :param rate_list: List with all the rates. Considering all the processes for all exciton
    :return: Process duration
    """
    r = rd.random()
    return (-np.log(r))/np.sum(rate_list)
