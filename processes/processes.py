import numpy as np
from scipy.spatial import distance
from numpy import pi


def get_transfer_rates(center, neighbour_index, system, rate_memory, molecular_memory):
    """
    :param center: Index of the studies excited molecule
    :param neighbour_index: Index of a nearby molecule (possible acceptor of the exciton)
    :param system: Dictionary with the list of molecules and additional physical information
    :return: Dictionary with the possible transfer processes between the two involved molecules
    as keys and its rates as arguments.
    """
    molecules = system['molecules']
    transfer_rates = {}
    molecular_info = []

    if molecules[center].type == 1:
        if molecules[neighbour_index].type == 1:

            molecular_info.append(molecules[center].type)
            molecular_info.append(molecules[neighbour_index].type)

            k = system['conditions']['orientational_factor_1']
            n = system['conditions']['refractive_index']
            alfa = 1.15  # short distances correction
            sigma = system['conditions']['a_e_spectra_deviation'] / 27.211  # atomic units
            delta = system['conditions']['delta'] / 27.211

            if molecules[center].state == 1:
                if molecules[neighbour_index].state == 0:
                    process = 'Singlet_transfer'

                    u = system['conditions']['transition_dipole']
                    center_position = np.array(molecules[center].coordinates)
                    neighbour_position = np.array(molecules[neighbour_index].coordinates)
                    inter_distance = distance.euclidean(center_position, neighbour_position)/0.053    # atomic units
                    """
                    Per ara el rate sols depen de la distància intermolecular, doncs la resta de paràmetres venen
                    fixats per les condicions del material. 
                    Per això estaria bé que guardar el rate calculat per una determinada distància a fi de salvar cost
                    computacional.
                    En una mateixa simulació en general tendrem paràmetres
                    del material que no variaran. Sols variaran el moment de transició dipolar (en quan a modul) i 
                    les posicions relatives. 
                    """

                    molecular_info.append(molecules[center].state)
                    molecular_info.append(molecules[neighbour_index.state])
                    molecular_info.append(u)
                    molecular_info.append(inter_distance)

                    if molecular_info in molecular_memory:
                        index = molecular_memory.index(molecular_info)
                        rate = rate_memory[index]
                        transfer_rates[process] = rate

                    else:
                        factor_1 = k**2 * pi**(3/2) * u**4 * np.exp(-(delta**2) / (2*sigma)**2)
                        factor_2 = n**4 * (alfa*u + inter_distance)**6 * sigma

                        rate = factor_1 / factor_2      # atomic units
                        transfer_rates[process] = rate / (2.4189*10**-8)      # in ns⁻¹

                        molecular_memory.append(molecular_info)
                        rate_memory.append(rate)

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






