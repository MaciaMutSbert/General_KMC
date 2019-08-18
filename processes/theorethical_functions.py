import numpy as np
from numpy import pi
from molecules import Molecule
from conversion_functions import from_ev_to_au, from_ns_to_au


def theorethical_diffusion_values(system_information):
    """
    :param system_information: dictionary with the system information (parameter to parameter)
    :return: dictionary with the theorethical lifetime, diffusion constant and length. (In the cases that they can be
    defined).
    """

    possible_cases = [['ordered', 'parallel'], ['2d_ordered', 'parallel'],['3d_ordered', 'parallel'],
                      ['ordered', 'antiparallel'], ['2d_ordered', 'antiparallel'],['3d_ordered', 'antiparallel']]

    boltzmann_constant = 8.617333 * 10 ** (-5)  # Boltzmann constant in eV * K^(-1)
    k = 2
    r = system_information['lattice_parameter']
    T = system_information['conditions']['temperature']
    n = system_information['conditions']['refractive_index']
    relax = system_information['relaxation_energies']['s_1']
    transition_moment_vector = system_information['transition_moment']
    transition_moment = np.linalg.norm(transition_moment_vector)

    generic_molecule = Molecule(state_energies=system_information['state_energies'],
                                state='s_1', reorganization_energies=system_information['relaxation_energies'],
                                transition_moment=transition_moment_vector)
    decay_rates = generic_molecule.decay_rates()
    decay_sum = 0
    for decay in decay_rates:
        decay_sum =+ decay_rates[decay]
    life_time = 1 / decay_sum

    d = len(system_information['dimensions'])


    if [system_information['order'], system_information['orientation']] in possible_cases:

        factor_1 = 2 * pi * (transition_moment**2 * k**2 / ((r/0.053)**3 * n**2 ))**2
        factor_2 = np.sqrt(1 / (4 * pi * relax * boltzmann_constant * T)) * np.exp(- relax / (4*boltzmann_constant*T))
        factor_2_transformed = from_ev_to_au(factor_2, 'inverse')

        rate_au = factor_1 * factor_2_transformed

        rate = from_ns_to_au(rate_au, 'direct')
        D = rate * r**2
        Ld = np.sqrt(2 * d * rate * r**2 * life_time)

    else:
        D = None
        Ld = None




    return {'life_time': life_time, 'diffusion_constant':  D,
            'diffusion_length': Ld}
