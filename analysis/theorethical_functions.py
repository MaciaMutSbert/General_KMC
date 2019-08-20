import numpy as np
from numpy import pi
from molecules import Molecule
from conversion_functions import from_ev_to_au, from_ns_to_au


def theoretical_diffusion_values(system_information, molecule_information):
    """
    :param system_information: dictionary with the system information (parameter to parameter)
    :param molecule_information: dictionary with the internal properties of the molecule
    :return: dictionary with the theoretical lifetime, diffusion constant and length. (In the cases that they can be
    defined).
    """

    # system information
    boltzmann_constant = 8.617333 * 10 ** (-5)  # Boltzmann constant in eV * K^(-1)
    k = 2                         # orientational factor. Always 2 in the systems where this model is possible

    r = system_information['lattice']['lattice_parameter'][0]       # equall in all directions (nm)

    T = system_information['conditions']['temperature']             # temperature (K)
    n = system_information['conditions']['refractive_index']        # refractive index

    d = len(system_information['lattice']['dimensions'])                       # dimensionality of the system

    ############################################################################################################

    # molecule information.     Again a generic object of molecule is initialized

    excited_state = 's1'
    # in a diffusion study there will be only on exciton, so only one key

    reorganization = molecule_information['reorganization_energies'][excited_state]
    # reorganization energy of the excited electronic state. eV

    transition_moment_vector = molecule_information['transition_moment']          # transition dipole moment (a.u)
    transition_moment = np.linalg.norm(transition_moment_vector)                # modulus

    # initialization of a generic instance of class molecule
    generic_molecule = Molecule(state_energies=molecule_information['state_energies'],
                                state='s1', reorganization_energies=molecule_information['reorganization_energies'],
                                transition_moment=transition_moment_vector)

    ############################################################################################################

    # The life time of the exciton is taken as the inverse of the sum of all decay rates
    # A theoretical value for the life time will be always computed.
    decay_rates = generic_molecule.decay_rates()
    decay_sum = 0
    for decay in decay_rates:
        print(decay_rates[decay])
        decay_sum =+ decay_rates[decay]
    life_time = 1 / decay_sum
    print(life_time)

    ############################################################################################################

    # cases in which there is a theoretical model for the diffusion
    possible_cases = [['1d_ordered', 'parallel'], ['2d_ordered', 'parallel'], ['3d_ordered', 'parallel'],
                      ['1d_ordered', 'antiparallel'], ['2d_ordered', 'antiparallel'], ['3d_ordered', 'antiparallel']]

    # the theoretical values for D and LD will always be computed for the possible cases:
    if [system_information['type'], system_information['orientation']] in possible_cases:

        factor_1 = 2 * pi * (transition_moment**2 * k**2 / ((r/0.053)**3 * n**2))**2    # in atomic units
        factor_2 = np.sqrt(1 / (4 * pi * reorganization * boltzmann_constant * T)) * \
                   np.exp(- reorganization / (4*boltzmann_constant*T))                  # in eV
        factor_2_transformed = from_ev_to_au(factor_2, 'inverse')                       # in atomic units

        rate_au = factor_1 * factor_2_transformed                   # rate in atomic units (time⁻¹)
        rate = from_ns_to_au(rate_au, 'direct')                     # rate in ns⁻¹

        D = rate * r**2                                             # diffusion constant (nm² ns⁻¹)
        Ld = np.sqrt(2 * d * rate * r**2 * life_time)               # diffusion length (nm)

    else:                                                           # Ld and D taken as None if there is no
        D = None                                                    # theoretical model
        Ld = None

    # return the 3 parameters in a dictionary
    return {'lifetime': life_time, 'diffusion_constant':  D,
            'diffusion_length': Ld}
