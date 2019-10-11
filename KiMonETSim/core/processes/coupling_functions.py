import numpy as np
from KiMonETSim.conversion_functions import from_nm_to_au


##########################################################################################
#                                   COUOPLING FUNCTIONS
##########################################################################################

coupling_memory = {}


def compute_forster_coupling(donor, acceptor, conditions):
    """
    :param donor: excited molecules. Donor
    :param acceptor: neighbouring molecule. Possible acceptor
    :param conditions: dictionary with physical conditions
    :return: Förster coupling between both molecules. We don't implement any correction for short distances.
    """

    # definition of the parameters of the donor and acceptor needed in the calculation of the Förster coupling

    u_d = donor.get_transition_moment()                     # transition dipole moment (donor) a.u
    u_a = acceptor.get_transition_moment()                  # transition dipole moment (acceptor) a.u
    momentum_projection = np.dot(u_d, u_a)                  # inner product of both a.u

    r_vector = intermolecular_vector(donor, acceptor)       # position vector between donor and acceptor
    r = np.linalg.norm(r_vector)                            # inter molecular distance a.u

    n = conditions['refractive_index']                      # refractive index of the material

    info = str(hash((momentum_projection, r, n, 'förster')))
    # we define a compact string with the characteristic information of the coupling

    if info in coupling_memory:
        # use of the memory for computed situations
        forster_coupling = coupling_memory[info]

    else:
        k = orientational_factor(u_d, u_a, r_vector)            # orientational factor between molecules

        forster_coupling = k**2 * momentum_projection / (n**2 * r**3)       # electronic coupling a.u

        coupling_memory[info] = forster_coupling                            # memory update for new couplings

    return forster_coupling


##########################################################################################
#                               DICTIONARY COUPLINGS
##########################################################################################

couplings = {'s1_gs': compute_forster_coupling}


###############################################################################################
#                            AUXILIAR FUNCTIONS
###############################################################################################


def intermolecular_vector(molecule1, molecule2):
    """
    :param molecule1: donor
    :param molecule2: acceptor
    :return: the euclidean distance between the donor and the acceptor
    """
    position_d = molecule1.molecular_coordinates()
    position_a = molecule2.molecular_coordinates()
    r = position_d - position_a
    return from_nm_to_au(r, 'direct')


def orientational_factor(u_d, u_a, r):
    """
    :param u_d: dipole transition moment of the donor
    :param u_a: dipole transition moment of the acceptor
    :param r:  intermolecular_distance
    :return: the orientational factor between both molecules
    """
    nd = unity(u_d)
    na = unity(u_a)
    e = unity(r)
    return np.dot(nd, na) - 3*np.dot(e, nd)*np.dot(e, na)


def unity(vector):
    """
    :param vector:
    :return: computes a unity vector in the direction of vector
    """
    return vector / np.linalg.norm(vector)
