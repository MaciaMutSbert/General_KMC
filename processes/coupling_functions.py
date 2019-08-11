import numpy as np
from scipy.spatial import distance
from conversion_functions import from_nm_to_au


"""
COUPLING FUNCTIONS:
En aquest arxiu contruïm un diccionari amb les funcions per calcular els acoplaments electrònics entre el donant i
l'acceptor. Totes aquestes funcions rebran els mateixos arguments: molecule1 (donor), molecule2 (acceptor) i el diccionari amb 
les condicions físiques del sistema (conditions). Consta de 3 parts diferenciades:
        Definició de les funcions que calcularan els diferents acoplaments electrònics
        Definició del diccionari i les seves claus
        Definició d'una sèrie de funcions auxiliars.
A més definim un diccionari, coupling_memory, que anirà guardant els acoplaments calculats i serà usat en cas de que
ja estigui calculat.

Acoplaments considerats:
    Acoplament de Förster entre una molècula en l'estat singlet 1 i una en el singlet 0 (estat fonamental).
        Etiqueta: s_1_g_s
        Funció: compute_forster_coupling
"""

##########################################################################################
#                                   COUOPLING FUNCTIONS
##########################################################################################

coupling_memory = {}


def compute_forster_coupling(donor, acceptor, conditions):
    """
    :param donor: excited moleculen. Donor
    :param acceptor: neighbouring molecule. Possible acceptor
    :param conditions: dictionary with physical conditions
    :return: Förster coupling between both molecules. We don't implement any correction for short distances.
    """

    u_d = donor.get_transition_moment()
    print(u_d)
    u_a = acceptor.get_transition_moment()
    momentum_projection = np.dot(u_d, u_a)                  # a. u.
    r_vector = intermolecular_vector(donor, acceptor)       # a. u.
    r = np.linalg.norm(r_vector)                            # a. u.
    n = conditions['refractive_index']

    info = str(hash((momentum_projection, r, n, 'förster')))
    if info in coupling_memory:
        forster_coupling = coupling_memory[info]

    else:
        k = orientational_factor(u_d, u_a, r_vector)
        forster_coupling = k**2 * momentum_projection / (n**2 * r**3)
        coupling_memory[info] = forster_coupling

    return forster_coupling


##########################################################################################
#                               DICTIONARY COUPLINGS
##########################################################################################

couplings = {'s_1_g_s': compute_forster_coupling}


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
