import numpy as np
from scipy import distance
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


def compute_forster_coupling(molecule1, molecule2, conditions):
    """
    :param molecule1: excited moleculen. Donor
    :param molecule2: neighbouring molecule. Possible acceptor
    :param conditions: dictionary with physical conditions
    :return: Förster coupling between both molecules. We don't implement any correction for short distances.
    """

    u_D = molecule1.get_transition_moment()
    u_A = molecule2.get_transition_moment()
    r = intermolecular_distance(molecule1, molecule2)       # a. u.
    n = conditions['refractive_index']

    info = str(hash((u_D, u_A, r, n)))
    if info in coupling_memory:
        forster_coupling = coupling_memory[info]

    else:
        k = orientational_factor(u_D, u_A, r)
        forster_coupling = k**2 * np.dot(u_D, u_A) / (n**2 * r**3)
        coupling_memory[info] = forster_coupling

    return forster_coupling


##########################################################################################
#                               DICTIONARY COUPLINGS
##########################################################################################

couplings = {'s_1_g_s': compute_forster_coupling}



###############################################################################################
#                            AUXILIAR FUNCTIONS
###############################################################################################


def intermolecular_distance(molecule1, molecule2):
    """
    :param molecule1: donor
    :param molecule2: acceptor
    :return: the euclidean distance between the donor and the acceptor
    """

    position_D = molecule1.molecular_coordinates()
    position_A = molecule2.molecular_coordinates()
    r = distance.euclidean(position_D, position_A)
    return from_nm_to_au(r, 'direct')


def orientational_factor(u_D, u_A, r):
    """
    :param u_D: dipole transition moment of the donor
    :param u_A: dipole transition moment of the acceptor
    :param r:  intermolecular_distance
    :return: the orientational factor between both molecules
    """

    nd = unity(u_D)
    na = unity(u_A)
    e = unity(r)

    return np.dot(nd, na) - 3*np.dot(e, nd)*np.dot(e, na)


def unity(vector):
    """
    :param vector:
    :return: computes a unity vector in the direction of vector
    """
    return vector / np.linalg.norm(vector)