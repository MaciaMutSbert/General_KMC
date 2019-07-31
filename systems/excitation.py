import numpy as np

"""
Plantejam una sèrie de funcions per excitar el sistema. Estan dissenyades a partir del supòsit que 
la informació sobre els excitons ve donada per un diccionari tipus:
    clau: referència de l'estat excitat (han d'anar en consonància amb els estats vàlids on pot estar la molecula)
    argument: llista de longitud = número d'excitons amb les posicions de cada un.
"""


def excited_system(molecules, excitons, tolerance):
    """
    :param molecules: List of the defined molecules
    :param excitons: Information about the desired excitation
    :return: list with the index of the excited molecules
    The function modifies the list of molecules with the excitons.
    """
    centre_list = []

    for state in excitons:
        for position in excitons[state]:
            index = get_excited_index(position, centre_list, molecules, tolerance)

            molecules[index].change_state(state)     # check the correct method
            centre_list.append(index)

    return centre_list


def get_excited_index(position, centre_list, molecules, tolerance):
    """
    This function simply redirects the program to more specific functions
    """

    if position is 'random':
        return pick_random(centre_list, molecules)

    if position is 'first':
        return pick_first(centre_list, molecules)

    if position is 'last':
        return pick_last(centre_list, molecules)

    if position is 'centre':
        return pick_centre(centre_list, molecules, tolerance)


############################################


def pick_random(centre_list, molecules):
    """
    :param centre_list: list with the indexes of the excited molecules
    :param molecules: list with all the instances of class Molecule
    :return: a random position of the list molecules
    """
    picked = True
    while picked is True:
        index = np.random.randint(0, len(molecules))

        if index in centre_list:
            picked = True
        else:
            picked = False
    return index


def pick_first(centre_list, molecules):
    """
    :param centre_list: list with the indexes of the excited molecules
    :param molecules: list with all the instances of class Molecule
    :return: the first element of the list. If this element is already taken picks another randomly
    """
    index = 0
    if index in centre_list:
        index = pick_random(centre_list, molecules)
    return index


def pick_last(centre_list, molecules):
    """
    :param centre_list: list with the indexes of the excited molecules
    :param molecules: list with all the instances of class Molecule
    :return: the last element of the list. If this element is already taken picks another randomly
    """
    index = len(molecules)
    if index in centre_list:
        index = pick_random(centre_list, molecules)
    return index


def pick_centre(centre_list, molecules, tolerance):
    """
    :param centre_list: list with the indexes of the excited molecules
    :param molecules: list with all the instances of class Molecule
    :param tolerance: Since the centre may not coincide with the 0 point we give this extra parameter.
    The function looks for a molecule in a cercle of radius = tolerance and it is considered the centre.
    It the centre is already taken or there is not any molecule in the centre of the distribution,
    as it migth happen in an amorphous material, the function picks another index randomly.
    It is taken as the lattice_parameter / 2 in ordered systems and (COM S'AGAFA) in disordered systems.
    :return: the index of the excited molecule
    """
    index = None
    for i, molecule in enumerate(molecules):
        position = np.abs(molecule.molecular_coordinates())
        if position <= tolerance:
            index = i

    if index in centre_list:
        index = pick_random(centre_list, molecules)

    if index is None:
        print('Not centre found. Randomly located exciton')
        index = pick_random(centre_list, molecules)

    return index



