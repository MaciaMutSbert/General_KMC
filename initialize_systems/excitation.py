import numpy as np


def excited_system(molecules, initial_excitation, tolerance):
    """
    :param molecules: List of the defined molecules
    :param initial_excitation: Information about the desired excitation
    :param tolerance: (only used in pick_centre()).
    :return: list with the index of the excited molecules
    The function modifies the list of molecules with the excitons.
    """
    centre_list = []

    for state in initial_excitation:
        for position in initial_excitation[state]:
            index = get_excited_index(position, centre_list, molecules, tolerance)

            molecules[index].set_state(state)     # check the correct method
            centre_list.append(index)

    return centre_list


########################################################################################


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

    if position is 'furthest':
        return pick_furthest(centre_list, molecules)


########################################################################################


def pick_random(centre_list, molecules):
    """
    :param centre_list: list with the indexes of the excited molecules
    :param molecules: list with all the instances of class Molecule
    :return: a random position of the list molecules
    """

    picked = True                   # the loop will work as long as a dexcited molecule is not chosen. If the chosen
    counter = 1                     # molecule is excited then picked = True and the loop will start again.
    while picked is True:

        index = np.random.randint(0, len(molecules))            # random position of the list molecules

        if index in centre_list:
            picked = True

            counter =+ 1
            if counter == len(molecules):                   # flux control. Stops the while loop if all molecules are
                print('All molecules excited')              # excited
                break

        else:
            picked = False
            break

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
    It is taken as the average over the lattice parameters in ordered systems and as 0.2 nm in disordered systems
    :return: the index of the excited molecule
    """

    index = None
    for i, molecule in enumerate(molecules):

        position = np.linalg.norm(molecule.molecular_coordinates())
        if position <= tolerance:                   # checks if the chosen position if in a cercle of radius = tolerance
            index = i

    # if any other case occurs (the center is already excited or there is not a molecule in the nearby of (0,0,0)) a
    # random molecule from the list is chosen

    if index in centre_list:
        index = pick_random(centre_list, molecules)

    if index is None:
        print('Not centre found. Randomly located exciton')
        index = pick_random(centre_list, molecules)

    return index


def pick_furthest(centre_list, molecules):
    """
    :param centre_list: list with the indexes of the excited molecules
    :param molecules: list with all the instances of class Molecule
    :return: the furthest molecule from the centre.
    """

    furthest_position = 0
    index = None
    for i, molecule in enumerate(molecules):

        position = np.linalg.norm(molecule.molecular_coordinates())
        if position > furthest_position:            # checks if the chosen position is further than the furthest at
            furthest_position = position            # the moment.
            index = i

    # if the furthest molecule from the center in the system is already excited, picks another randomly.

    if index in centre_list:
        index = pick_random(centre_list, molecules)

    return index

