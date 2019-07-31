import numpy as np

"""
EXCITACIÓ DEL SISTEMA:
Plantejam una sèrie de funcions que a partir d'un diccionari que contengui la informació sobre els excitons, modificaran
la llista de molecules excitant les corresponents. Per un bon funcionament aquest diccionari ha de tenir la forma:
    :key: referència de l'estat excitat (han d'anar en consonància amb els estats vàlids on pot estar la molecula)
    :argument: llista de longitud = número d'excitons amb les posicions de cada un. Les posicions són: 'random', 
    'first', 'last' i 'centre'.

La funció general i que es cridarà en el programa que inicialitza el sistema és excited_system(). Arguments:
    molecules: llista de les instàncies tipus molecula
    excitons: diccionari amb l'informació sobre els excitons que volem generar (amb la forma esmentada)
    tolerance: paràmetre que usarem alhora de detectar molècules donades unes coordenades (s'explicarà més endavant)
En un doble bulce (sobre el diccionari i sobre la llista que té com argument) agafa un índex de la llista molècules i 
l'excita. En cada iteració actualitza la llista dels indexs de les molecules excitades que retornarà.

Funció get_excited_index(position, centre_list, molecules, tolerance). De fet aquesta funció sols redirigeix a altres 
funcions segons el paràmetre position. Retorna l'índex triat.

pick_random(centre_list, molecules). Tria un index aleatori de la llista molecules. Passant-li la centre_list 
asseguram que no es repeteixen indexes.

pick_first(centre_list, molecules) i l'anàlega pick_last. Prenen el primer o últim element de la llista. Si aquest ja 
hagués estat triat (suposem en un random) es tria un índex de manera aleatoria. Val a dir, si tenim un sistema ordenat
aquestes posicions seran cantons de la xarxa si el material és amorf el resultat no serà diferent de demanar una posició 
random. Pentura s'hi guanyaria eficiència.

pick_centre(centre_list, molecules, tolerance). Recórre tota la llista de molècules i agafa la que estigui situada en una
esfera de radi tolerance centrada en 0. Si tenim un sistema ordenat amb un paràmetre de xarxa donat, la tolerància n'és
la meitat. Per un sistema amorf la definim com la 30èssima part de la mitjana dels 3 costats.
Si no tenim cap molècula en aquest entorn (podria donar-se en un sistema amorf) o aquesta ja està excitada pren un
índex aleatòriament, tot printant un avís de que no tenim cap molècula al centre.

pick_furthest(centre_list, molecules). Recórre tota la llista de molècules i agafa la que estigui situada més lluny del 
punt 0. En el cas de que ja estigui excitada n'agafa una de manera aleatòria.
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

    if position is 'furthest':
        return pick_furthest(centre_list, molecules)


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
    It is taken as the lattice_parameter / 2 in ordered systems and as the 30th part of the average length of
     the system in disordered systems.
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


def pick_furthest(centre_list, molecules):
    """
    :param centre_list: list with the indexes of the excited molecules
    :param molecules: list with all the instances of class Molecule
    :return: the furthest molecule from the centre.
    """
    furthest_position = 0
    index = None
    for i, molecule in enumerate(molecules):
        position = np.abs(molecule.molecular_coordinates())

        if position > furthest_position:
            furthest_position = position
            index = i

    if index in centre_list:
        index = pick_random(centre_list, molecules)

    return index

