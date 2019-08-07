import numpy as np
from scipy.spatial import distance
from systems.excitation import excited_system
import copy


def get_disordered_system(conditions,  # External conditions of the system such as temperature
                          generic_molecule,
                          dimensions=[10, 10],
                          number_molecules=1000,
                          orientation='parallel',
                          reference_orientation=[1,0,0],
                          excitons={'s_1': ['centre']},):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param generic_molecule: generic instance of class Molecule with the intern information defined.
        In the function a position and orientation for the molecule are given
    :param dimensions: Number of molecules per side. A list with len = num_dimensions.
    :param orientation: String parameter. Indicates whether the orientation of the molecules is parallel,
    random or antiparallel (each molecule finds its 1sts neighbours antiparallely orientated)
    :param reference_orientation: 3 length string (vector), indicates the privileged orientation.
    :param number_molecules: Number of molecules in the system.
    :param excitons: Dictionary with the information of the excitons (type and a position for each).
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """
    intermolecular_distance = 0.1       # medium intermolecular distance, 0.1 nm
    physical_dimensions = np.array(dimensions) * intermolecular_distance
    capacity = get_capacity(dimensions)

    if orientation is 'antiparallel':
        print('Not antiparallel orientation considered for an amorphous material')
        return

    molecules = []
    molecule_count = 0
    while molecule_count <= number_molecules:

        too_close = True
        while too_close is True:
            coordinates = []
            for physical_dimension in physical_dimensions:
                x = physical_dimensions*np.random.random() - physical_dimension/2
                # We want the distribution center at (0,0,0)
                coordinates.append(x)

            too_close = distance_checking(coordinates, molecules)

        orientation_vector = get_orientation(orientation, reference_orientation, len(dimensions), pointing=1)
        molecule = copy.deepcopy(generic_molecule)
        molecule.initialize_coordinates(coordinates)
        molecule.initialize_orientation(orientation_vector)
        molecules.append(molecule)
        molecule_count += 1

        if molecule_count == capacity:
            break

    tolerance = intermolecular_distance * 2
    centre_indexes = excited_system(molecules, excitons, tolerance)
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions, 'centres': centre_indexes}

    return system


def get_1d_ordered_system(conditions,
                          generic_molecule,
                          dimensions=[10],
                          lattice_parameter=1,
                          orientation='parallel',
                          reference_orientation=[1, 0, 0],
                          excitons={'s_1': ['centre']},):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param generic_molecule: generic instance of class Molecule with the intern information defined.
        In the function a position and orientation for the molecule are given
    :param dimensions: Number of molecules per side. A list with len = num_dimensions.
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param orientation: String parameter. Indicates whether the orientation of the molecules is parallel,
    random or antiparallel (each molecule finds its 1sts neighbours antiparallely orientated)
    :param reference_orientation: 3 length string (vector), indicates the privileged orientation.
    :param excitons: Dictionary with the information of the excitons (type and a position for each).
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """
    physical_dimensions = np.array(dimensions) * lattice_parameter    # nm

    if check_lattice(lattice_parameter, generic_molecule) is False:
        print('Lattice parameter is smaller than molecular characteristic length')
        return

    molecules = []
    x_max = physical_dimensions[0]/2
    pointing = 1
    symmetry = get_symmetry[orientation]
    # We want the distribution center at 0
    for x in np.arange(-x_max, x_max, lattice_parameter):
        coordinates = [x, 0, 0]
        orientation_vector = get_orientation(orientation, reference_orientation, len(dimensions), pointing)

        molecule = copy.deepcopy(generic_molecule)
        molecule.initialize_coordinates(coordinates)
        molecule.initialize_orientation(orientation_vector)
        molecules.append(molecule)

        pointing = pointing * symmetry[0]

    centre_indexes = excited_system(molecules, excitons, lattice_parameter / 2)
    conditions['lattice_parameter'] = lattice_parameter
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions, 'centres': centre_indexes}

    return system


def get_2d_ordered_system(conditions,
                          generic_molecule,
                          dimensions=[10, 10],
                          lattice_parameter=1,
                          orientation='parallel',
                          reference_orientation=[1, 0, 0],
                          excitons={'s_1': ['centre']},):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param generic_molecule: generic instance of class Molecule with the intern information defined.
        In the function a position and orientation for the molecule are given    :param dimensions: Dimensions of the system. A list with len = 3. By default [10, 10, 0]
    :param dimensions: Number of molecules per side. A list with len = num_dimensions.
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param orientation: String parameter. Indicates whether the orientation of the molecules is parallel,
    random or antiparallel (each molecule finds its 1sts neighbours antiparallely orientated)
    :param reference_orientation: 3 length string (vector), indicates the privileged orientation.
    :param excitons: Dictionary with the information of the excitons (type and a position for each).
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """
    physical_dimensions = np.array(dimensions) * lattice_parameter    # nm

    if check_lattice(lattice_parameter, generic_molecule) is False:
        print('Lattice parameter is smaller than molecular characteristic length')
        return

    molecules = []
    x_max = physical_dimensions[0]/2
    y_max = physical_dimensions[1]/2
    # We want the center of the distribution at (0,0)
    symmetry = get_symmetry[orientation]
    x_count = 0
    for x in np.arange(-x_max, x_max, lattice_parameter):
        pointing = (-1)**x_count
        for y in np.arange(-y_max, y_max, lattice_parameter):
            coordinates = [x, y, 0]
            orientation_vector = get_orientation(orientation, reference_orientation, len(dimensions), pointing)

            molecule = copy.deepcopy(generic_molecule)
            molecule.initialize_coordinates(coordinates)
            molecule.initialize_orientation(orientation_vector)
            molecules.append(molecule)

            pointing = pointing * symmetry[0]
        x_count = x_count + symmetry[1]

    centre_indexes = excited_system(molecules, excitons, lattice_parameter / 2)
    conditions['lattice_parameter'] = lattice_parameter
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions, 'centres': centre_indexes}

    return system


def get_3d_ordered_system(conditions,
                          generic_molecule,
                          dimensions=[10, 10, 10],
                          lattice_parameter=1,
                          orientation='parallel',
                          reference_orientation=[1, 0, 0],
                          excitons={'s_1': ['centre']},):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param generic_molecule: generic instance of class Molecule with the intern information defined.
        In the function a position and orientation for the molecule are given    :param dimensions: Dimensions of the system. A list with len = num_dimensions. By deafult [10, 10, 10]
    :param dimensions: Number of molecules per side. A list with len = num_dimensions.
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Dictionary with the information of the excitons (type and a position for each).
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """
    physical_dimensions = np.array(dimensions) * lattice_parameter    # nm

    if check_lattice(lattice_parameter, generic_molecule) is False:
        print('Lattice parameter is smaller than molecular characteristic length')
        return

    molecules = []
    x_max = physical_dimensions[0]/2
    y_max = physical_dimensions[1]/2
    z_max = physical_dimensions[2]/2
    step = lattice_parameter
    # We want the center of the distribution at (0,0,0)

    x_count = 0
    symmetry = get_symmetry[orientation]
    for x in np.arange(-x_max, x_max, step):
        x_pointing = (-1)**x_count
        y_count = 0
        for y in np.arange(-y_max, y_max, step):
            pointing = x_pointing * (-1)**y_count
            for z in np.arange(-z_max, z_max, step):
                coordinates = [x, y, z]
                orientation_vector = get_orientation(orientation, reference_orientation, len(dimensions), pointing)

                molecule = copy.deepcopy(generic_molecule)
                molecule.initialize_coordinates(coordinates)
                molecule.initialize_orientation(orientation_vector)
                molecules.append(molecule)

                pointing = pointing*symmetry[0]
            y_count = y_count + symmetry[1]
        x_count = x_count + symmetry[1]

    centre_indexes = excited_system(molecules, excitons, lattice_parameter / 2)
    conditions['lattice_parameter'] = lattice_parameter
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions, 'centres': centre_indexes}

    return system


######################################################################################
#                           DICTIONARIES
######################################################################################

get_ordered_system = {'1': get_1d_ordered_system,
                      '2': get_2d_ordered_system,
                      '3': get_3d_ordered_system}

get_homogeneous_system = {'ordered': get_ordered_system,
                          'disordered': get_disordered_system}


######################################################################################
#                              AUXILIAR FUNCTIONS


def distance_checking(coordinates, molecules):
    """
    :param coordinates: coordinates of the studied molecule
    :param molecules: list of the molecules already defined
    :return: Boolean. Indicates if the new molecule is too close to some of the already defined.
    """
    coordinates = np.array(coordinates)

    for molecule in molecules:
        if distance.euclidean(coordinates, molecule.molecular_coordinates()) < molecule.characteristic_length:
            return True

    return False


def get_capacity(dimensions):
    """
    :param dimensions: Number of molecules per side
    :param number_molecules: total number of molecules given
    :return: Integer. Maximus number of molecules that could fit in the system.
    """
    capacity = 1
    for dimension in dimensions:
        capacity = capacity*dimension
    return capacity




def check_lattice(lattice_parameter, generic_molecule):
    """
    :param lattice_parameter: Lattice parameter
    :param generic_molecule: instance of class molecule with all its natural parameters defined
    :return: Boolean. Checks if the lattice parameter chosen is smaller than the molecular characteristic length.
    """
    if lattice_parameter < generic_molecule.characteristic_length:
        return False
    else:
        return True


"""
ANTIPARAL·LELISME. COM L'IMPLEMENTAM?   
Hem inclòs en la inicialització dels sistemes ordenats la possibilitat de què les molècules estiguin ordenades
antiparal·lelament. L'entenem com un antiparal·lelisme entre veïns primers. O sigui, donada una molècula amb orientació
'+', els 2, 4 o 6 veïns primers de la xarxa tendran orientació '-'.
Val a dir que per un material amorf no es contempla l'antiparal·lelisme. Si es donàssim com a paràmetres:
    order = 'disordered'
    orientation = 'antiparallel'
Sortiria un avís per pantalla indicant la no viabilitat d'aquesta distribució i no s'inicialitzaria cap sistema.
"""

get_symmetry = {'parallel': [1, 2], 'antiparallel': [-1, 1], 'random': [0,0]}


def get_orientation(orientation, reference_orientation, dimensionality, pointing):

    if orientation is 'random':
        if dimensionality == 3:
            phi = 2*np.pi() * np.random.rand()
            theta = np.pi * np.random.rand()
            return [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]

        else:
            phi = 2*np.pi() * np.random.rand()
            return [np.cos(phi), np.sin(phi), 0]

    else:
        return reference_orientation * pointing
