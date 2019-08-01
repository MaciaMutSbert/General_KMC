import numpy as np
from scipy.spatial import distance
from systems.molecules import Molecule
from systems.excitation import excited_system

"""
def get_system(morphology, excitons):

    if morphology['Homogeneity'] is True:
        return get_homogeneous_system(morphology, excitons)

    else:
        return get_inhomogeneous_system(morphology, excitons)
"""


def get_homogeneous_system(conditions,
                           generic_molecule,
                           order='ordered',
                           dimensionality=2,
                           dimensions=[10, 10, 0],
                           num_molecules=0,             # Only if order = Disordered
                           lattice_parameter=0.1,       # Only if order = Ordered
                           orientation='parallel',
                           reference_orientation=[1, 0, 0],
                           excitons={'s_1': ['centre']}):
    """
    PER COMENTARIS GUARDAR ESCRITS DIPC
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param generic_molecule: generic instance of class Molecule with the nature parameters already defined.
    Only position and orientation have to be defined
    :param num_molecules: Number of molecules in the system. Only needed for a disordered system.
    :param dimensionality: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = num_dimensions
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :param order: String argument. Indicates if we want an ordered system (crystal) or an amorphous system.
    :return: A dictionary with a list of molecules, an updated dictionary with the physical conditions
    and a list with the indexes of the excited molecules.
    """

    if order == 'ordered':
        return get_ordered_system(conditions, generic_molecule, dimensionality, dimensions, lattice_parameter,
                                  orientation, reference_orientation, excitons)

    elif order == 'disordered':
        return get_disordered_system(conditions, generic_molecule, dimensionality, dimensions, num_molecules,
                                     orientation, reference_orientation, excitons)

    else:
        print("A valid 'order' parameter is needed")


def get_disordered_system(conditions,  # External conditions of the system such as temperature
                          generic_molecule,
                          dimensionality=2,
                          dimensions=[10, 10, 0],
                          number_molecules=1000,
                          orientation='parallel',
                          reference_orientation=[1,0.0],
                          excitons={'s_1': ['centre']},):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensionality: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = num_dimensions
    :param number_molecules: Number of molecules in the system. Only needed for a disordered system.
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """
    capacity = get_capacity(dimensionality, dimensions, generic_molecule)
    if number_molecules > capacity:
        print('Only %3d molecules could be fitted' % capacity)

    if orientation is 'antiparallel':
        print('Not antiparallel orientation considered for an amorphous material')
        return

    molecules = []
    molecule_count = 0
    while molecule_count <= number_molecules:

        too_close = True
        while too_close is True:
            coordinates = []
            for dimension in dimensions:
                x = dimension*np.random.random() - dimension/2
                # We want the distribution center at (0,0,0)
                coordinates.append(x)

            too_close = distance_checking(coordinates, molecules)

        orientation_vector = get_orientation(orientation, reference_orientation, dimensionality, pointing=1)
        molecule = generic_molecule.copy()
        molecule.initialize_coordinates(coordinates)
        molecule.initialize_orientation(orientation_vector)
        molecules.append(molecule)
        molecule_count += 1


        if molecule_count == capacity:
            break

    tolerance = (dimensions[0]+dimensions[1]+dimensions[2]) / (3*30)
    centre_indexes = excited_system(molecules, excitons, tolerance)
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions, 'centres': centre_indexes}

    return system


def get_ordered_system(conditions,
                       generic_molecule,
                       dimensionality=2,
                       dimensions=[10, 10, 0],
                       lattice_parameter=0.1,
                       orientation='parallel',
                       reference_orientation=[1, 0, 0],
                       excitons={'s_1': ['centre']}, ):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensionality: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = 3. By default = [10, 10, 0]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if dimensionality == 1:
        return get_1d_ordered_system(conditions, generic_molecule, dimensionality, dimensions, lattice_parameter,
                                     orientation, reference_orientation, excitons)

    elif dimensionality == 2:
        return get_2d_ordered_system(conditions, generic_molecule, dimensionality, dimensions, lattice_parameter,
                                     orientation, reference_orientation, excitons)

    elif dimensionality == 3:
        return get_3d_ordered_system(conditions, generic_molecule, dimensionality, dimensions, lattice_parameter,
                                     orientation, reference_orientation, excitons)

    else:
        print('A number of dimensions must be provided')


def get_1d_ordered_system(conditions,
                          generic_molecule,
                          dimensionality=1,
                          dimensions=[10, 0, 0],
                          lattice_parameter=0.1,
                          orientation='parallel',
                          reference_orientation=[1, 0, 0],
                          excitons={'s_1': ['centre']},):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param generic_molecule: #veure escrit a DIPC
    :param dimensionality: number of dimensions of the system (=1)
    :param dimensions: Dimensions of the system. A list with len = 3. By default = [10, 0, 0]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param orientation: String parameter. Indicates whether the orientation of the molecules is parallel,
    random or antiparallel (each molecule finds its 1sts neighbours antiparallely orientated)
    :param reference_orientation: 3 length string (vector), indicates the privileged orientation.
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if check_lattice(lattice_parameter, generic_molecule) is False:
        print('Lattice parameter smaller than molecular characteristic length')
        return

    molecules = []
    x_max = dimensions[0]/2
    pointing = 1
    symmetry = get_symmetry(orientation)
    # We want the distribution center at 0
    for x in np.arange(-x_max, x_max, lattice_parameter):
        coordinates = [x, 0, 0]
        orientation_vector = get_orientation(orientation, reference_orientation, dimensionality, pointing)

        molecule = generic_molecule.copy()
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
                          dimensionality,
                          dimensions=[10, 10, 0],
                          lattice_parameter=0.1,
                          orientation = 'parallel',
                          reference_orientation = [1, 0 ,0],
                          excitons={'s_1': ['centre']},):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param generic_molecule: #veure escrit a DIPC
    :param dimensionality: number of dimensions of the system (=2)
    :param dimensions: Dimensions of the system. A list with len = 3. By default [10, 10, 0]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param orientation: String parameter. Indicates whether the orientation of the molecules is parallel,
    random or antiparallel (each molecule finds its 1sts neighbours antiparallely orientated)
    :param reference_orientation: 3 length string (vector), indicates the privileged orientation.
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if check_lattice(lattice_parameter, generic_molecule) is False:
        print('Lattice parameter smaller than molecular characteristic length')
        return

    molecules = []
    x_max = dimensions[0]/2
    y_max = dimensions[1]/2
    # We want the center of the distribution at (0,0)
    symmetry = get_symmetry(orientation)
    x_count = 0
    for x in np.arange(-x_max, x_max, lattice_parameter):
        pointing = (-1)**x_count
        for y in np.arange(-y_max, y_max, lattice_parameter):
            coordinates = [x, y, 0]
            orientation_vector = get_orientation(orientation, reference_orientation, dimensionality, pointing)

            molecule = generic_molecule.copy()
            molecule.initialize_coordinates(coordinates)
            molecule.innitialize_orientation(orientation_vector)
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
                          dimensionality=3,
                          dimensions=[10, 10, 10],
                          lattice_parameter=0.1,
                          orientation='parallel',
                          reference_orientation=[1, 0, 0],
                          excitons={'s_1': ['centre']},):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensions: Dimensions of the system. A list with len = num_dimensions. By deafult [10, 10, 10]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if check_lattice(lattice_parameter, generic_molecule) is False:
        print('Lattice parameter smaller than molecular characteristic length')
        return

    molecules = []
    x_max = dimensions[0]/2
    y_max = dimensions[1]/2
    z_max = dimensions[2]/2
    step = lattice_parameter
    # We want the center of the distribution at (0,0,0)

    x_count = 0
    symmetry = get_symmetry(orientation)
    for x in np.arange(-x_max, x_max, step):
        x_pointing = (-1)**x_count
        y_count = 0
        for y in np.arange(-y_max, y_max, step):
            pointing = x_pointing * (-1)**y_count
            for z in np.arange(-z_max, z_max, step):
                coordinates = [x, y, z]
                orientation_vector = get_orientation(orientation, reference_orientation, dimensionality, pointing)

                molecule = generic_molecule.copy()
                molecule.initialize_coordinates(coordinates)
                molecule.innitialize_orientation(orientation_vector)
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
def distance_checking(coordinates, molecules):
    """
    :param coordinates: coordinates of the studied molecule
    :param molecules: list of the molecules already defined
    :return: Boolean. Indicates if the new molecule is too close to some of the already defined.
    """
    coordinates = np.array(coordinates)

    for molecule in molecules:
        molecular_coordinates = np.array(molecule.coordinates)

        if distance.euclidean(coordinates, molecular_coordinates) < molecule.molecular_volume:
            return True

    return False


def get_capacity(num_dimensions, dimensions, generic_molecule):
    """
    :param num_dimensions: Dimensionality of the system (1, 2, 3)
    :param dimensions: Size of the system
    :param generic_molecule: instance of class molecule with all its natural parameters defined
    :return: Integer. Maximus number of molecules that could fit in the system.
    """
    elemental_site = generic_molecule.characteristic_length ** num_dimensions

    total_volume = 1
    for dimension in dimensions:
        total_volume = total_volume*dimension

    return int(total_volume/elemental_site)


def check_lattice(lattice_parameter, generic_molecule):
    """
    :param lattice_parameter: Lattice parameter
    :param generic_molecule: instance of class molecule with all its natural parameters defined
    :return: Boolean. Checks if the lattice parameter chosen is smaller than the molecular characteristic length.
    """
    if lattice_parameter <  generic_molecule.characteristic_length:
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


def get_symmetry(orientation):
    if orientation is 'parallel':
        return [1, 2]

    if orientation is 'antiparallel':
        return [-1, 1]

    else:
        return [0, 0]


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
