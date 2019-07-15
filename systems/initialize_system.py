import numpy as np
from scipy.spatial import distance
from systems.molecules import Molecule


"""
def get_system(morphology, excitons):

    if morphology['Homogeneity'] is True:
        return get_homogeneous_system(morphology, excitons)

    else:
        return get_inhomogeneous_system(morphology, excitons)
"""


def get_homogeneous_system(conditions,
                           num_molecules=0,           # Only if order = Disordered
                           num_dimensions=2,
                           dimensions=[10, 10],
                           lattice_parameter=0.1,   # Only if order = Ordered
                           excitons={'number': 1, 'positions': [4950]},
                           order='Ordered'):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param num_molecules: Number of molecules in the system. Only needed for a disordered system.
    :param num_dimensions: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = num_dimensions
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :param order: String argument. Indicates if we want an ordered system (crystal) or an amorphous system.
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if order == 'Ordered':
        return get_ordered_system(conditions, num_dimensions, dimensions, lattice_parameter, excitons)

    elif order == 'Disordered':
        return get_disordered_system(conditions, num_dimensions, dimensions, num_molecules, excitons)

    else:
        print('A valid Order parameter is needed')


def get_disordered_system(conditions,               # External conditions of the system such as temperature
                          num_dimensions=2,
                          dimensions=[10, 10],
                          number_molecules=1000,
                          excitons={'number': 1, 'positions': [5000]}):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param num_dimensions: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = num_dimensions
    :param number_molecules: Number of molecules in the system. Only needed for a disordered system.
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """
    capacity = get_capacity(num_dimensions, dimensions)
    if number_molecules > capacity:
        print('Only %1d molecules could be fitted' % capacity)

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

        molecules.append(Molecule(coordinates=coordinates))
        molecule_count += 1

        if molecule_count == capacity:
            break

    molecules = excited_system(molecules, excitons)
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions}

    return system


def get_ordered_system(conditions,
                       num_dimensions=2,
                       dimensions=[10, 10],
                       lattice_parameter=0.1,
                       excitons={'number': 1, 'positions': [4950]}):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param num_dimensions: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = num_dimensions. By default = [10, 10]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if num_dimensions == 1:
        return get_1d_ordered_system(conditions, dimensions, lattice_parameter, excitons)

    elif num_dimensions == 2:
        return get_2d_ordered_system(conditions, dimensions, lattice_parameter, excitons)

    elif num_dimensions == 3:
        return get_3d_ordered_system(conditions, dimensions, lattice_parameter, excitons)

    else:
        print('A number of dimensions must be provided')


def get_1d_ordered_system(conditions,
                          dimension=10,
                          lattice_parameter=0.1,
                          excitons={'number': 1, 'positions': [50]}):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensions: Dimensions of the system. A list with len = num_dimensions. By default = 10
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if check_lattice(lattice_parameter) is False:
        print('Lattice parameter smaller than molecular characteristic length')
        return

    molecules = []
    x_max = dimension/2
    # We want the distribution center at 0
    for x in np.arange(-x_max, x_max, lattice_parameter):
        molecules.append(Molecule([x]))

    molecules = excited_system(molecules, excitons)
    conditions['lattice_parameter'] = lattice_parameter
    conditions['dimensions'] = dimension
    system = {'molecules': molecules, 'conditions': conditions}

    return system


def get_2d_ordered_system(conditions,
                          dimensions=[10, 10],
                          lattice_parameter=0.1,
                          excitons={'number': 1, 'positions': [4950]}):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensions: Dimensions of the system. A list with len = num_dimensions. By default [10, 10]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if check_lattice(lattice_parameter) is False:
        print('Lattice parameter smaller than molecular characteristic length')
        return

    molecules = []
    x_max = dimensions[0]/2
    y_max = dimensions[1]/2
    # We want the center of the distribution at (0,0)
    for x in np.arange(-x_max, x_max, lattice_parameter):
        for y in np.arange(-y_max, y_max, lattice_parameter):
            molecules.append(Molecule([x, y]))

    molecules = excited_system(molecules, excitons)
    conditions['lattice_parameter'] = lattice_parameter
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions}

    return system


def get_3d_ordered_system(conditions,
                          dimensions=[10, 10, 0],
                          lattice_parameter=0.1,
                          excitons={'number': 1, 'positions': [4950]}):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensions: Dimensions of the system. A list with len = num_dimensions. By deafult [10, 10, 10]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if check_lattice(lattice_parameter) is False:
        print('Lattice parameter smaller than molecular characteristic length')
        return

    molecules = []
    x_max = dimensions[0]/2
    y_max = dimensions[1]/2
    z_max = dimensions[2]/2
    step = lattice_parameter
    # We want the center of the distribution at (0,0,0)
    for x in np.arange(-x_max, x_max, step):
        for y in np.arange(-y_max, y_max, step):
            for z in np.arange(-z_max, z_max, step):
                molecules.append(Molecule([x, y, z]))

    molecules = excited_system(molecules, excitons)
    conditions['lattice_parameter'] = lattice_parameter
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions}

    return system


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


def excited_system(molecules, excitons):
    """
    :param molecules: List of the defined molecules
    :param excitons: Information about the desired excitation
    :return: Excited system
    It could be defined in the file where all processes and excitations will be defined (??)
    """

    if excitons['positions'] == 'random':
        for i in range(excitons['number']):
            molecules[np.random.randint(0, len(molecules))].state = 1

    else:
        for position in excitons['positions']:
            molecules[position].state = 1

    return molecules


def get_capacity(num_dimensions, dimensions):
    """
    :param num_dimensions: Dimensionality of the system (1, 2, 3)
    :param dimensions: Size of the system
    :return: Integer. Maximus number of molecules that could fit in the system.
    """
    molecule = Molecule(0, 0)
    elemental_site = molecule.characteristic_length ** num_dimensions

    total_volume = 1
    for dimension in dimensions:
        total_volume = total_volume*dimension

    return int(total_volume/elemental_site)


def check_lattice(lattice_parameter):
    """
    :param lattice_parameter: Lattice parameter
    :return: Boolean. Checks if the lattice parameter chosen is smaller than the molecular characteristic length.
    """
    molecule = Molecule(0, 0)
    if lattice_parameter < molecule.characteristic_length:
        return False
    else:
        return True
