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
                           num_molecules=0,  # Only if order = Disordered
                           dimensionality=2,
                           dimensions=[10, 10, 0],
                           lattice_parameter=0.1,  # Only if order = Ordered
                           excitons={'number': 1, 'positions': [85]},
                           order='Ordered'):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param num_molecules: Number of molecules in the system. Only needed for a disordered system.
    :param dimensionality: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = num_dimensions
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :param order: String argument. Indicates if we want an ordered system (crystal) or an amorphous system.
    :return: A dictionary with a list of molecules, an updated dictionary with the physical conditions
    and a list with the indexes of the excited molecules.
    """

    if order == 'Ordered':
        return get_ordered_system(conditions, dimensionality, dimensions, lattice_parameter, excitons)

    elif order == 'Disordered':
        return get_disordered_system(conditions, dimensionality, dimensions, num_molecules, excitons)

    else:
        print('A valid Order parameter is needed')


def get_disordered_system(conditions,  # External conditions of the system such as temperature
                          dimensionality=2,
                          dimensions=[10, 10, 0],
                          number_molecules=1000,
                          excitons={'number': 1, 'positions': [500]}):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensionality: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = num_dimensions
    :param number_molecules: Number of molecules in the system. Only needed for a disordered system.
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """
    capacity = get_capacity(dimensionality, dimensions)
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

        e_s = conditions['singlet_energy']
        u = conditions['transition_dipole']
        molecules.append(Molecule(coordinates=coordinates, transition_dipole=u, singlet_excitation_energy=e_s))
        molecule_count += 1

        if molecule_count == capacity:
            break

    centre_indexes = excited_system(molecules, excitons)
    conditions['dimensions'] = dimensions
    system = {'molecules': molecules, 'conditions': conditions, 'centres': centre_indexes}

    return system


def get_ordered_system(conditions,
                       dimensionality=2,
                       dimensions=[10, 10, 0],
                       lattice_parameter=0.1,
                       excitons={'number': 1, 'positions': [85]}):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensionality: Dimensionality of the system (1, 2, 3). Default = 2
    :param dimensions: Dimensions of the system. A list with len = 3. By default = [10, 10, 0]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if dimensionality == 1:
        return get_1d_ordered_system(conditions, dimensions, lattice_parameter, excitons)

    elif dimensionality == 2:
        return get_2d_ordered_system(conditions, dimensions, lattice_parameter, excitons)

    elif dimensionality == 3:
        return get_3d_ordered_system(conditions, dimensions, lattice_parameter, excitons)

    else:
        print('A number of dimensions must be provided')


def get_1d_ordered_system(conditions,
                          generic_molecule,
                          dimensionality=1,
                          dimensions=[10, 0, 0],
                          lattice_parameter=0.1,
                          orientation='parallel',
                          reference_orientation=[1, 0, 0],
                          excitons={'number': 1, 'positions': [50]}):
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

    if check_lattice(lattice_parameter, conditions) is False:
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

    centre_indexes = excited_system(molecules, excitons)
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
                          excitons={'number': 1, 'positions': [85]}):
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

    if check_lattice(lattice_parameter, conditions) is False:
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

    centre_indexes = excited_system(molecules, excitons)
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
                          excitons={'number': 1, 'positions': [85]}):
    """
    :param conditions: A dictionary with the pysical conditions of the problem such as temperature.
    :param dimensions: Dimensions of the system. A list with len = num_dimensions. By deafult [10, 10, 10]
    :param lattice_parameter: Parameter of the lattice we want to construct. Default = 0.1
    :param excitons: Number of excitons in the system. Default  1 in the 5000th molecule of the list (...)
    :return: A dictionary with a list of molecules and updated dictionary with the physical conditions.
    """

    if check_lattice(lattice_parameter, conditions) is False:
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


    centre_indexes = excited_system(molecules, excitons)
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


def excited_system(molecules, excitons):
    """
    :param molecules: List of the defined molecules
    :param excitons: Information about the desired excitation
    :return: Excited system, list with the index of the excited molecules
    It could be defined in the file where all processes and excitations will be defined (??)
    """
    centre_list = []

    if excitons['positions'] == 'random':
        for i in range(excitons['number']):
            centre = np.random.randint(0, len(molecules))
            molecules[centre].state = 1
            centre_list.append(centre)
    else:
        for position in excitons['positions']:
            molecules[position].state = 1
            centre_list.append(position)

    return centre_list


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


def check_lattice(lattice_parameter, conditions):
    """
    :param lattice_parameter: Lattice parameter
    :param conditions:
    :return: Boolean. Checks if the lattice parameter chosen is smaller than the molecular characteristic length.
    """
    characteristic_length = conditions['characteristic_length']
    if lattice_parameter < characteristic_length:
        return False
    else:
        return True


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
            return [np.sin(thehta)*np.cos(phi), np.sin(thetha)*np.sin(phi), np.cos(theta)]

        else:
            phi = 2*np.pi() * np.random.rand()
            return [np.cos(phi), np.sin(phi), 0]

    else:
        return reference_orientation * pointing