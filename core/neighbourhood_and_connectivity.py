import numpy as np
from scipy.spatial import distance

# neighbourhood is defined in this new file because it surely will be expanded when defining connectivities
# for some systems


def neighbourhood(centre, system, radius=1.05):
    """
    :param centre: Index of an excited Molecule object
    :param system: dictionary with all the information of the system
    :param radius: Effective distance where interaction may be considerable. Default 0.11
    :return: List of indexes of molecules in a neighbourhood of center
    If there is not any neighbours in the defined neighbourhood an alert is printed.
    """

    molecules = system['molecules']
    center_position = np.array(molecules[centre].molecular_coordinates())
    # position of the excited molecule (as a numpy array)

    # some connectivity tricks have been defined for specific cases: 1d, 2d ordered systems (squared systems).

    if system['type'] == '1d_ordered':
        neighbours = [centre-1, centre+1]

        if centre == 0:                      # removes the previous if the center is in the first position of the list
            neighbours.remove(centre-1)
        elif centre == len(molecules)-1:     # removes the next if the center is in the last position of the list
            neighbours.remove(centre+1)

    if system['type'] == '2d_ordered':
        dimensions = system['conditions']['dimensions']
        lat_param = system['conditions']['lattice_parameter']

        neighbours = [centre-dimensions[1], centre-1, centre+1, centre+dimensions[1]]
        # the previous and the next molecule are taken. Besides, the molecules of the adjacent columns need to be taken
        # The number of molecules per column is used to take these last positions

        # if the excited molecule is on the border of the system the extra molecules are removed.
        if center_position[0] == -dimensions[0]*lat_param / 2:                  # **
            neighbours.remove(centre-dimensions[1])
        elif center_position[1] == -dimensions[1]*lat_param / 2:                # **
            neighbours.remove(centre-1)
        elif center_position[0] == dimensions[0] * lat_param / 2 - lat_param:   # **
            neighbours.remove(centre + dimensions[1])
        elif center_position[1] == dimensions[1] * lat_param / 2 - lat_param:   # **
            neighbours.remove(centre+1)
        # **: bordering coordinates for each side
        # the way they are chosen and removed has been tested.

    # for different type systems the neighbours are searched by hand. The program looks for molecules around the centre
    # closer than a radius.
    else:
        neighbours = []
        for i, molecule in enumerate(molecules):
            coordinates = np.array(molecule.coordinates)

            # checks if the distance is smaller than the given radius and bigger than 0 (avoid taking the center as
            # neighbour) If the molecule satisfies both its index is taken.
            if 0 < distance.euclidean(center_position, coordinates) < radius:
                neighbours.append(i)

    if len(neighbours) == 0:
        print('No neighbours found. Check neighbourhood radius')

    return neighbours
