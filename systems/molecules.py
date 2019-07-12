class Molecule:

    def __init__(self, coordinates, state=0):
        """
        :param coordinates: Coordinates of the molecule
        :param state: state of the molecule. For now an integer: 0-ground state, 1-excited singlet (...)
        :characteristic_length: Dimensions (length) of the molecule. We consider it finite.
        """
        self.coordinates = coordinates
        self.state = state
        self.characteristic_length = 0.000000001

    def decay_rate(self, state):         # Static method ???
        """
        :param state: Indicates the excitonic state of the molecule
        :return: A dictionary with the possible decay processes as keys and its rates as arguments.
        """
        if state == 1:
            decay = 1/3
            return {'Singlet_radiative_decay_rate': decay}
