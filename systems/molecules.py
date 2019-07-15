class Molecule:

    def __init__(self, coordinates, state=0):
        """
        :param coordinates: Coordinates of the molecule
        :param state: state of the molecule. For now an integer: 0-ground state, 1-excited singlet (...)
        Think about state structure.
        :characteristic_length: Dimensions (length) of the molecule. We consider it finite.
        We'll use atomic units (taking the mass and charge of the electron as unity)
        """
        self.coordinates = coordinates
        self.state = state
        self.characteristic_length = 0.000000001
        self.transition_dipole = 6      # atomic units
        self.singlet_excitation_energy = 2.0      # eV

    def decay_rate(self):         # Static method ???
        """
        :param state: Indicates the excitonic state of the molecule
        :return: A dictionary with the possible decay processes as keys and its rates as arguments.
        (State argument structure may change)
        """
        c = 137      # atomic units
        if self.state == 1:
            decay = (4*(self.singlet_excitation_energy**2)*self.transition_dipole**2) / (3*c**3)
            return {'Singlet_radiative_decay_rate': decay}

        """
        Some type of correction must be included in case that excited_time != 0.
        """
