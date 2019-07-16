class Molecule:

    def __init__(self, coordinates,
                 transition_dipole,
                 singlet_excitation_energy,
                 molecule_type,
                 state=0):
        """
        :param coordinates: position of the molecule in nm
        :param transition_dipole in atomic units
        :param singlet_excitation_energy in eV
        :param state: integer. Indicates the excited state of the molecule (1-first singlet...)
        """
        self.coordinates = coordinates
        self.state = state
        self.type = molecule_type
        self.characteristic_length = 0.00000001
        self.transition_dipole = transition_dipole      # atomic units
        self.singlet_excitation_energy = singlet_excitation_energy              # eV

    def decay_rate(self):         # Static method ???
        """
        :return: A dictionary with the possible decay processes as keys and its rates as arguments.
        (State argument structure may change)
        """
        c = 137      # light speed in atomic units
        if self.state == 1:
            e_s = self.singlet_excitation_energy / 27.211       # atomic units

            decay = (4 * e_s**3 * self.transition_dipole**2) / (3 * c**3)
            return {'Singlet_decay_rate': decay}

