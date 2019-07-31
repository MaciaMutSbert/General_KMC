import numpy as np
from conversion_functions import from_ev_to_au, from_ns_to_au


class Molecule:

    def __init__(self,
                 states_energies,
                 state,
                 transition_moment,
                 characteristic_length=10**(-8),
                 coordinates=[0, 0, 0],
                 orientation=[1, 0, 0]):
        """
        :param states_energies: dictionary {'state': energy}
        :param state: sting of the name of the state
        The name of the state should coincide with some key of the dictionary in order to identify the state with
        its energy.
        Names of 'state' would be: g_s (ground state), s_1 (first singlet), t_1 (first triplet), etc.
        Energies should be given with eV.

        :param transition_moment: Dipole transition moment vector (3d). The vector is given in respect to the RS
        of the molecule. So for all molecules of a same type if will be equal.
        This dipole moment is given in atomic units..

        :param characteristic_length: We consider a finite size molecule. The simplified shape of the molecule
        is longitudinal, squared or cubic and is defined with this characteristic length. Units: nm

        :param coordinates: 3d vector. Gives the position of the molecule in the system (in general the 0 position
        will coincide with the center of the distribution). Units: nm. If the system has less than 3 dimensions,
        the extra coordinates will be taken as 0.
        :param orientation: 3d unit vector. Gives the orientation of the molecule in the global reference system.
        """

        self.state_energies = states_energies
        self.state = state
        self.transition_moment = transition_moment
        self.characteristic_length = characteristic_length
        self.coordinates = coordinates
        self.orientation = orientation

    def initialize_coordinates(self, coordinate_list):
        """
        :param coordinate_list: List [x, y, z] with the coordinates of the molecule. Units: nm
        Changes self.coordinates to this new position. Format: numpy array.
        """
        self.coordinates = np.array(coordinate_list)

    def molecular_coordinates(self):
        """
        :return: Array with the molecular coordinates.
        """
        return self.coordinates

    def initialize_orientation(self, orientation):
        """
        :param orientation: list with the coordinates of the orientation vector
        Changes self.orientation to this new orientation. Format: numpy array
        """
        self.orientation = np.array(orientation)

    def molecular_orientation(self):
        """
        :return: Array with the molecular orientation
        """
        return self.orientation


    """
    VEURE AL DIPC QUINES FUNCIONS MÉS TENIM
    """

    def desexcitation_energies(self):
        """
        D'ENTRADA NO SERÀ NECESSÀRIA.
        Given an electronic state, calculates the possible desexcitation energy. Generates and sorts
        a list with the energies, then calculates the possible desexcitation energies (the energy difference
        between the given state and the less energetic states).
        :return: Dictionary with the decay processes as key, e.g. 'State1_to_State0', and the energy as argument
        """
        desexcitations = {}

        for state_key in self.state_energies:
            if self.state_energies[self.state] > self.state_energies[state_key]:
                decay_process = 'from_'+self.state+'_to_'+state_key
                energy_gap = self.state_energies[self.state] - self.state_energies[state_key]
                desexcitations[decay_process] = energy_gap

        return desexcitations

    def decay_rates(self):
        """
        :return: Dictionary with the decay processes as key and the decay rates as arguments.
        """
        decay_rates = {}
        if self.state == 's_1':
            "We consider only the radiative decay to singlet 0 state"
            decay_process = 'Singlet_radiative_decay'
            desexcitation_energy = self.state_energies[self.state] - self.state_energies['s_0']
            desexcitation_energy = from_ev_to_au(desexcitation_energy, 'direct')            # energy in atomic units
            u = np.linalg.norm(self.transition_moment)                                      # transition moment norm.

            c = 137         # light speed in atomic units

            rate = 4 * desexcitation_energy**3 * u**2 /(3 * c**3)
            decay_rates[decay_process] = from_ns_to_au(rate, 'inverse')

        return decay_rates

    def get_transition_moment(self, reference_orientation):
        """
        This method computes a basis transformation in order to get the transition dipole moment in a global
        reference coordinate system. This system is described by the reference_orientation vector.
        When the molecule is orientated in this way the transition_moment in the RS system of the molecule
        and in the global RS coincides. Else, we can say that the molecule is rotated and the rotation angle
        is given by the inner product between the reference_orientation and molecule_orientation vectors
        :param reference_orientation: Orientation in which the t.d.m coincides in both the molecular and global
        reference system. UNIT  3D VECTOR
        :return: The t.d.m in the global reference system
        """

        # compute the director cosinus of the rotation (inner product)
        cos_director = np.dot(reference_orientation, self.orientation)

        # construct the rotation matrix
        sin_director = np.sqrt(1-cos_director**2)
        rotation_matrix = np.array([[cos_director, -sin_director],
                                    [sin_director, cos_director]])

        # basis transform (matrix product)
        return np.dot(rotation_matrix, self.transition_moment)
