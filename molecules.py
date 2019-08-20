import numpy as np
from conversion_functions import from_ev_to_au, from_ns_to_au


"""
COMENTARI INICIAL: Els estats (fonamental i excitats) s'indicaran per mitjà d'etiquetes tipus string. Els estats possibles 
i les seves respectives referències hauran de venir indicades a l'inici del programa, e.g, quan s'inicialitzi la 
molècula genèrica. Aquestes referències s'usaran per caracteritzar la variable estat i per simplicitat i eficiència
s'usaran com a claus dels diccionaris que contenguin informació sobre els excitons. A saber: state_energies, 
relaxation_energies, exciton_energies, (REFERENTS A RATES...).

INICIALITZACIÓ DE LA MOLÈCULA GENÈRICA:
    -state_energies: diccionari amb les energies d'excitació de cada estat de la molècula:
                :key: energia de l'estat excitat
                :argument:  energia de l'estat
    -relaxation_energies: diccionari amb les energies de relaxació de cada estat. Per ara en prenem una fitxa i què depèn
    sols de cada estat.
                :key: energia de l'estat excitat
                :argument:  energia de relaxació de l'estat
    -transition_moment: moment de transició dipolar de la molècula. Donat com un vector (llista de 3 elements) en el 
    sistema de referència de la molècula.
    
Per defecte venen donats:
    -characteristic_length: defineix les dimensions finites de la molècula. Aquesta s'aproxima com una línia, quadrat o 
    cub i aquest paràmetre en defineix el costat. Alternativament podriem fer una aproximació esfèrica i que en sigui 
    el radi.
    -coordinates. D'entrada la suposam en l'origen. 
    -orientation. D'entrada la suposam orientada segons [1, 0, 0]. Aquest vector ve donat en un SR extern que 
    anomenarem global.
    Aquests dos darrers paràmetres no són estrictament necessaris per estudiar la naturalesa del tipus de molècula.
    La classe inclou 2 mètodes per cada un d'aquests 2 darrers paràmetres. Un per inicialitzar-los (defineix de manera
    la posició/orientació com un 3-array) i un per cridar-los alhora d'operar (per assegurar que no s'alteren en el procés).
    Noms:
    initialize_coordinates(coordinates)                  initialize_orientation(orientation)
    molecular_coordinates()                              molecular_orientation()

Mètodes de la molècula:
Apart dels 4 ja comentats la classe molècula inclou:
    - electronic_state. Mètode que retorna l'estat electrònic de la molecula
    - get_relaxation_state_energy. Mètode que dóna l'energia de relaxació de l'estat en què es troba la molècula
    - change_state(new_state): Canvia l'estat de la molècula pel nou donat.
    - decay_rates: mètode que retorna un diccionari amb els possibles rates de decaïment {'decay process': rate}
    - get_transition_moment(reference_orientation). Necessita com argument un vector de referència. L'orientació de la
        molècula en la qual el moment de transició en el SR de la molècula i en el SR global coindideixen. 
        Aleshores, donada aquesta referència i l'orientació de la molècula, aquest mètode fa un canvi de base i retorna
        el moment de transició dipolar en el SR global.
"""


class Molecule:

    def __init__(self,
                 state_energies,
                 reorganization_energies,
                 transition_moment,
                 dipole_moment_direction=[1, 0, 0],                 # unity vector
                 state='gs',
                 characteristic_length=10**(-8),
                 coordinates=[0, 0, 0],
                 orientation=[1, 0, 0]):                            # unity vector
        """
        :param states_energies: dictionary {'state': energy}
        :param state: sting of the name of the state
        The name of the state should coincide with some key of the dictionary in order to identify the state with
        its energy.
        :param reorganization_energies: dictionary {'state': relaxation energy of the state}
        Names of 'state' would be: g_s (ground state), s_1 (first singlet), t_1 (first triplet), etc.
        Energies should be given with eV.


        :param transition_moment: Dipole transition moment vector (3d). The vector is given in respect to the RS
        of the molecule. So for all molecules of a same type if will be equal.
        This dipole moment is given in atomic units.

        :param dipole_moment_direction: Gives a referene direction This reference direction is the molecular orientation
        in which the transition dipole moment in the molecular reference system coincides with the itself in the
        global reference system. Given by default as [1,0,0]

        :param characteristic_length: We consider a finite size molecule. The simplified shape of the molecule
        is longitudinal, squared or cubic and is defined with this characteristic length. Units: nm

        :param coordinates: 3d vector. Gives the position of the molecule in the system (in general the 0 position
        will coincide with the center of the distribution). Units: nm. If the system has less than 3 dimensions,
        the extra coordinates will be taken as 0.
        :param orientation: 3d unit vector. Gives the orientation of the molecule in the global reference system.
        """

        self.state_energies = state_energies
        self.state = state
        self.reorganization_energies = reorganization_energies
        self.transition_moment = np.array(transition_moment)
        self.dipole_moment_direction = np.array(dipole_moment_direction)        # unity vector
        self.characteristic_length = characteristic_length
        self.coordinates = np.array(coordinates)
        self.orientation = np.array(orientation)                                # unity vector

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

    def get_reorganization_state_energy(self):
        return self.reorganization_energies[self.state]

    def electronic_state(self):
        """
        :return: the electronic state of the molecule
        """
        return self.state

    def set_state(self, new_state):
        """
        :param new_state:
        No return method. Only changes the molecular state when the exciton is transferred.
        """
        self.state = new_state

    def desexcitation_energies(self):
        """
        IS NOT USED (19/08/2019).
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
        :return: A list of two elements: list of the possible decay processes and another with the respective rates
        for a given electronic state.

        More if(s) entrances shall be added if more electronic states are considered.
        """

        decay_rates = {}

        if self.state == 's1':

            desexcitation_energy = self.state_energies[self.state] - self.state_energies['gs']      # energy in eV
            desexcitation_energy = from_ev_to_au(desexcitation_energy, 'direct')                    # energy in a.u.

            u = np.linalg.norm(self.transition_moment)              # transition moment norm.
            c = 137                                                 # light speed in atomic units

            rate = 4 * desexcitation_energy**3 * u**2 / (3 * c**3)
            decay_process = 'Singlet_radiative_decay'
            # for a first singlet state only radiative decay is considered.

            decay_rates[decay_process] = from_ns_to_au(rate, 'direct')

        return decay_rates

    def get_transition_moment(self):
        """
        This method computes a basis transformation in order to get the transition dipole moment in a global
        reference coordinate system. This system is described by the dipole_moment_direction
        When the molecule is orientated in this way the transition_moment in the RS system of the molecule
        and in the global RS coincides. Else, we can say that the molecule is rotated and the rotation angle
        is given by the inner product between the reference_orientation and molecule_orientation vectors
        :return: The t.d.m in the global reference system
        """

        # compute the director cosinus of the rotation (inner product)
        cos_director = np.dot(self.dipole_moment_direction, self.orientation)

        # construct the rotation matrix
        sin_director = np.sqrt(1-cos_director**2)

        rotation_matrix = np.array([[cos_director, -sin_director, 0],
                                    [sin_director, cos_director, 0],
                                    [0,                 0,       0]])

        # basis transform (matrix product)
        return np.dot(rotation_matrix, self.transition_moment)


