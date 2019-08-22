import json
import numpy as np


def update_trajectory(trajectory, change_step, time_step, system):
    """
    :param trajectory: dictionary with the system trajectory (time, number of excitons, positions and process occurred)
    :param change_step: process occured during time_step: {donor, process, acceptor}
    :param time_step: duration of the chosen process
    :param system: dictionary with the information of the system
    No return function, only updates the dictionary trajectory
    """

    # time update:
    trajectory['time'].append(trajectory['time'][-1] + time_step)

    # number of excitons update:
    trajectory['n'].append(len(system['centres']))


    # excitons positions update:
    exciton_positions = []  # list with all the positions of the excited molecules (simultaneous)

    if len(system['centres']) == 0:
        # if there is not any excited molecule saves the position of the last excited molecule.

        last_excited = change_step['acceptor']
        exciton_coordinates = list(system['molecules'][last_excited].molecular_coordinates())   # cartesian coordinates
        excited_state = system['molecules'][last_excited].electronic_state()                    # electronic state

        # a tetra vector is built ([x,y,z], state)
        exciton_point = [exciton_coordinates, excited_state]

        # collector of the positions of all excitons transferring.
        exciton_positions.append(exciton_point)

    else:
        for centre in system['centres']:

            exciton_coordinates = list(system['molecules'][centre].molecular_coordinates())      # cartesian coordinates
            excited_state = system['molecules'][centre].electronic_state()                       # electronic state

            # a tetra vector is built ([x,y,z], state)
            exciton_point = [exciton_coordinates, excited_state]

            # collector of the positions of all excitons transferring.
            exciton_positions.append(exciton_point)

    trajectory['positions'].append(exciton_positions)

    # process occurred:
    trajectory['process'].append(change_step['process'])

    # No return function: only updates trajectory.


def merge_json_files(file1, file2):
    """
    :param file1: name (string) of the first file in .json format
    :param file2: the same for the second file
    :return:
        if the headings of both files (system_information and molecule_information) are equal, then
        the trajectories are merged. A single file with the same heading and with the merged trajectories is returned.
        in other case, None is returned and a a warning message is printed.
    """

    # reading of the files
    with open(file1, 'r') as read_file1:
        data_1 = json.load(read_file1)

    with open(file2, 'r') as read_file2:
        data_2 = json.load(read_file2)

    heading_1 = data_1['system_information']
    heading_2 = data_2['system_information']

    if heading_1 == heading_2:
        trajectories_1 = data_1['trajectories']
        trajectories_2 = data_2['trajectories']

        merged_trajectories = trajectories_1 + trajectories_2

        data_1['trajectories'] = merged_trajectories

        # data_1 has in its 'trajectories' entrance the merged information of both files,
        # while keeping the system information (that is equal in both files).
        with open(file1, 'w') as write_file:
            json.dump(data_1, write_file)

        print('Trajectories has been merged in ' + file1)

        return

    else:
        print('Trajectories taken in different conditions.')
        return


def get_l2_t_from_file(diffusion_file):
    """
    :param diffusion_file: json file with diffusion results
    :return: The function takes from the file the experimental values of D and lifetime, It returns a list
    with values of time from 0 to lifetime. With D computes the respective values of l².
    The aim of this function is to get this sets of data to plot several of them and be able to compare the slopes.
    As well, the value of the lifetime and diffusion length will be seen in this plot (the last point of each line).
    """

    with open(diffusion_file, 'r') as readfile:
        data = json.load(readfile)

    diffusion_constant = data['experimental']['diffusion_constant']
    lifetime = data['experimental']['lifetime']

    time_set = np.arange(0, lifetime, 0.1)

    squared_distance_set = (2 * diffusion_constant * time_set)
    # the 2 factor stands for the double of the dimensionality

    changed_parameter = data['changed_parameter']

    return time_set, squared_distance_set, diffusion_constant, changed_parameter


def get_diffusion_vs_mu(diffusion_file):
    """
    :param diffusion_file: json file with diffusion results
    :return: Reads the file and returns:
        D, lifetime, L_d, mu
    """

    with open(diffusion_file, 'r') as readfile:
        data = json.load(readfile)

    diffusion_constant = data['experimental']['diffusion_constant']
    lifetime = data['experimental']['lifetime']
    diffusion_length = data['experimental']['diffusion_length']

    mu = np.linalg.norm(np.array(data['conditions']['transition_moment']))

    return diffusion_constant, lifetime, diffusion_length, mu

