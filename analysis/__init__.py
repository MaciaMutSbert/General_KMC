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
    if len(system['centres']) == 0:
        # if there is not any excited molecule saves exciton positions as None.
        exciton_positions = None

    else:
        exciton_positions = []             # list with all the positions of the excited molecules (simultaneous)
        for centre in system['centres']:

            exciton_coordinates =(system['molecules'][centre].molecular_coordinates())       # cartesian coordinates
            excited_state = system['molecules'][centre].electronic_state()                  # electronic state

            # a tetra vector is built (x,y,z, state)
            exciton_point = [exciton_coordinates, excited_state]
            exciton_positions.append(exciton_point)

    trajectory['positions'].append(exciton_positions)

    # process occurred:
    trajectory['process'].append(change_step['process'])

    # No return function: only updates trajectory.


