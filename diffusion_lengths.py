import numpy as np
from scipy.spatial import distance


def final_position_distance(path_collector, system):
    index_i = path_collector[0]['donor']
    initial_position = system['molecules'][index_i]

    index_f = path_collector[-1]['acceptor']
    final_position = system['molecules'][index_f]

    return distance.euclidean(initial_position, final_position)


def furthest_position_distance(path_collector, system, max_distance):
    max_distance = 0

    for i, path in enumerate(path_collector):
        path_part = path_collector[0:i]
        if final_position_distance(path_part, system) >= max_distance:
            max_distance = final_position_distance(path_collector, system)

        else:
            max_distance = max_distance

    return max_distance


def diffusion_length_time_residence(time_list, lattice_parameter):
    life_time = 1 / 0.01662         # ns

    ld_list = []
    for time in time_list:
        ld = lattice_parameter * np.sqrt(life_time / (2*time))
        ld_list.append(ld)

    return np.average(ld_list)