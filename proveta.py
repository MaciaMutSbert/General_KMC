from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system, check_finish
from result_analysis import initial_position, final_position, exciton_shift, x_y_splitter
import matplotlib.pyplot as plt

conditions = {'temperature': 273.15}

initial_positions = []
final_positions = []
time_list = []

for i in range(100):
    system = get_homogeneous_system(conditions)
    total_time = 0
    path_list = []
    finished = False

    while finished is False:
        path, time = update_system(system)
        total_time += time
        path_list.append(path)

        if check_finish(path_list) is True:
            break

    initial_positions.append(initial_position(path_list, system))
    final_positions.append(final_position(path_list, system))
    time_list.append(total_time)

print(initial_positions)
exciton_distance = exciton_shift(initial_positions, final_positions)
final_x_list, final_y_list = x_y_splitter(final_positions)


"""
                    Results
"""
plt.plot(final_x_list, final_y_list, 'ro', label='Position')
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Final positions')
plt.show()
