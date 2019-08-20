import json
import numpy as np
import matplotlib.pyplot as plt
from analysis.diffusion import statistical_diffusion_study, diffusion_parameters
from analysis.theorethical_functions import theoretical_diffusion_values


input_file_name = '1d_simulation_trajectories.json'             # name of the file with the simulation data (.json)

with open(input_file_name, 'r') as read_file:                   # reading of the file (.json)
    simulation_data = json.load(read_file)


#   Information about the simulation conditions (defined molecule and system)
system_information = simulation_data['system_information']
molecule_information = simulation_data['molecules']

#   Data from the simulation: trajectories of the system
trajectories = simulation_data['trajectories']


################################################################################################################

# theoretical diffusion values in a dictionary (only computed if the system type allows it)
diffusion_theoretical = theoretical_diffusion_values(system_information, molecule_information)

################################################################################################################

# results obtained from the simulation (last positions of each trajectory)
# offers 3 plots:
#       convergence of the diffusion length and exciton values to the theoretical values
#       histogram of the final distances of the exciton
diffusion_experimental = diffusion_parameters(trajectories, diffusion_theoretical, system_information)

################################################################################################################

# calculation of the diffusion constant from the statistical study of the points of the trajectory
# uses the theoretical value of the lifetime to compute the diffusion length
# plot offered:
#       points (<r²>, <t>) over the trajectories at each step
#       lineal regression of these set of values: the slope of the line is the diffusion constant.

diffusion_statistical = statistical_diffusion_study(trajectories, diffusion_theoretical, system_information)

################################################################################################################

diffusion_results_1d = {'theorical': diffusion_theoretical, 'statistical': diffusion_statistical,
                        'experimental': diffusion_experimental}

with open('diffusion_results_1d.json', 'w') as write_file:
    json.dump(diffusion_results_1d, write_file)


r_square = diffusion_statistical_study['mean_square_distances']
lifetime = diffusion_statistical_study['life_times']


"Gràfica r**2 vs t (lifetime), ha de ser una recta de pendent D"
plt.plot(lifetime, r_square, 'ro')
plt.xlabel('Exciton lifetime (ns)')
plt.ylabel('Mean square distances (nm^2)')
plt.title('Diffusion constant, r^2 vs t')
plt.show()

"Histograma amb les posicions finals"
plt.hist(diffusion_experimental['max_distances'])
plt.show()




