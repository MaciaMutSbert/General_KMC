import json
import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import statistical_diffusivity, diffusion_parameters
from processes.theorethical_functions import theorethical_diffusion_values


input_file_name = '1d_simulation_trajectories_random.json'

with open(input_file_name, 'r') as read_file:
    data = json.load(read_file)

trajectories = data['trajectories']

dimensionality = len(data['system_information']['dimensions'])
trajectory_steps = data['trajectory_steps']

"Theorical values for the diffusion parameters"
system_information = data['system_information']
diffusion_theorical = theorethical_diffusion_values(system_information)

print('Theorical values')
print(diffusion_theorical)

"Dictionary with a list of diffusion parameters computed for every iteration"
diffusion_statistical_study = statistical_diffusivity(trajectories, trajectory_steps, dimensionality)

stad_diffusion_constant = np.average(np.array(diffusion_statistical_study['diffusion_constants']))
stad_diffusion_length = np.sqrt(2* dimensionality* stad_diffusion_constant * diffusion_theorical['life_time'])

diffusion_statistical_values = {'diffusion_constant': stad_diffusion_constant,
                                'diffusion_length': stad_diffusion_length}

print('\n Statistical values:')
print(diffusion_statistical_values)

"Dictionary with the life time and diffusion length using only the last point"
diffusion_experimental_study = diffusion_parameters(trajectories, dimensionality)

print('\n Experimental values')
print(diffusion_experimental_study)


diffusion_results_1d = {'theorical': diffusion_theorical, 'statistical': diffusion_statistical_values,
                        'experimental': diffusion_experimental_study, 'comment': 'Valors teòrics no vàlids'}

with open('diffusion_results_random.json', 'w') as write_file:
    json.dump(diffusion_results_1d, write_file)


r_square = diffusion_statistical_study['mean_square_distances']
lifetime = diffusion_statistical_study['life_times']


"Gràfica r**2 vs t (lifetime), ha de ser una recta de pendent D"
plt.plot(lifetime, r_square, 'ro')
plt.xlabel('Exciton lifetime (ns)')
plt.ylabel('Mean square distances (nm^2)')
plt.title('Diffusion constant, r^2 vs t')
plt.show()






