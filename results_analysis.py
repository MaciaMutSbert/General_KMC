import json
import numpy as np
from analysis_functions import stadistical_diffusivity, diffusion_parameters


with open('results_file_2d.json', 'r') as read_file:
    data = json.load(read_file)

trajectories = data['trajectories']

dimensionality = len(data['system_information']['dimensions'])
steps = data['system_information']['steps']

diffusion_constant_list = stadistical_diffusivity(trajectories, steps, dimensionality)

diffusion_constant = np.average(diffusion_constant_list)

theoretical_lifetime = 60           # ns
theoretical_diffusion_constant = 3.25277
theoretical_diffusion_length = np.sqrt(2*dimensionality*theoretical_diffusion_constant*theoretical_lifetime)

diffusion_length = np.sqrt(2*dimensionality*diffusion_constant*theoretical_lifetime)

print('Diffusion constant (average) %.5f' % diffusion_constant)
print('Experimental diffusion length %.5f' % diffusion_length)

print('--------------------------')
print('Theoretical diffusion constant %.5f' % theoretical_diffusion_constant)
print('Theoretical diffusion length %.5f' % theoretical_diffusion_length)


