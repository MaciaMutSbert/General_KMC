import json
import numpy as np
import matplotlib.pyplot as plt

"""
In this script the user can get the plot of l² vs t for different diffusion situations.
As default there is a reference set of values corresponding to the values introduced in reference_simulation in 
example_simulation directory.

Given the exciton lifetime a set of equidistant time points between 0 and the lifetime is
constructed. Then with the diffusion constant we compute a set for l². We plot both sets.
The slope of the line will be D and the final point will be (lifetime, diffusion_length).

As example a comparison between the diffusion of a singlet exciton in a molecular lattice with transition dipole
moment = [1,0,0] and another with [2,0,0] is presented.
"""

# REFERENCE DATA
reference_file = 'diffusion_results_1d.json'

with open(reference_file, 'r') as readfile:
    reference_data = json.load(readfile)

ref_diffusion_constant = reference_data['experimental']['diffusion_constant']
ref_lifetime = reference_data['experimental']['lifetime']

ref_time_set = np.arange(0, ref_lifetime, 0.1)

ref_squared_distance_set = 2 * ref_diffusion_constant * ref_time_set
# the 2 factor stands for the double of the dimensionality

################################################

# NEW DATA. TRANSITION DIPOLE MOMENT = [2,0,0]
file_1 = 'diffusion_results_1d_mu_1-1.json'

with open(file_1, 'r') as readfile:
    data_1 = json.load(readfile)

diffusion_constant_1 = data_1['experimental']['diffusion_constant']
lifetime_1 = data_1['experimental']['lifetime']

time_set_1 = np.arange(0, lifetime_1, 0.1)

squared_distance_set_1 = 2 * diffusion_constant_1 * time_set_1


##########################################################################


# PLOTS l² vs t for both cases
plt.plot(ref_time_set, ref_squared_distance_set, color='b', label='$ mu = 1, D= %.3f $' % ref_diffusion_constant)
plt.plot(time_set_1, squared_distance_set_1, color='r', label='$ mu = 1.1, D= %.3f $' % diffusion_constant_1)
plt.xlabel('$ Time, ns $')
plt.ylabel('$ l^2, nm $')
plt.legend()
plt.show()




