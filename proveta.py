from systems.initialize_system import get_homogeneous_system
from update_functions.update_file import update_system

conditions = {'temperature': 273.15}
system = get_homogeneous_system(conditions)

print(len(system['molecules']))
print(conditions)
