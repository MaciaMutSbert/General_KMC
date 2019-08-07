import json
import numpy as np
from analysis_functions import get_root_mean_square, initial_final_distance, duration, stadistics


with open('results_file.json', 'r') as read_file:
    data = json.load(read_file)

trajectories = data['trajectories']

"""
Per cada trajectòria hem de calcular el root_mean_square, la distància entre el punt inicial i el final (pel
que fa a longituds de difussió). A més de calcular el temps de vida de l'excitó.
"""
rms_collector = []
exciton_shift_collector = []
exciton_lifetime_collector = []

for trajectory in trajectories:
    rms_collector.append(get_root_mean_square(trajectory))
    exciton_shift_collector.append(initial_final_distance(trajectory))
    exciton_lifetime_collector.append(duration(trajectory))

rms_average, rms_deviation = stadistics(rms_collector)

exciton_shift_average,  exciton_shift_deviation = stadistics(exciton_shift_collector)

exciton_lifetime_average, exciton_lifetime_deviation = stadistics(exciton_lifetime_collector)

output = {'rms': {'average': rms_average, 'deviation': rms_deviation},
          'exciton_shift': {'average': exciton_shift_average, 'deviation': exciton_shift_deviation},
          'exciton_lifetime': {'average': exciton_lifetime_average, 'deviation': exciton_lifetime_deviation}}

with open('analysis_results.json', 'w') as write_file:
    json.dump(output, write_file)


