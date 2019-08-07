from proveta import diffusion_length, diffusion_length_time_residence, final_positions, exciton_shift, ld_rms,\
    ld_rms_average, life_time
from result_analysis import x_y_splitter
import matplotlib.pyplot as plt
import numpy as np


df_deviation = np.float(np.std(exciton_shift))
df_residence_time_deviation = np.float(np.std(ld_rms_average))
lt_deviation = np.float(np.std(np.array(life_time)))


print('Average distance: %.5f +- %.5f' % (diffusion_length[-1], df_deviation))
print('Average distance (computed with the residence time): %.5f +- %.5f' % (ld_rms_average[-1],
                                                                             df_residence_time_deviation))

print('Average time %.5f +- %.5f' % (life_time[-1], lt_deviation))

final_x_list, final_y_list = x_y_splitter(final_positions)

plt.plot(final_x_list, final_y_list, 'bo')
plt.xlabel('final x positions')
plt.ylabel('final y positions')
plt.title('Final positions')
plt.show()

iterations = np.arange(0, 3000, 1)

plt.plot(iterations, diffusion_length, 'ro')
plt.xlabel('Number of iterations')
plt.ylabel('Diffusion length')
plt.show()

plt.hist(exciton_shift, bins=30)
plt.title('histograma amb les distàncies finals')
plt.show()

plt.plot(iterations, ld_rms_average, 'ro')
plt.xlabel('Number of iterations')
plt.ylabel('Diffusion length. Ld_teo = 24,7274')
plt.show()

plt.plot(iterations, life_time, 'ro')
plt.xlabel('Number of iterations, Ld_teo = 24,7274')
plt.ylabel('Life_time')
plt.show()

