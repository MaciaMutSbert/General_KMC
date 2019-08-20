import matplotlib.pyplot as plt
import numpy as np


def final_distances_histogram(squared_distances, bins=15):
    """
    :param squared_distances:
    :param bins:
    :return: Plots an histogram with the final positions of the excitons
    """
    plt.hist(np.sqrt(squared_distances), bins=bins)
    plt.title('Final distances histogram')
    plt.xlabel('Final position (nm)')
    plt.ylabel('# of occurrences')
    plt.show()


def diffusion_length_and_lifetime_convergence(root_mean_square, lifetimes, theoretical_values):
    """
    :param root_mean_square:
    :param lifetimes:
    :param theoretical_values:
    :return: Plots the convergence of the simulation values to the theoretical values with the number of
    trajectories averaged.
    """

    iterations = np.arange(len(root_mean_square))

    # diffusion length convergence:
    theo_diff_length = np.ones((len(iterations), 1)) * theoretical_values['diffusion_length']

    plt.plot(iterations, theo_diff_length, label='theoretical value')
    plt.plot(iterations, root_mean_square, 'ro', label='simulation values')
    plt.xlabel('# of trajectories averaged')
    plt.ylabel('Diffusion length (nm)')
    plt.title('Convergence of diffusion length')
    plt.legend()
    plt.show()

    # exciton lifetime convergence:
    theo_lifetime = np.ones((len(iterations), 1)) * theoretical_values['lifetime']

    plt.plot(iterations, theo_lifetime, label='theoretical value')
    plt.plot(iterations, lifetimes, 'ro', label='simulation values')
    plt.xlabel('# of trajectories averaged')
    plt.ylabel('Lifetime (ns)')
    plt.title('Convergence of exciton lifetime')
    plt.legend()
    plt.show()
