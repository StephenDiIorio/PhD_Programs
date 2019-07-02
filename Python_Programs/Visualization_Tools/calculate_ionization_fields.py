from matplotlib import use
use('Agg')

import glob
import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
from matplotlib import rc
from matplotlib.ticker import FuncFormatter

from PlottingTools import get_si_prefix

# rc('mathtext', default='regular')


def ion_intens(energy, ion_charge):
    # leave ion_charge as at least 1 to get intensity to fully ionize
    # can think of this as the number of electrons left to ionize
    return 4E9 * (energy**4 / ion_charge**2) * 1E4


def main():
    title = 'Nitrogen and Oxygen Ionization Intensities'
    N_energies = [14.53414, 29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 667.046]
    O_energies = [13.61806, 35.11730, 54.9355, 77.41353, 113.8990, 138.1197, 739.29, 871.4101]

    N_intens = []
    O_intens = []

    for i, val in enumerate(N_energies):
        N_intens.append(ion_intens(val, i + 1))
    for i, val in enumerate(O_energies):
        O_intens.append(ion_intens(val, i + 1))

    fig = plt.figure()
    plt.title(title, y=1.12)

    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    ax = fig.add_subplot(111)
    ax.set_xlabel(r'Electric Field $(V/m)$')

    left_lim, right_lim = [1e10, 1e14]
    ax.set_xscale('log')
    ax.set_xlim(left_lim, right_lim)
    scale_left_lim = np.power(left_lim, 2) * sc.speed_of_light * sc.epsilon_0 / 2.0
    scale_right_lim = np.power(right_lim, 2) * sc.speed_of_light * sc.epsilon_0 / 2.0

    ax2 = ax.twiny()
    ax2.set_xscale('log')
    ax2.set_xlim(scale_left_lim, scale_right_lim)
    ax2.set_xlabel(r'Intensity $(W/m^2)$')
    for val in O_intens:
        l1 = ax2.axvline(x=val, c='r', label='O Intensities')
    for val in N_intens:
        l2 = ax2.axvline(x=val, c='b', label='N Intensities')

    ls = [l1, l2]
    labels = [l.get_label() for l in ls]
    ax.legend(ls, labels, loc=0)

    plt.tight_layout()
    plt.savefig(title + '.png', dpi=600)
    return


if __name__ == "__main__":
    main()
