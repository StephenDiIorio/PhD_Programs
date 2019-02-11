import glob

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
from matplotlib import rc, use
from matplotlib.ticker import FuncFormatter

use('Agg')
rc('mathtext', default='regular')


def get_si_prefix(scale):
    scale = abs(scale)
    mult = 1
    sym = ''

    scale = scale * mult
    if scale <= 0:
        pwr = 0
    else:
        pwr = (-np.floor(np.log10(scale)))
    mult = mult * np.power(10.0, pwr)
    if np.rint(pwr) != 0:
        sym = "x10^{%.0f}" % (-pwr) + sym

    return mult, sym


def ion_intens(energy, ion_charge):
    return 4E9 * (energy**4 / ion_charge**2) * 1E4


def main():
    filepath = "Ion_Results/10000fs/2Levels_N"
    files = glob.glob(filepath + "/*.npz")
    # print(files)
    # print(len(files))

    title = 'Nitrogen Ionization Last 2 Ionization Levels(10000fs), w BSI, wo Multiphoton'
    # ion_energies = [13.59844]
    #ion_energies = [24.58741, 54.41778]
    ion_energies = [14.53414, 29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 667.046]

    e_nums = []
    e_field = []
    calc_intens = []

    for f in files:
        data = np.load(f)
        e_nums.append(data['arr_0'][0])
        e_field.append(data['arr_1'][0])
        data.close()

    e_field, e_nums = (list(t) for t in zip(*sorted(zip(e_field, e_nums))))

    for i, val in enumerate(ion_energies):
        calc_intens.append(ion_intens(val, i + 1))

    # print(e_nums)
    # print(e_field)
    # print(calc_intens)

    fig = plt.figure()
    plt.title(title, y=1.12)

    ax = fig.add_subplot(111)
    l1, = ax.semilogx(e_field, e_nums, label='EPOCH Results')
    ax.set_xlabel(r'Electric Field $(V/m)$')
    ymult, ysym = get_si_prefix(np.max(e_nums) - np.min(e_nums))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult)))
    ax.set_ylabel(r'Derived Electron Number Density $(' + ysym + 'm^{-3})$')

    left_lim, right_lim = ax.get_xlim()
    scale_left_lim = np.power(left_lim, 2) * sc.speed_of_light * sc.epsilon_0 / 2.0
    scale_right_lim = np.power(right_lim, 2) * sc.speed_of_light * sc.epsilon_0 / 2.0

    ax2 = ax.twiny()
    ax2.set_xscale('log')
    ax2.set_xlim(scale_left_lim, scale_right_lim)
    ax2.set_xlabel(r'Intensity $(W/m^2)$')
    for val in calc_intens:
        l2 = ax2.axvline(x=val, c='r', label='Calculated Intensities')

    ax3 = ax.twinx()
    ax3.set_ylim(0, len(ion_energies))
    ax3.set_ylabel('Ionization per Atom')

    ls = [l1, l2]
    labels = [l.get_label() for l in ls]
    ax.legend(ls, labels, loc=0)

    plt.tight_layout()
    plt.savefig(title + '.png')
    return


if __name__ == "__main__":
    main()
