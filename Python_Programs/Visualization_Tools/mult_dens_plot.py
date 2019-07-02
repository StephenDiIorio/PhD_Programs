from matplotlib import use
use('Agg')

import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

import sdf
from PlottingTools import get_si_prefix, get_var_range


plt.rc('font', size=20)        # controls default text sizes
plt.rc('axes', titlesize=18)   # fontsize of the axes title
plt.rc('axes', labelsize=15)   # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)  # fontsize of the tick labels
plt.rc('ytick', labelsize=12)  # fontsize of the tick labels
plt.rc('legend', fontsize=16)  # legend fontsize


def get_varname(varname, species=None):
    if species is None:
        return varname
    else:
        return varname + "_" + species


def main():
    path = "/Users/stephendiiorio/Documents/Research/1D_Sim/"
    fnums = ["00056", "00060", "00070", "00088"]
    fname = []
    for n in fnums:
        fname.append(path + n + ".sdf")

    fig, axarr = plt.subplots(2, 2, sharex='col', sharey='row')
    fig.set_facecolor("w")
    txt = fig.suptitle('Number Density Evolution', verticalalignment='bottom')

    # only have a single file to plot, so we so this little hack since
    # axarr does not come out in an array in this case
    if not isinstance(axarr, np.ndarray):
        axarr = np.array([axarr])

    e_data = []
    ar_data = []
    ar8_data = []
    t_data = []

    for i in range(len(fname)):
        sdfdata = sdf.read(fname[i])

        e_data.append(sdfdata.Derived_Number_Density_Electron)
        ar_data.append(sdfdata.Derived_Number_Density_Argon)
        ar8_data.append(sdfdata.Derived_Number_Density_Argon8)
        t_data.append(sdfdata.Header['time'])
        x_axis = e_data[i].grid_mid.data[0]

    tmult0, tsym0 = get_si_prefix(t_data[0])
    tmult1, tsym1 = get_si_prefix(t_data[1])
    tmult2, tsym2 = get_si_prefix(t_data[2])
    tmult3, tsym3 = get_si_prefix(t_data[3])

    l1, = axarr[0, 0].plot(x_axis, e_data[0].data, 'k', label='Electron')
    l2, = axarr[0, 0].plot(x_axis, ar_data[0].data, 'r', label='Argon')
    l3, = axarr[0, 0].plot(x_axis, ar8_data[0].data, 'b', label='Argon$^{+8}$')
    axarr[0, 0].set_title(('t={0:.2f} $' + tsym0 + 's$').format(t_data[0] * tmult0))
    axarr[0, 1].plot(x_axis, e_data[1].data, 'k')
    axarr[0, 1].plot(x_axis, ar_data[1].data, 'r')
    axarr[0, 1].plot(x_axis, ar8_data[0].data, 'b')
    axarr[0, 1].set_title(('t={0:.2f} $' + tsym1 + 's$').format(t_data[1] * tmult1))
    axarr[1, 0].plot(x_axis, e_data[2].data, 'k')
    axarr[1, 0].plot(x_axis, ar_data[2].data, 'r')
    axarr[1, 0].plot(x_axis, ar8_data[0].data, 'b')
    axarr[1, 0].set_title(('t={0:.2f} $' + tsym2 + 's$').format(t_data[2] * tmult2))
    axarr[1, 1].plot(x_axis, e_data[3].data, 'k')
    axarr[1, 1].plot(x_axis, ar_data[3].data, 'r')
    axarr[1, 1].plot(x_axis, ar8_data[0].data, 'b')
    axarr[1, 1].set_title(('t={0:.2f} $' + tsym3 + 's$').format(t_data[3] * tmult3))

    xmin, xmax = get_var_range(x_axis)
    xmult, xsym = get_si_prefix(xmax - xmin)
    axarr[1, 0].xaxis.set_major_formatter(FuncFormatter(lambda x, y: '{0:g}'.format(x * xmult)))
    axarr[1, 1].xaxis.set_major_formatter(FuncFormatter(lambda x, y: '{0:g}'.format(x * xmult)))

    ymult1, ysym1 = get_si_prefix(e_data[0].data.max() - e_data[0].data.min())
    ymult2, ysym2 = get_si_prefix(e_data[2].data.max() - e_data[2].data.min())
    axarr[0, 0].yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult1)))
    axarr[1, 0].yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult2)))

    axarr[0, 0].set(ylabel='Number Density' + ' $(' + ysym1 + 'm^{-3})$')
    axarr[1, 0].set(xlabel='x' + ' $(' + xsym + 'm)$', ylabel='Number Density' + ' $(' + ysym2 + 'm^{-3})$')
    axarr[1, 1].set(xlabel='x' + ' $(' + xsym + 'm)$')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axarr.flat:
        ax.label_outer()

    ls = [l1, l2, l3]
    labels = [l.get_label() for l in ls]
    lgd = axarr[0, 1].legend(ls, labels, bbox_to_anchor=(1.0, 1.045), loc=2)#loc='upper right')

    plt.savefig('numdens.png', dpi=600, bbox_extra_artists=(lgd,txt,), bbox_inches = "tight")


if __name__ == "__main__":
    main()
