from matplotlib import use
use('Agg')

import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import FuncFormatter

import sdf
from PlottingTools import get_si_prefix, get_var_range

try:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except:
    try:
        # Workaround for broken macOS installation
        import os
        import sys
        import matplotlib
        sys.path.append(os.path.join(matplotlib.__path__[0],
                                     '..', 'mpl_toolkits'))
        from axes_grid1 import make_axes_locatable
    except:
        pass


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
    path = "/scratch/lsa_flux/diiorios/Exp_Ar_32618/Full_smaller/"
    fnums = ["00016", "00060"]
    fname = []
    for n in fnums:
        fname.append(path + n + ".sdf")

    fig, axarr = plt.subplots(1, 2, sharey='row')
    plt.set_cmap(cm.coolwarm)
    fig.set_facecolor("w")
    # txt = fig.suptitle('Laser Focus and Initial Electric Field', verticalalignment='bottom')

    # only have a single file to plot, so we so this little hack since
    # axarr does not come out in an array in this case
    if not isinstance(axarr, np.ndarray):
        axarr = np.array([axarr])

    x_data = []
    y_data = []
    e_data = []
    t_data = []
    i0 = 1
    i1 = 0

    for i in range(len(fname)):
        sdfdata = sdf.read(fname[i])

        e_data.append(sdfdata.Electric_Field_Ey)
        x_data.append(e_data[i].grid.data[i0])
        y_data.append(e_data[i].grid.data[i1])
        t_data.append(sdfdata.Header['time'])

    tmult0, tsym0 = get_si_prefix(t_data[0])
    tmult1, tsym1 = get_si_prefix(t_data[1])

    axarr[0].pcolormesh(x_data[0], y_data[0], e_data[0].data)
    axarr[0].set_title(('t={0:.2f} $' + tsym0 + 's$').format(t_data[0] * tmult0))
    im = axarr[1].pcolormesh(x_data[1], y_data[1], e_data[1].data)
    axarr[1].set_title(('t={0:.2f} $' + tsym1 + 's$').format(t_data[1] * tmult1))

    xmin1, xmax1 = get_var_range(x_data[0])
    xmult1, xsym1 = get_si_prefix(xmax1 - xmin1)
    xmin2, xmax2 = get_var_range(x_data[1])
    xmult2, xsym2 = get_si_prefix(xmax2 - xmin2)

    ymin1, ymax1 = get_var_range(y_data[0])
    ymult1, ysym1 = get_si_prefix(ymax1 - ymin1)
    ymin2, ymax2 = get_var_range(y_data[1])
    ymult2, ysym2 = get_si_prefix(ymax2 - ymin2)

    axarr[0].xaxis.set_major_formatter(FuncFormatter(lambda x, y: '{0:g}'.format(x * xmult1)))
    axarr[1].xaxis.set_major_formatter(FuncFormatter(lambda x, y: '{0:g}'.format(x * xmult2)))

    axarr[0].yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult1)))

    axarr[0].set(xlabel='x' + ' $(' + xsym1 + 'm)$', ylabel='y' + ' $(' + ysym1 + 'm)$')
    axarr[1].set(xlabel='x' + ' $(' + xsym2 + 'm)$')

    axarr[0].set_ylim([ymin1 - (7 / ymult1), ymax1 - (7 / ymult1)])
    axarr[1].set_ylim([ymin2 - (7 / ymult1), ymax2 - (7 / ymult1)])

    plt.set_cmap(cm.coolwarm)

    divider = make_axes_locatable(axarr[1])
    cax = divider.append_axes("right", "5%", pad="15%")
    cbar = fig.colorbar(im, cax=cax, ax=axarr[1],
                        label='Electric Field Direction',
                        ticks=[e_data[1].data.min(), e_data[1].data.max()])
    cbar.ax.set_yticklabels(['$(-)$', '$(+)$'])

    plt.savefig('2defield.png', dpi=600, # bbox_extra_artists=(txt,),
        bbox_inches="tight")


if __name__ == "__main__":
    main()
