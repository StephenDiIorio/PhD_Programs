import sdf_helper as sh
import numpy as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.cm as cm
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
plt.ion()


plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=15)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize


def get_si_prefix(scale, full_units=False):
    scale = abs(scale)
    mult = 1
    sym = ''

    if scale < 1e-24:
        full_units = True
    elif scale < 1e-21:
        # yocto
        mult = 1e24
        sym = 'y'
    elif scale < 1e-19:
        # zepto
        mult = 1e21
        sym = 'z'
    elif scale < 1e-16:
        # atto
        mult = 1e18
        sym = 'a'
    elif scale < 1e-13:
        # femto
        mult = 1e15
        sym = 'f'
    elif scale < 1e-10:
        # pico
        mult = 1e12
        sym = 'p'
    elif scale < 1e-7:
        # nano
        mult = 1e9
        sym = 'n'
    elif scale < 1e-4:
        # micro
        mult = 1e6
        sym = '{\mu}'
    elif scale < 1e-1:
        # milli
        mult = 1e3
        sym = 'm'
    elif scale >= 1e27:
        full_units = True
    elif scale >= 1e24:
        # yotta
        mult = 1e-24
        sym = 'Y'
    elif scale >= 1e21:
        # zetta
        mult = 1e-21
        sym = 'Z'
    elif scale >= 1e18:
        # exa
        mult = 1e-18
        sym = 'E'
    elif scale >= 1e15:
        # peta
        mult = 1e-15
        sym = 'P'
    elif scale >= 1e12:
        # tera
        mult = 1e-12
        sym = 'T'
    elif scale >= 1e9:
        # giga
        mult = 1e-9
        sym = 'G'
    elif scale >= 1e6:
        # mega
        mult = 1e-6
        sym = 'M'
    elif scale >= 1e3:
        # kilo
        mult = 1e-3
        sym = 'k'

    if full_units:
        scale = scale * mult
        if scale <= 0:
            pwr = 0
        else:
            pwr = (-np.floor(np.log10(scale)))

        mult = mult * np.power(10.0, pwr)
        if np.rint(pwr) != 0:
            sym = "(10^{%.0f})" % (-pwr) + sym

    return mult, sym


def get_var_range(sdfdata):
    """Get a the data range for a given variable across an entire run"""

    vmin = float("inf")
    vmax = -float("inf")

    var_min = sdfdata.min()
    var_max = sdfdata.max()
    if var_min < vmin:
        vmin = var_min
    if var_max > vmax:
        vmax = var_max

    print vmin, vmax
    return vmin, vmax


def get_varname(varname, species=None):
    if species is None:
        return varname
    else:
        return varname + "_" + species


def main():
    path = "/Users/stephendiiorio/Desktop/"
    fnums = ["00016", "00060"]
    fname = []
    for n in fnums:
        fname.append(path + n + ".sdf")

    fig, axarr = plt.subplots(1, 2, sharey='row')
    plt.set_cmap(cm.coolwarm)
    fig.set_facecolor("w")
    txt = fig.suptitle('Laser Focus and Initial Electric Field', verticalalignment='bottom')

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

    for i in xrange(len(fname)):
        sdfdata = sh.getdata(fname[i])
        # sh.list_variables(sdfdata)

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

    axarr[0].set_ylim([ymin1, ymax1])
    axarr[1].set_ylim([ymin2, ymax2])

    plt.set_cmap(cm.coolwarm)

    divider = make_axes_locatable(axarr[1])
    cax = divider.append_axes("right", "5%", pad="15%")
    cbar = fig.colorbar(im, cax=cax, ax=axarr[1],
                        label='Electric Field Direction',
                        ticks=[e_data[1].data.min(), e_data[1].data.max()])
    cbar.ax.set_yticklabels(['$(-)$', '$(+)$'])

    plt.savefig('2defield.png', dpi=600, bbox_extra_artists=(txt,), bbox_inches="tight")


if __name__ == "__main__":
    main()
