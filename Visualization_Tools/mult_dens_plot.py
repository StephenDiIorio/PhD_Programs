import sdf_helper as sh
import numpy as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
plt.ion()


plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=15)     # fontsize of the x and y labels
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

    for i in xrange(len(fname)):
        sdfdata = sh.getdata(fname[i])
        # sh.list_variables(sdfdata)

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
