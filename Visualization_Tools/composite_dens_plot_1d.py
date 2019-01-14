import sdf
import glob
import numpy as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FuncFormatter
# plt.ion()


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


def get_var_range(file_list, varname):
    """Get a the data range for a given variable across an entire run"""

    vmin = float("inf")
    vmax = -float("inf")

    for f in file_list:
        try:
            data = sdf.read(f, mmap=0)
            var = data.__dict__[varname].data
            var_min = var.min()
            var_max = var.max()
            if var_min < vmin:
                vmin = var_min
            if var_max > vmax:
                vmax = var_max
        except:
            pass

    if verbose > 0:
        print vmin, vmax

    return vmin, vmax


def get_files(wkdir='Data'):
    """Get a list of SDF filenames belonging to the same run"""
    flist = glob.glob(wkdir + "/[0-9]*.sdf")
    flist = sorted(flist)

    return flist


def clean_file_list(file_list, varname):
    new_file_list = []

    for f in file_list:
        try:
            data = sdf.read(f, mmap=0)
            data.__dict__[varname]
            new_file_list.append(f)
        except KeyError:
            pass
    return new_file_list


def calculate_aspect(shape, extent):
    dx = (extent[1] - extent[0]) / float(shape[1])
    dy = (extent[3] - extent[2]) / float(shape[0])
    if dx / dy > 0:
        return dx / dy
    else:
        return -dx / dy


def composite_field_plot(varname, vmin=None, vmax=None, directory='Data'):
    global verbose, dpi

    file_list = get_files(wkdir=directory)
    file_list = clean_file_list(file_list, varname)

    file_list.remove(directory + '00000.sdf')
    file_list.remove(directory + '00002.sdf')
    file_list.remove(directory + '00003.sdf')
    file_list.remove(directory + '00004.sdf')
    file_list.remove(directory + '00005.sdf')
    file_list.remove(directory + '00006.sdf')
    file_list.remove(directory + '00007.sdf')
    file_list.remove(directory + '00008.sdf')
    file_list.remove(directory + '00009.sdf')
    file_list.remove(directory + '00010.sdf')
    file_list.remove(directory + '00011.sdf')
    file_list.remove(directory + '00012.sdf')
    file_list.remove(directory + '00013.sdf')
    file_list.remove(directory + '00014.sdf')
    file_list.remove(directory + '00015.sdf')
    file_list.remove(directory + '00016.sdf')
    file_list.remove(directory + '00017.sdf')
    file_list.remove(directory + '00018.sdf')
    file_list.remove(directory + '00019.sdf')
    file_list.remove(directory + '00020.sdf')
    file_list.remove(directory + '00021.sdf')
    file_list.remove(directory + '00022.sdf')
    file_list.remove(directory + '00023.sdf')
    file_list.remove(directory + '00024.sdf')
    file_list.remove(directory + '00025.sdf')
    file_list.remove(directory + '00027.sdf')
    file_list.remove(directory + '00028.sdf')
    file_list.remove(directory + '00029.sdf')
    file_list.remove(directory + '00030.sdf')
    file_list.remove(directory + '00031.sdf')
    file_list.remove(directory + '00032.sdf')
    file_list.remove(directory + '00033.sdf')
    file_list.remove(directory + '00034.sdf')
    file_list.remove(directory + '00035.sdf')
    file_list.remove(directory + '00037.sdf')
    file_list.remove(directory + '00038.sdf')
    file_list.remove(directory + '00039.sdf')
    file_list.remove(directory + '00040.sdf')
    file_list.remove(directory + '00041.sdf')
    file_list.remove(directory + '00042.sdf')
    file_list.remove(directory + '00043.sdf')
    file_list.remove(directory + '00044.sdf')
    file_list.remove(directory + '00045.sdf')
    file_list.remove(directory + '00047.sdf')
    file_list.remove(directory + '00048.sdf')
    file_list.remove(directory + '00049.sdf')
    file_list.remove(directory + '00050.sdf')
    file_list.remove(directory + '00051.sdf')
    file_list.remove(directory + '00052.sdf')
    file_list.remove(directory + '00053.sdf')
    file_list.remove(directory + '00054.sdf')
    file_list.remove(directory + '00055.sdf')

    if verbose > 0:
        print('Found {} files to plot'.format(len(file_list)))

    data = []
    for f in file_list:
        d = sdf.read(f, mmap=0)
        var = d.__dict__[varname]
        data.append(var.data)
    data = np.asarray(data)
    data = data.T

    tmin = sdf.read(file_list[0], mmap=0).Header['time']
    tmax = sdf.read(file_list[-1], mmap=0).Header['time']
    grid = var.grid_mid
    xmin = np.min(grid.data[0])
    xmax = np.max(grid.data[0])

    shape = data.shape
    extent = [tmin, tmax, xmax, xmin]

    xmult, xsym = get_si_prefix(xmax - xmin)  # y axis
    tmult, tsym = get_si_prefix(tmax - tmin)  # x axis

    if vmin is None and vmax is None:
        vmin, vmax = get_var_range(file_list, varname)
    elif vmin is None:
        vmin = get_var_range(file_list, varname)[0]
    elif vmax is None:
        vmax = get_var_range(file_list, varname)[1]
    mult, sym = get_si_prefix(vmax - vmin)

    fig, ax = plt.subplots()
    im = ax.imshow(data, extent=extent, aspect=calculate_aspect(shape, extent), interpolation='none', cmap=cm.plasma, vmin=vmin, vmax=vmax)

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * tmult)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * xmult)))
    plt.xlabel('t $(' + tsym + 's)$')
    plt.ylabel(grid.labels[0] + ' $(' + xsym + grid.units[0] + ')$')
    # data_label = var.name + ' $(' + sym + var.units + ')$'
    data_label = 'Argon$^{+8}$ Number Density $(' + sym + var.units + ')$'
    plt.title('Electron Density Evolution')

    cbar = fig.colorbar(im, label=data_label, format=FuncFormatter(lambda x, y: x * mult))
    plt.tight_layout()
    plt.savefig('electron_comp_thermal_nocoll.png', dpi=600, bbox_inches="tight")
    # plt.show()


if __name__ == "__main__":
    import argparse
    global verbose, dpi

    verbose = 0
    dpi = 300

    varname = 'Electric_Field_Ex'
    vmin = None
    vmax = None

    parser = argparse.ArgumentParser(description='''
    This generates a movie from a series of SDF files
    ''')
    # parser.add_argument('variable', type=str, default=varname, nargs='?',
    #                     help="Variable to plot")
    # EDIT: TAKE IN MULTIPLE VARIABLES TO PLOT
    parser.add_argument('variable', type=str, default=varname, nargs='*',
                        help="Variable(s) to plot")
    parser.add_argument('-v', '--verbose', action="count", default=0,
                        help="Increase verbosity")
    # EDIT: OPTION TO SEND THE DIRECTORY NAME WITH DATA FILES
    parser.add_argument('-d', '--directory', type=str, default='Data', nargs=1,
                        help="Directory for one of the simulation dumps")
    parser.add_argument('-n', '--noneg', action="count", default=0,
                        help="Makes vmin positive")
    parser.add_argument('-i', '--vmin', type=float, default=vmin, nargs=1,
                        help="Vmin for data display (change min of colorbar)")
    parser.add_argument('-a', '--vmax', type=float, default=vmax, nargs=1,
                        help="Vmax for data display (change max of colorbar)")

    args = parser.parse_args()

    verbose = args.verbose

    if isinstance(args.variable, list):
        var = args.variable[0]
    else:
        var = args.variable
    if isinstance(args.vmin, list):
        if args.noneg:
            vmin = args.vmin[0]
        else:
            vmin = -args.vmin[0]
    else:
        vmin = args.vmin
    if isinstance(args.vmax, list):
        vmax = args.vmax[0]
    else:
        vmax = args.vmax

    # EDIT: CALL PLOT_FIGURES WITH DIRECTORY ARGUMENT
    # take first element of directory because parser gives args in a list
    composite_field_plot(var, vmin=vmin, vmax=vmax, directory=args.directory[0])
