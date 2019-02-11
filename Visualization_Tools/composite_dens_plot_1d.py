import glob
import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import use
from matplotlib.ticker import FuncFormatter

import sdf
from PlottingTools import get_si_prefix, get_var_range_from_sdf_files, calculate_aspect

use('Agg')


plt.rc('font', size=20)        # controls default text sizes
plt.rc('axes', titlesize=18)   # fontsize of the axes title
plt.rc('axes', labelsize=15)   # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)  # fontsize of the tick labels
plt.rc('ytick', labelsize=12)  # fontsize of the tick labels
plt.rc('legend', fontsize=16)  # legend fontsize


def get_files(wkdir='Data'):
    """
    Get a list of SDF filenames belonging to the same run
    """
    flist = glob.glob(wkdir + "/[0-9]*.sdf")
    flist = sorted(flist)

    return flist


def clean_file_list(file_list, varname):
    new_file_list = []

    for f in file_list:
        try:
            data = sdf.read(f)
            data.__dict__[varname]
            new_file_list.append(f)
        except KeyError:
            pass
    return new_file_list


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
        d = sdf.read(f)
        var = d.__dict__[varname]
        data.append(var.data)
    data = np.asarray(data)
    data = data.T

    tmin = sdf.read(file_list[0]).Header['time']
    tmax = sdf.read(file_list[-1]).Header['time']
    grid = var.grid_mid
    xmin = np.min(grid.data[0])
    xmax = np.max(grid.data[0])

    shape = data.shape
    extent = [tmin, tmax, xmax, xmin]

    xmult, xsym = get_si_prefix(xmax - xmin)  # y axis
    tmult, tsym = get_si_prefix(tmax - tmin)  # x axis

    if vmin is None and vmax is None:
        vmin, vmax = get_var_range_from_sdf_files(file_list, varname)
    elif vmin is None:
        vmin = get_var_range_from_sdf_files(file_list, varname)[0]
    elif vmax is None:
        vmax = get_var_range_from_sdf_files(file_list, varname)[1]
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
