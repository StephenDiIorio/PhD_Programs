from matplotlib import use
use('Agg')

import glob
import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

from CoordinateTransforms import cart2polar, reproject_image_into_polar
from PlottingTools import get_si_prefix, calculate_aspect
import sdf


plt.rc('font', size=20)        # controls default text sizes
plt.rc('axes', titlesize=18)   # fontsize of the axes title
plt.rc('axes', labelsize=15)   # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)  # fontsize of the tick labels
plt.rc('ytick', labelsize=12)  # fontsize of the tick labels
plt.rc('legend', fontsize=16)  # legend fontsize


def get_files(wkdir='Data'):
    """Get a list of SDF filenames belonging to the same run"""
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


def electric_radial_average(sdf_data):
    ex = sdf_data.__dict__['Electric_Field_Ex']
    ey = sdf_data.__dict__['Electric_Field_Ey']

    x_grid = ex.grid_mid.data[0]
    y_grid = ex.grid_mid.data[1]
    ex = ex.data
    ey = ey.data

    X, Y = np.meshgrid(x_grid, y_grid)
    R, T = cart2polar(X, Y)

    er = np.multiply(ex, np.cos(T)) + np.multiply(ey, np.sin(T))
    et = np.multiply(ey, np.cos(T)) - np.multiply(ex, np.sin(T))

    o, r, t = reproject_image_into_polar(er, origin=None, Jacobian=False, dr=1, dt=None)

    o_ma = np.ma.masked_equal(o, 0.)
    avg = np.average(o_ma, axis=1)

    return avg, R, T


def composite_field_plot(varname, vmin=None, vmax=None, directory='Data'):
    global verbose, dpi

    file_list = get_files(wkdir=directory)
    file_list = clean_file_list(file_list, varname)

    file_list.remove(directory + '0000.sdf')
    file_list.remove(directory + '0002.sdf')
    file_list.remove(directory + '0003.sdf')
    file_list.remove(directory + '0004.sdf')
    file_list.remove(directory + '0005.sdf')
    file_list.remove(directory + '0006.sdf')
    file_list.remove(directory + '0007.sdf')
    file_list.remove(directory + '0008.sdf')
    file_list.remove(directory + '0009.sdf')
    file_list.remove(directory + '0011.sdf')
    file_list.remove(directory + '0012.sdf')
    file_list.remove(directory + '0013.sdf')
    file_list.remove(directory + '0014.sdf')
    file_list.remove(directory + '0015.sdf')
    file_list.remove(directory + '0016.sdf')
    file_list.remove(directory + '0017.sdf')
    file_list.remove(directory + '0018.sdf')
    file_list.remove(directory + '0019.sdf')
    file_list.remove(directory + '0021.sdf')
    file_list.remove(directory + '0022.sdf')
    file_list.remove(directory + '0023.sdf')
    file_list.remove(directory + '0024.sdf')
    file_list.remove(directory + '0025.sdf')
    file_list.remove(directory + '0026.sdf')
    file_list.remove(directory + '0027.sdf')
    file_list.remove(directory + '0028.sdf')
    file_list.remove(directory + '0029.sdf')
    file_list.remove(directory + '0031.sdf')
    file_list.remove(directory + '0032.sdf')
    file_list.remove(directory + '0033.sdf')
    file_list.remove(directory + '0034.sdf')
    file_list.remove(directory + '0035.sdf')
    file_list.remove(directory + '0036.sdf')
    file_list.remove(directory + '0037.sdf')
    file_list.remove(directory + '0038.sdf')
    file_list.remove(directory + '0039.sdf')
    file_list.remove(directory + '0041.sdf')
    file_list.remove(directory + '0042.sdf')
    file_list.remove(directory + '0043.sdf')
    file_list.remove(directory + '0044.sdf')
    file_list.remove(directory + '0045.sdf')
    file_list.remove(directory + '0046.sdf')
    file_list.remove(directory + '0047.sdf')
    file_list.remove(directory + '0048.sdf')
    file_list.remove(directory + '0049.sdf')
    file_list.remove(directory + '0051.sdf')
    file_list.remove(directory + '0052.sdf')
    file_list.remove(directory + '0053.sdf')
    file_list.remove(directory + '0054.sdf')
    file_list.remove(directory + '0055.sdf')
    file_list.remove(directory + '0056.sdf')
    file_list.remove(directory + '0057.sdf')
    file_list.remove(directory + '0058.sdf')
    file_list.remove(directory + '0059.sdf')
    file_list.remove(directory + '0061.sdf')
    file_list.remove(directory + '0062.sdf')
    file_list.remove(directory + '0063.sdf')
    file_list.remove(directory + '0064.sdf')
    file_list.remove(directory + '0065.sdf')
    file_list.remove(directory + '0066.sdf')
    file_list.remove(directory + '0067.sdf')
    file_list.remove(directory + '0068.sdf')
    file_list.remove(directory + '0069.sdf')
    file_list.remove(directory + '0071.sdf')
    file_list.remove(directory + '0072.sdf')
    file_list.remove(directory + '0073.sdf')
    file_list.remove(directory + '0074.sdf')
    file_list.remove(directory + '0075.sdf')
    file_list.remove(directory + '0076.sdf')
    file_list.remove(directory + '0077.sdf')
    file_list.remove(directory + '0078.sdf')
    file_list.remove(directory + '0079.sdf')
    file_list.remove(directory + '0081.sdf')
    file_list.remove(directory + '0082.sdf')
    file_list.remove(directory + '0083.sdf')
    file_list.remove(directory + '0084.sdf')
    file_list.remove(directory + '0085.sdf')
    file_list.remove(directory + '0086.sdf')
    file_list.remove(directory + '0087.sdf')
    file_list.remove(directory + '0088.sdf')
    file_list.remove(directory + '0089.sdf')
    file_list.remove(directory + '0091.sdf')
    file_list.remove(directory + '0092.sdf')
    file_list.remove(directory + '0093.sdf')
    file_list.remove(directory + '0094.sdf')
    file_list.remove(directory + '0095.sdf')
    file_list.remove(directory + '0096.sdf')
    file_list.remove(directory + '0097.sdf')
    file_list.remove(directory + '0098.sdf')
    file_list.remove(directory + '0099.sdf')

    if verbose > 0:
        print('Found {} files to plot'.format(len(file_list)))

    data = []
    for f in file_list:
        d = sdf.read(f)
        avg, r, t = electric_radial_average(d)
        data.append(avg)
    var = d.__dict__[varname]  # for units later
    data = np.asarray(data)
    data = data.T

    tmin = sdf.read(file_list[0]).Header['time']
    tmax = sdf.read(file_list[-1]).Header['time']
    rmin = np.min(r)
    rmax = np.max(r)
    print(rmin, rmax)

    shape = data.shape
    extent = [tmin, tmax, rmin, rmax]

    rmult, rsym = get_si_prefix(rmax - rmin)  # y axis
    tmult, tsym = get_si_prefix(tmax - tmin)  # x axis

    if vmin is None and vmax is None:
        vmin = np.min(data)
        vmax = np.max(data)
    elif vmin is None:
        vmin = np.min(data)
    elif vmax is None:
        vmax = np.max(data)
    mult, sym = get_si_prefix(vmax - vmin)

    fig, ax = plt.subplots()
    im = ax.imshow(data, extent=extent, aspect=calculate_aspect(shape, extent), interpolation='none', cmap=cm.coolwarm, vmin=vmin, vmax=vmax, origin='lower')

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * tmult)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * rmult)))
    plt.xlabel('t $(' + tsym + 's)$')
    plt.ylabel('r $(' + rsym + 'm)$')
    # data_label = var.name + ' $(' + sym + var.units + ')$'
    data_label = 'Radial Electric Field $(' + sym + var.units + ')$'
    plt.title('Radial Electric Field Evolution')

    cbar = fig.colorbar(im, label=data_label, format=FuncFormatter(lambda x, y: x * mult))
    plt.tight_layout()
    plt.savefig('rad_test.png', dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    import argparse
    global verbose, dpi

    verbose = 0
    dpi = 300

    varname = 'Electric_Field_Ex'
    vmin = None
    vmax = None

    parser = argparse.ArgumentParser(description='''
    This composite plot of variable evolution
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

    composite_field_plot(var, vmin=vmin, vmax=vmax, directory=args.directory[0])
