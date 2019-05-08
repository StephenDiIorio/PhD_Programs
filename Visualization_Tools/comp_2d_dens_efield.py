from matplotlib import use
use('Agg')

import glob
import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

from CoordinateTransforms import cart2polar, reproject_image_into_polar
from PlottingTools import get_si_prefix, get_var_range
import sdf


# plt.rc('font', size=20)        # controls default text sizes
# plt.rc('axes', titlesize=18)   # fontsize of the axes title
# plt.rc('axes', labelsize=15)   # fontsize of the x and y labels
# plt.rc('xtick', labelsize=12)  # fontsize of the tick labels
# plt.rc('ytick', labelsize=12)  # fontsize of the tick labels
# plt.rc('legend', fontsize=16)  # legend fontsize


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


def density_radial_average(sdf_data, varname):
    dens = sdf_data.__dict__[varname]

    x_grid = dens.grid_mid.data[0]
    y_grid = dens.grid_mid.data[1]
    dens = dens.data

    X, Y = np.meshgrid(x_grid, y_grid)
    R, T = cart2polar(X, Y)

    o, r, t = reproject_image_into_polar(dens, origin=None, Jacobian=False, dr=1, dt=None)

    o_ma = np.ma.masked_equal(o, 0.)
    avg = np.average(o_ma, axis=1)

    return avg, R, T


def main():
    directory = '/scratch/lsa_flux/diiorios/2d_run/'
    field = 'Electric_Field_Ex'
    spec = 'Electron'
    den = 'Derived_Number_Density_' + spec

    file_list = get_files(wkdir=directory)
    e_file_list = clean_file_list(file_list, field)
    d_file_list = clean_file_list(file_list, den)

    e_file_list.remove(directory + '0000.sdf')
    e_file_list.remove(directory + '0002.sdf')
    e_file_list.remove(directory + '0003.sdf')
    e_file_list.remove(directory + '0004.sdf')
    e_file_list.remove(directory + '0005.sdf')
    e_file_list.remove(directory + '0006.sdf')
    e_file_list.remove(directory + '0007.sdf')
    e_file_list.remove(directory + '0008.sdf')
    e_file_list.remove(directory + '0009.sdf')
    e_file_list.remove(directory + '0011.sdf')
    e_file_list.remove(directory + '0012.sdf')
    e_file_list.remove(directory + '0013.sdf')
    e_file_list.remove(directory + '0014.sdf')
    e_file_list.remove(directory + '0015.sdf')
    e_file_list.remove(directory + '0016.sdf')
    e_file_list.remove(directory + '0017.sdf')
    e_file_list.remove(directory + '0018.sdf')
    e_file_list.remove(directory + '0019.sdf')
    e_file_list.remove(directory + '0021.sdf')
    e_file_list.remove(directory + '0022.sdf')
    e_file_list.remove(directory + '0023.sdf')
    e_file_list.remove(directory + '0024.sdf')
    e_file_list.remove(directory + '0025.sdf')
    e_file_list.remove(directory + '0026.sdf')
    e_file_list.remove(directory + '0027.sdf')
    e_file_list.remove(directory + '0028.sdf')
    e_file_list.remove(directory + '0029.sdf')
    e_file_list.remove(directory + '0031.sdf')
    e_file_list.remove(directory + '0032.sdf')
    e_file_list.remove(directory + '0033.sdf')
    e_file_list.remove(directory + '0034.sdf')
    e_file_list.remove(directory + '0035.sdf')
    e_file_list.remove(directory + '0036.sdf')
    e_file_list.remove(directory + '0037.sdf')
    e_file_list.remove(directory + '0038.sdf')
    e_file_list.remove(directory + '0039.sdf')
    e_file_list.remove(directory + '0041.sdf')
    e_file_list.remove(directory + '0042.sdf')
    e_file_list.remove(directory + '0043.sdf')
    e_file_list.remove(directory + '0044.sdf')
    e_file_list.remove(directory + '0045.sdf')
    e_file_list.remove(directory + '0046.sdf')
    e_file_list.remove(directory + '0047.sdf')
    e_file_list.remove(directory + '0048.sdf')
    e_file_list.remove(directory + '0049.sdf')
    e_file_list.remove(directory + '0051.sdf')
    e_file_list.remove(directory + '0052.sdf')
    e_file_list.remove(directory + '0053.sdf')
    e_file_list.remove(directory + '0054.sdf')
    e_file_list.remove(directory + '0055.sdf')
    e_file_list.remove(directory + '0056.sdf')
    e_file_list.remove(directory + '0057.sdf')
    e_file_list.remove(directory + '0058.sdf')
    e_file_list.remove(directory + '0059.sdf')
    e_file_list.remove(directory + '0061.sdf')
    e_file_list.remove(directory + '0062.sdf')
    e_file_list.remove(directory + '0063.sdf')
    e_file_list.remove(directory + '0064.sdf')
    e_file_list.remove(directory + '0065.sdf')
    e_file_list.remove(directory + '0066.sdf')
    e_file_list.remove(directory + '0067.sdf')
    e_file_list.remove(directory + '0068.sdf')
    e_file_list.remove(directory + '0069.sdf')
    e_file_list.remove(directory + '0071.sdf')
    e_file_list.remove(directory + '0072.sdf')
    e_file_list.remove(directory + '0073.sdf')
    e_file_list.remove(directory + '0074.sdf')
    e_file_list.remove(directory + '0075.sdf')
    e_file_list.remove(directory + '0076.sdf')
    e_file_list.remove(directory + '0077.sdf')
    e_file_list.remove(directory + '0078.sdf')
    e_file_list.remove(directory + '0079.sdf')
    e_file_list.remove(directory + '0081.sdf')
    e_file_list.remove(directory + '0082.sdf')
    e_file_list.remove(directory + '0083.sdf')
    e_file_list.remove(directory + '0084.sdf')
    e_file_list.remove(directory + '0085.sdf')
    e_file_list.remove(directory + '0086.sdf')
    e_file_list.remove(directory + '0087.sdf')
    e_file_list.remove(directory + '0088.sdf')
    e_file_list.remove(directory + '0089.sdf')
    e_file_list.remove(directory + '0091.sdf')
    e_file_list.remove(directory + '0092.sdf')
    e_file_list.remove(directory + '0093.sdf')
    e_file_list.remove(directory + '0094.sdf')
    e_file_list.remove(directory + '0095.sdf')
    e_file_list.remove(directory + '0096.sdf')
    e_file_list.remove(directory + '0097.sdf')
    e_file_list.remove(directory + '0098.sdf')
    e_file_list.remove(directory + '0099.sdf')

    d_file_list.remove(directory + '0000.sdf')
    d_file_list.remove(directory + '0002.sdf')
    d_file_list.remove(directory + '0003.sdf')
    d_file_list.remove(directory + '0004.sdf')
    d_file_list.remove(directory + '0005.sdf')
    d_file_list.remove(directory + '0006.sdf')
    d_file_list.remove(directory + '0007.sdf')
    d_file_list.remove(directory + '0008.sdf')
    d_file_list.remove(directory + '0009.sdf')
    d_file_list.remove(directory + '0011.sdf')
    d_file_list.remove(directory + '0012.sdf')
    d_file_list.remove(directory + '0013.sdf')
    d_file_list.remove(directory + '0014.sdf')
    d_file_list.remove(directory + '0015.sdf')
    d_file_list.remove(directory + '0016.sdf')
    d_file_list.remove(directory + '0017.sdf')
    d_file_list.remove(directory + '0018.sdf')
    d_file_list.remove(directory + '0019.sdf')
    d_file_list.remove(directory + '0021.sdf')
    d_file_list.remove(directory + '0022.sdf')
    d_file_list.remove(directory + '0023.sdf')
    d_file_list.remove(directory + '0024.sdf')
    d_file_list.remove(directory + '0025.sdf')
    d_file_list.remove(directory + '0026.sdf')
    d_file_list.remove(directory + '0027.sdf')
    d_file_list.remove(directory + '0028.sdf')
    d_file_list.remove(directory + '0029.sdf')
    d_file_list.remove(directory + '0031.sdf')
    d_file_list.remove(directory + '0032.sdf')
    d_file_list.remove(directory + '0033.sdf')
    d_file_list.remove(directory + '0034.sdf')
    d_file_list.remove(directory + '0035.sdf')
    d_file_list.remove(directory + '0036.sdf')
    d_file_list.remove(directory + '0037.sdf')
    d_file_list.remove(directory + '0038.sdf')
    d_file_list.remove(directory + '0039.sdf')
    d_file_list.remove(directory + '0041.sdf')
    d_file_list.remove(directory + '0042.sdf')
    d_file_list.remove(directory + '0043.sdf')
    d_file_list.remove(directory + '0044.sdf')
    d_file_list.remove(directory + '0045.sdf')
    d_file_list.remove(directory + '0046.sdf')
    d_file_list.remove(directory + '0047.sdf')
    d_file_list.remove(directory + '0048.sdf')
    d_file_list.remove(directory + '0049.sdf')
    d_file_list.remove(directory + '0051.sdf')
    d_file_list.remove(directory + '0052.sdf')
    d_file_list.remove(directory + '0053.sdf')
    d_file_list.remove(directory + '0054.sdf')
    d_file_list.remove(directory + '0055.sdf')
    d_file_list.remove(directory + '0056.sdf')
    d_file_list.remove(directory + '0057.sdf')
    d_file_list.remove(directory + '0058.sdf')
    d_file_list.remove(directory + '0059.sdf')
    d_file_list.remove(directory + '0061.sdf')
    d_file_list.remove(directory + '0062.sdf')
    d_file_list.remove(directory + '0063.sdf')
    d_file_list.remove(directory + '0064.sdf')
    d_file_list.remove(directory + '0065.sdf')
    d_file_list.remove(directory + '0066.sdf')
    d_file_list.remove(directory + '0067.sdf')
    d_file_list.remove(directory + '0068.sdf')
    d_file_list.remove(directory + '0069.sdf')
    d_file_list.remove(directory + '0071.sdf')
    d_file_list.remove(directory + '0072.sdf')
    d_file_list.remove(directory + '0073.sdf')
    d_file_list.remove(directory + '0074.sdf')
    d_file_list.remove(directory + '0075.sdf')
    d_file_list.remove(directory + '0076.sdf')
    d_file_list.remove(directory + '0077.sdf')
    d_file_list.remove(directory + '0078.sdf')
    d_file_list.remove(directory + '0079.sdf')
    d_file_list.remove(directory + '0081.sdf')
    d_file_list.remove(directory + '0082.sdf')
    d_file_list.remove(directory + '0083.sdf')
    d_file_list.remove(directory + '0084.sdf')
    d_file_list.remove(directory + '0085.sdf')
    d_file_list.remove(directory + '0086.sdf')
    d_file_list.remove(directory + '0087.sdf')
    d_file_list.remove(directory + '0088.sdf')
    d_file_list.remove(directory + '0089.sdf')
    d_file_list.remove(directory + '0091.sdf')
    d_file_list.remove(directory + '0092.sdf')
    d_file_list.remove(directory + '0093.sdf')
    d_file_list.remove(directory + '0094.sdf')
    d_file_list.remove(directory + '0095.sdf')
    d_file_list.remove(directory + '0096.sdf')
    d_file_list.remove(directory + '0097.sdf')
    d_file_list.remove(directory + '0098.sdf')
    d_file_list.remove(directory + '0099.sdf')


    print('Found {} files to plot'.format(len(file_list)))

    e_data = []
    d_data = []
    for f in e_file_list:
        d = sdf.read(f)
        avg, r, dummy_t = electric_radial_average(d)
        e_data.append(avg)
    for f in d_file_list:
        d = sdf.read(f)
        avg, r, dummy_t = density_radial_average(d, den)
        d_data.append(avg)
    e_var = d.__dict__[field]  # for units later
    d_var = d.__dict__[den]  # for units later

    e_data = np.asarray(e_data)
    e_data = e_data.T
    d_data = np.asarray(d_data)
    d_data = d_data.T

    tmin = sdf.read(e_file_list[0]).Header['time']
    tmax = sdf.read(e_file_list[-1]).Header['time']
    t_data = np.linspace(tmin, tmax, np.shape(e_data)[1])
    rmin = np.min(r)
    rmax = np.max(r)
    r_data = np.linspace(rmin, rmax, np.shape(e_data)[0])
    print(rmin, rmax)

    # shape = e_data.shape
    # extent = [tmin, tmax, rmin, rmax]

    rmult, rsym = get_si_prefix(rmax - rmin)  # y axis
    tmult, tsym = get_si_prefix(tmax - tmin)  # x axis

    emin, emax = get_var_range(e_data)
    emax = 3e9
    emult, esym = get_si_prefix(emax - emin)
    dmin, dmax = get_var_range(d_data)
    # dmin = max(d_data.min(), 0.1)
    dmax = 0.5e23
    dmult, dsym = get_si_prefix(dmax - dmin)

    fig, axarr = plt.subplots(2, 1, sharex='col')
    # plt.set_cmap(cm.coolwarm)
    fig.set_facecolor("w")

    e_im = axarr[0].pcolormesh(t_data, r_data, e_data,
                            #    norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=emin, vmax=emax),
                               cmap=cm.coolwarm,
                               vmin=emin, vmax=emax)
    axarr[0].set_title('(a)', loc='left')
    d_im = axarr[1].pcolormesh(t_data, r_data, d_data,
                            #    norm=colors.LogNorm(vmin=dmin, vmax=dmax),
                               cmap=cm.plasma,
                               vmin=dmin, vmax=dmax)
    axarr[1].set_title('(b)', loc='left')

    e_label = '$E_{r}$ $(' + esym + e_var.units + ')$'
    d_label = 'Radial $n_{e}$ $(' + dsym + d_var.units + ')$'
    fig.colorbar(e_im, ax=axarr[0], label=e_label, format=FuncFormatter(lambda x, y: x * emult))
    fig.colorbar(d_im, ax=axarr[1], label=d_label, format=FuncFormatter(lambda x, y: x * dmult))


    axarr[1].xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * tmult)))
    axarr[0].yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * rmult)))
    axarr[1].yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * rmult)))

    axarr[0].set(ylabel='r $(' + rsym + 'm)$')
    axarr[1].set(xlabel='t $(' + tsym + 's)$', ylabel='r $(' + rsym + 'm)$')

    axarr[0].set_ylim([rmin, 40e-6])
    axarr[1].set_ylim([rmin, 40e-6])

    axarr[0].get_xaxis().set_tick_params(direction='in')
    axarr[1].get_xaxis().set_tick_params(direction='in')
    axarr[0].get_yaxis().set_tick_params(direction='in')
    axarr[1].get_yaxis().set_tick_params(direction='in')


    # plt.tight_layout()
    plt.savefig('time_rad_evo.png', dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    main()
