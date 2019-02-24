from matplotlib import use
use('Agg')

import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

import sdf
from PlottingTools import shiftedColorMap
from Weighting import first_order_weight_2d


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


def in_domain(p, x, y):
    return p[0] >= np.amin(x) and p[0] <= np.amax(x) and p[1] >= np.amin(y) and p[1] <= np.amax(y)


def main():
    species = "Electron"

    path = "/scratch/lsa_flux/diiorios/2d_run/"
    fnums = ["0200"]
    fname = []
    for n in fnums:
        fname.append(path + n + ".sdf")

    # x_axis_num = 0
    y_axis_num = 1
    # z_axis_num = 2

    fig, axarr = plt.subplots(len(fname), sharex=True, sharey=True)
    # only have a single file to plot, so we so this little hack since
    # axarr does not come out in an array in this case
    if not isinstance(axarr, np.ndarray):
        axarr = np.array([axarr])

    axarr[0].set_title(species + " files " + str(fnums))

    for i in range(len(fname)):
        sdfdata = sdf.read(fname[i])
        print(sdfdata.Header['time'])

        grid = sdfdata.__dict__["Grid_Grid_mid"].data

        x = grid[0]
        y = grid[1]

        dx = grid[0][1] - grid[0][0]
        dy = grid[1][1] - grid[1][0]

        p_pos = sdfdata.__dict__[get_varname("Grid_Particles", species)].data
        p_list = list(zip(p_pos[0], p_pos[1]))

        p_list[:] = [p for p in p_list if in_domain(p, x, y)]

        w = sdfdata.__dict__[get_varname("Particles_Weight", species)].data

        num_dens = first_order_weight_2d(x, y, dx, dy, p_list, weight=w)
        real_num_dens = sdfdata.__dict__[get_varname("Derived_Number_Density", species)].data

        d = np.subtract(num_dens, real_num_dens)

    plt.figure()
    plt.pcolormesh(d, cmap=cm.coolwarm)
    cbar = plt.colorbar()
    plt.savefig('den.png')

    return


if __name__ == "__main__":
    main()
