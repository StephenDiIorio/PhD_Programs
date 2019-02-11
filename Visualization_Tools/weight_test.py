import glob
import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

import sdf
from PlottingTools import shiftedColorMap
from Weighting import first_order_weight_2d


def get_varname(varname, species=None):
    if species is None:
        return varname
    else:
        return varname + "_" + species


def in_domain(p, x, y):
    if p[0] < np.amin(x) or p[0] > np.amax(x):
        return False
    elif p[1] < np.amin(y) or p[1] > np.amax(y):
        return False
    else:
        return True


def main():
    species = "Argon1"

    path = "/Users/stephendiiorio/Desktop/"
    fnums = ["restart0010"]
    fname = []
    for n in fnums:
        fname.append(path + n + ".sdf")

    for i in range(len(fname)):
        sdfdata = sdf.read(fname[i])

        grid = sdfdata.__dict__["Grid_Grid"].data

        x = grid[0]
        y = grid[1]

        dx = grid[0][1] - grid[0][0]
        dy = grid[1][1] - grid[1][0]

        p_pos = sdfdata.__dict__[get_varname("Grid_Particles", species)].data
        p_list = list(zip(p_pos[0], p_pos[1]))

        p_list[:] = [p for p in p_list if in_domain(p, x, y)]

        w = sdfdata.__dict__[get_varname("Particles_Weight", species)].data
        v = sdfdata.__dict__[get_varname("Particles_Px", species)].data

        grid = first_order_weight_2d(x, y, dx, dy, p_list, values=v, weight=w)

    # p = (2 - (-2)) * np.random.rand(500,2) - 2#np.array([[0., 0.],[0.5,0.5],[-1,-1],[-2,-2],[-1.8,-1.5]])

    # grid = first_order_weight_2d(x, y, dx, dy, p, values=None, weight=None)

    cmap = shiftedColorMap(cm.coolwarm, np.amin(grid), np.amax(grid))
    plt.figure()
    plt.pcolormesh(grid, cmap=cmap)
    cbar = plt.colorbar()
    plt.savefig('w.png')
    return


if __name__ == "__main__":
    main()
