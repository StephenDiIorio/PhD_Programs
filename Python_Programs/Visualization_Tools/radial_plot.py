from matplotlib import use
use('Agg')

import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

from CoordinateTransforms import cart2polar, reproject_image_into_polar
import sdf


def main():
    path = "/Users/stephendiiorio/Desktop/"
    fnums = ["restart0010"]
    fname = []
    for n in fnums:
        fname.append(path + n + ".sdf")

    for i in range(len(fname)):
        sdfdata = sdf.read(fname[i])

        ex = sdfdata.__dict__['Electric_Field_Ex'].data
        ey = sdfdata.__dict__['Electric_Field_Ey'].data

        x_grid = sdfdata.__dict__['Electric_Field_Ex'].grid_mid.data[0]
        y_grid = sdfdata.__dict__['Electric_Field_Ex'].grid_mid.data[1]

        # See https://www.mathworks.com/matlabcentral/answers/95758-how-can-i-convert-vector-fields-from-cartesian-to-polar-coordinate-system-in-matlab-7-10-r2010a for a discussion on how to create a radial vector field from a cartesian one.
        X, Y = np.meshgrid(x_grid, y_grid)

        R, T = cart2polar(X, Y)
        print(np.max(R), np.min(R))

        er = np.multiply(ex, np.cos(T)) + np.multiply(ey, np.sin(T))
        et = np.multiply(ey, np.cos(T)) - np.multiply(ex, np.sin(T))

        # See https://stackoverflow.com/questions/2164570/reprojecting-polar-to-cartesian-grid and the python package Ansel, where these functions are taken from, to see how to transform a cartesian image into polar coordinates.
        o, r, t = reproject_image_into_polar(er, origin=None, Jacobian=False, dr=1, dt=None)

        o_ma = np.ma.masked_equal(o, 0.)
        avg = np.average(o_ma, axis=1)

        fig, ax = plt.subplots()
        # im = ax.pcolormesh(o,  cmap=cm.coolwarm,  vmin=-8e8, vmax=8e8)
        # cb = fig.colorbar(im)
        im = ax.plot(avg)
        plt.ylim(-8e8, 8e8)
        plt.savefig('asdf.png', dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    main()
