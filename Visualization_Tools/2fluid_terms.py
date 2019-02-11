import glob
import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
from matplotlib import use

import sdf

use('Agg')


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


def __Temp_Grad(sdfdata, axis, species=None):
    temp_key = get_varname("Derived_Temperature", species=species)
    num_dens_key = get_varname("Derived_Number_Density", species=species)

    # Can also get grid for each enty via:
    # sdfdata.__dict__[varname].grid.data
    grid = sdfdata.__dict__["Grid_Grid"]

    # get the spatial spacing in grid to perform gradient (deltax,y)
    dx = abs(grid.data[0][0]) / len(grid.data[0])
    dy = abs(grid.data[1][0]) / len(grid.data[1])

    temp = sdfdata.__dict__[temp_key]
    num_dens = sdfdata.__dict__[num_dens_key]
    return np.gradient(np.multiply(temp.data, num_dens.data), dx, dy)[axis]


# TODO: Change e_num_dens to match the rest of the function styles
# TODO: See which num density you should use to calculate this
def Temp_Field(sdfdata, axis, species=None):
    g = __Temp_Grad(sdfdata, axis, species=species)
    const = -sc.k / sc.e

    num_dens_key = get_varname("Derived_Number_Density", species=species)
    num_dens = sdfdata.__dict__[num_dens_key]

    return num_dens.grid, \
           const * np.divide(g, num_dens.data,
                             out=np.zeros_like(g),
                             where=num_dens.data != 0)


def __Cross_Product(x_comp1, y_comp1, z_comp1, x_comp2, y_comp2, z_comp2):
    cross_x = np.multiply(y_comp1, z_comp2) - np.multiply(z_comp1, y_comp2)
    cross_y = -(np.multiply(x_comp1, z_comp2) - np.multiply(z_comp1, x_comp2))
    cross_z = np.multiply(x_comp1, y_comp2) - np.multiply(y_comp1, x_comp2)
    return cross_x, cross_y, cross_z


# TODO: Check if currents should be for a particular species
def Hall_Field(sdfdata, axis):
    const = 1. / sc.e

    current_x_key = "Current_Jx"
    current_y_key = "Current_Jy"
    current_z_key = "Current_Jz"
    mag_x_key = "Magnetic_Field_Bx"
    mag_y_key = "Magnetic_Field_By"
    mag_z_key = "Magnetic_Field_Bz"

    c = __Cross_Product(sdfdata.__dict__[current_x_key].data,
                        sdfdata.__dict__[current_y_key].data,
                        sdfdata.__dict__[current_z_key].data,
                        sdfdata.__dict__[mag_x_key].data,
                        sdfdata.__dict__[mag_y_key].data,
                        sdfdata.__dict__[mag_z_key].data)

    e_num_dens_key = get_varname("Derived_Number_Density", species="Electron")
    e_num_dens = sdfdata.__dict__[e_num_dens_key]

    return e_num_dens.grid, \
           const * np.divide(c[axis], e_num_dens.data,
                             out=np.zeros_like(e_num_dens.data),
                             where=e_num_dens.data != 0)


def Ideal_MHD_Field(sdfdata, axis, species=None):
    vel_x_key = get_varname("Particles_Velocity_Vx", species)
    vel_y_key = get_varname("Particles_Velocity_Vy", species)
    vel_z_key = get_varname("Particles_Velocity_Vz", species)
    mag_x_key = "Magnetic_Field_Bx"
    mag_y_key = "Magnetic_Field_By"
    mag_z_key = "Magnetic_Field_Bz"

    c = __Cross_Product(sdfdata.__dict__[vel_x_key].data,
                        sdfdata.__dict__[vel_y_key].data,
                        sdfdata.__dict__[vel_z_key].data,
                        sdfdata.__dict__[mag_x_key].data,
                        sdfdata.__dict__[mag_y_key].data,
                        sdfdata.__dict__[mag_z_key].data)
    return sdfdata.__dict__[vel_x_key].grid, c[axis]


# TODO: check if collision rate is good approximation
def __Coll_Rate(sdfdata):
    temp_key = get_varname("Derived_Temperature", species="Electron")
    e_num_dens_key = get_varname("Derived_Number_Density", species="Electron")

    temp = sdfdata.__dict__[temp_key]
    e_num_dens = sdfdata.__dict__[e_num_dens_key]

    # calculate rate assuming that electron-ion rate is approx electron
    # collision rate, using equation from Plasma Formulary (pg 28)
    # take coulomb log to be 20
    # convert to cm^-3
    # use conversion rate of 1eV=kb*11600K
    const = 2.91E-6 * 1E-6 * 20
    return const * np.divide(e_num_dens.data,
                             np.power(temp.data / 11600., 1.5),
                             where=temp.data != 0)


# TODO: Make sure that you look at what the collision constant is here
def Resistive_MHD_Field(sdfdata, axis):
    const = -sc.m_e / sc.e / sc.e

    e_num_dens_key = get_varname("Derived_Number_Density", species="Electron")
    e_num_dens = sdfdata.__dict__[e_num_dens_key]

    if axis is 0:
        current_key = "Current_Jx"
    elif axis is 1:
        current_key = "Current_Jy"
    elif axis is 2:
        current_key = "Current_Jz"
    else:
        raise ValueError('Undefined axis {}.'.format(axis))

    current = sdfdata.__dict__[current_key]

    return e_num_dens.grid, \
           const * np.multiply(__Coll_Rate(sdfdata),
                               np.divide(current.data, e_num_dens.data,
                                         out=np.zeros_like(current.data),
                                         where=e_num_dens.data != 0))


def main():
    species = "Argon"
    path = "/Users/stephendiiorio/Documents/Research/"
    fnums = ["0514", "0515", "0516"]
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

    limit = 10E9

    for i in range(len(fname)):
        sdfdata = sdf.read(fname[i])
        # sh.list_variables(sdfdata)

        e_var = sdfdata.Electric_Field_Ey
        midpoint = int(len(e_var.data[0]) / 2)

        temp_grid, temp_data = Temp_Field(sdfdata,
                                          y_axis_num,
                                          species='Electron')
        res_mhd_grid, res_mhd_data = Resistive_MHD_Field(sdfdata, y_axis_num)

        x_axis = e_var.grid_mid.data[1]

        # fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
        l1, = axarr[i].plot(x_axis,
                            e_var.data[midpoint],
                            'k-',
                            label='Ey')
        # plt.plot(x_axis, e_var.data[midpoint], 'k-', label='Ey')
        l2, = axarr[i].plot(x_axis,
                            np.clip(temp_data[midpoint], -limit, limit),
                            'r--',
                            label='Temp_Ey')
        # plt.plot(x_axis, temp_data[midpoint], 'r--', label='Temp_Ey')
        l3, = axarr[i].plot(x_axis,
                            np.clip(res_mhd_data[midpoint], -limit, limit),
                            'b-.',
                            label='Res_MHD_Ey')

    # fig.legend([l1, l2, l3],
    #            ['Ey', 'Temp_Ey', 'Res_MHD_Ey'],
    #            loc='upper right')
    ls = [l1, l2, l3]
    labels = [l.get_label() for l in ls]
    fig.legend(ls, labels, loc='upper right')
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    # plt.show()

    plt.savefig('test.png')


if __name__ == "__main__":
    main()
