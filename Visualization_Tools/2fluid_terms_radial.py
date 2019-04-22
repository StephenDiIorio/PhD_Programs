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
from matplotlib.ticker import FuncFormatter

import sdf
from CoordinateTransforms import cart2polar, reproject_image_into_polar
from PlottingTools import get_si_prefix, get_var_range
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


def __Temp_Grad(sdfdata, axis, species=None):
    temp_key = get_varname("Derived_Temperature", species=species)
    num_dens_key = get_varname("Derived_Number_Density", species=species)

    # Can also get grid for each enty via:
    # sdfdata.__dict__[varname].grid.data
    grid = sdfdata.__dict__["Grid_Grid_mid"].data

    # get the spatial spacing in grid to perform gradient (deltax,y)
    dx = grid[0][1] - grid[0][0]
    dy = grid[1][1] - grid[1][0]

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
                             where=num_dens.data != 0.0)


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
                             where=e_num_dens.data != 0.0)


def __in_domain(p, x, y):
    return p[0] >= np.amin(x) and p[0] <= np.amax(x) and p[1] >= np.amin(y) and p[1] <= np.amax(y)


def __map_particle_to_grid(sdfdata, p_var_name, species=None):
    grid = sdfdata.__dict__["Grid_Grid_mid"].data

    x = grid[0]
    y = grid[1]

    dx = grid[0][1] - grid[0][0]
    dy = grid[1][1] - grid[1][0]

    p_pos = sdfdata.__dict__[get_varname("Grid_Particles", species)].data
    p_list = list(zip(p_pos[0], p_pos[1]))

    w = sdfdata.__dict__[get_varname("Particles_Weight", species)].data

    p_var = sdfdata.__dict__[p_var_name].data

    return first_order_weight_2d(x, y, dx, dy, p_list,
                                 values=p_var, weight=w)


def Ideal_MHD_Field(sdfdata, axis, species=None):
    vel_x_key = get_varname("Particles_Vx", species)
    vel_y_key = get_varname("Particles_Vy", species)
    vel_z_key = get_varname("Particles_Vz", species)
    vx = __map_particle_to_grid(sdfdata, vel_x_key, species)
    vy = __map_particle_to_grid(sdfdata, vel_y_key, species)
    vz = __map_particle_to_grid(sdfdata, vel_z_key, species)

    mag_x_key = "Magnetic_Field_Bx"
    mag_y_key = "Magnetic_Field_By"
    mag_z_key = "Magnetic_Field_Bz"

    c = __Cross_Product(vx,
                        vy,
                        vz,
                        sdfdata.__dict__[mag_x_key].data,
                        sdfdata.__dict__[mag_y_key].data,
                        sdfdata.__dict__[mag_z_key].data)
    return sdfdata.__dict__[mag_x_key].grid, c[axis]


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
                             where=temp.data != 0.0)


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
                                         where=e_num_dens.data != 0.0))


def Generalized_Ohm(sdfdata, axis, species=None):
    grid = sdfdata.__dict__["Grid_Grid_mid"].data
    # grid = sdfdata.__dict__["Grid_Grid"].data

    x = grid[0]
    y = grid[1]

    dx = grid[0][1] - grid[0][0]
    dy = grid[1][1] - grid[1][0]

    p_pos = sdfdata.__dict__[get_varname("Grid_Particles", species)].data
    vx = sdfdata.__dict__[get_varname("Particles_Vx", species)].data
    vy = sdfdata.__dict__[get_varname("Particles_Vy", species)].data
    vz = sdfdata.__dict__[get_varname("Particles_Vz", species)].data
    w = sdfdata.__dict__[get_varname("Particles_Weight", species)].data

    v_3 = (np.sqrt(vx**2 + vy**2 + vz**2))**3
    v_5 = (np.sqrt(vx**2 + vy**2 + vz**2))**5

    p_list = list(zip(p_pos[0], p_pos[1]))


    v_3_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_3)
    v_5_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_5)

    const = -sc.m_e / (6 * sc.e)
    grad_num = np.gradient(v_5_dist, dx, dy)[axis]
    # grad_num_y = np.gradient(v_5_dist, dx, dy)[1]

    return const * np.divide(grad_num, v_3_dist, where=v_3_dist!=0.0)
    # term_y = const * np.divide(grad_num_y, v_3_dist, where=v_3_dist!=0.0)


def Generalized_Ohm_radavg(sdfdata, species=None):
    # grid = sdfdata.__dict__["Grid_Grid_mid"].data
    grid = sdfdata.__dict__["Grid_Grid"].data

    x = grid[0]
    y = grid[1]

    dx = grid[0][1] - grid[0][0]
    dy = grid[1][1] - grid[1][0]

    p_pos = sdfdata.__dict__[get_varname("Grid_Particles", species)].data
    vx = sdfdata.__dict__[get_varname("Particles_Vx", species)].data
    vy = sdfdata.__dict__[get_varname("Particles_Vy", species)].data
    vz = sdfdata.__dict__[get_varname("Particles_Vz", species)].data
    w = sdfdata.__dict__[get_varname("Particles_Weight", species)].data

    v_3 = (np.sqrt(vx**2 + vy**2 + vz**2))**3
    v_5 = (np.sqrt(vx**2 + vy**2 + vz**2))**5

    p_list = list(zip(p_pos[0], p_pos[1]))


    v_3_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_3)
    v_5_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_5)

    # ax = plt.subplot(1, 2, 1)
    # im = ax.pcolormesh(v_3_dist, cmap=cm.coolwarm)

    # ax = plt.subplot(1, 2, 2)
    # im = ax.pcolormesh(v_5_dist,  cmap=cm.coolwarm)

    # plt.savefig('v5.png', dpi=600, bbox_inches="tight")

    avg_3, dummy_r, dummy_t = reproject_image_into_polar(v_3_dist, origin=None, Jacobian=False, dr=1, dt=None)
    o_ma = np.ma.masked_equal(avg_3, 0.)
    avg_3 = np.average(o_ma, axis=1)

    avg_5, dummy_r, dummy_t = reproject_image_into_polar(v_5_dist, origin=None, Jacobian=False, dr=1, dt=None)
    o_ma = np.ma.masked_equal(avg_5, 0.)
    avg_5 = np.average(o_ma, axis=1)

    const = -sc.m_e / (6 * sc.e)
    grad_num = np.gradient(avg_5, dx) #TODO: what to use for grid spacing radially

    return const * np.divide(grad_num, avg_3, where=avg_3!=0.0)
    # term_y = const * np.divide(grad_num_y, v_3_dist, where=v_3_dist!=0.0)


def __radial_average(x_grid, y_grid, field_x, field_y):
    X, Y = np.meshgrid(x_grid, y_grid)
    R, T = cart2polar(X, Y)

    er = np.multiply(field_x, np.cos(T)) + np.multiply(field_y, np.sin(T))
    dummy_et = np.multiply(field_y, np.cos(T)) - np.multiply(field_x, np.sin(T))

    o, dummy_r, dummy_t = reproject_image_into_polar(er, origin=None, Jacobian=False, dr=1, dt=None)

    o_ma = np.ma.masked_equal(o, 0.)
    avg = np.average(o_ma, axis=1)

    return avg, R, T


def higher_order(sdfdata, species=None):
    # grid = sdfdata.__dict__["Grid_Grid_mid"].data
    grid = sdfdata.__dict__["Grid_Grid"].data

    x = grid[0]
    y = grid[1]

    dx = grid[0][1] - grid[0][0]
    dy = grid[1][1] - grid[1][0]

    p_pos = sdfdata.__dict__[get_varname("Grid_Particles", species)].data
    vx = sdfdata.__dict__[get_varname("Particles_Vx", species)].data
    vy = sdfdata.__dict__[get_varname("Particles_Vy", species)].data
    vz = sdfdata.__dict__[get_varname("Particles_Vz", species)].data
    w = sdfdata.__dict__[get_varname("Particles_Weight", species)].data

    p_list = list(zip(p_pos[0], p_pos[1]))


    limit = 1e7

    v_r = np.sqrt(vx**2 + vy**2)# + vz**2) #NOTE: might need to get rid of vz. can then use cart2polar for v_r and v_t
    v_r_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_r)
    avg_r, r, dummy_t = reproject_image_into_polar(v_r_dist, origin=None, Jacobian=False, dr=1, dt=None)
    o_ma = np.ma.masked_equal(avg_r, 0.)
    avg_r = np.average(o_ma, axis=1)

    fig1, ax1 = plt.subplots(2, sharex=True)
    ax1[0].plot(avg_r)
    ax1[1].plot(np.clip(avg_r, -limit, limit))
    fig1.savefig('vr.png', dpi=600, bbox_inches="tight")


    v_t = np.arctan2(vx, vy)
    v_t_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_t)
    avg_t, r, t = reproject_image_into_polar(v_t_dist, origin=None, Jacobian=False, dr=1, dt=None)
    o_ma = np.ma.masked_equal(avg_t, 0.)
    avg_t = np.average(o_ma, axis=1)

    fig2, ax2 = plt.subplots(2, sharex=True)
    ax2[0].plot(avg_t)
    ax2[1].plot(np.clip(avg_t, -limit, limit))
    fig2.savefig('vt.png', dpi=600, bbox_inches="tight")


    v_3 = (np.sqrt(vx**2 + vy**2 + vz**2))**3
    v_3_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_3)
    avg_3, r, t = reproject_image_into_polar(v_3_dist, origin=None, Jacobian=False, dr=1, dt=None)
    o_ma = np.ma.masked_equal(avg_3, 0.)
    avg_3 = np.average(o_ma, axis=1)

    fig3, ax3 = plt.subplots(2, sharex=True)
    ax3[0].plot(avg_3)
    ax3[1].plot(np.clip(avg_3, -limit, limit))
    fig3.savefig('avg_v3.png', dpi=600, bbox_inches="tight")


    r = np.linspace(0.0, np.max(np.sqrt(x**2 + y**2)), num=avg_r.size)
    print(r)

    num_1 = np.multiply(r, avg_3)
    num_1 = np.multiply(avg_r, num_1)
    num_1 = np.multiply(avg_r, num_1)

    fig4, ax4 = plt.subplots(2, sharex=True)
    ax4[0].plot(num_1)
    ax4[1].plot(np.clip(num_1, -limit, limit))
    fig4.savefig('num1_before_grad', dpi=600, bbox_inches="tight")

    num_1 = np.gradient(num_1, dx) #TODO: what to use for grid spacing radially

    fig5, ax5 = plt.subplots(2, sharex=True)
    ax5[0].plot(num_1)
    ax5[1].plot(np.clip(num_1, -limit, limit))
    fig5.savefig('num1.png', dpi=600, bbox_inches="tight")


    num_2 = np.multiply(avg_r, avg_t)
    num_2 = np.multiply(avg_3, num_2)

    fig6, ax6 = plt.subplots(2, sharex=True)
    ax6[0].plot(num_2)
    ax6[1].plot(np.clip(num_2, -limit, limit))
    fig6.savefig('num2_before_grad.png', dpi=600, bbox_inches="tight")

    num_2 = np.gradient(num_2, max(x.size, y.size)) #TODO: what to use for grid spacing theta

    fig7, ax7 = plt.subplots(2, sharex=True)
    ax7[0].plot(num_2)
    ax7[1].plot(np.clip(num_2, -limit, limit))
    fig7.savefig('num2.png', dpi=600, bbox_inches="tight")


    # num = np.divide(num_1 + num_2, r, where=r!=0.0)
    num = np.divide(num_1, r, where=r!=0.0)
    fig8, ax8 = plt.subplots(2, sharex=True)
    ax8[0].plot(num)
    ax8[1].plot(np.clip(num, -limit, limit))
    fig8.savefig('num.png', dpi=600, bbox_inches="tight")

    const = -sc.m_e / (2 * sc.e)
    final = const * np.divide(num, avg_3, where=avg_3!=0.0)
    fig9, ax9 = plt.subplots(2, sharex=True)
    ax9[0].plot(final)
    ax9[1].plot(np.clip(final, -limit, limit))
    fig9.savefig('final.png', dpi=600, bbox_inches="tight")

    return final#const * np.divide(num, avg_3, where=avg_3!=0.0)


def cart_higher_order_x(sdfdata, species=None):
    # grid = sdfdata.__dict__["Grid_Grid_mid"].data
    grid = sdfdata.__dict__["Grid_Grid"].data

    x = grid[0]
    y = grid[1]

    dx = grid[0][1] - grid[0][0]
    dy = grid[1][1] - grid[1][0]

    p_pos = sdfdata.__dict__[get_varname("Grid_Particles", species)].data
    vx = sdfdata.__dict__[get_varname("Particles_Vx", species)].data
    vy = sdfdata.__dict__[get_varname("Particles_Vy", species)].data
    vz = sdfdata.__dict__[get_varname("Particles_Vz", species)].data
    w = sdfdata.__dict__[get_varname("Particles_Weight", species)].data

    v_3 = (np.sqrt(vx**2 + vy**2 + vz**2))**3

    p_list = list(zip(p_pos[0], p_pos[1]))


    v_3_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_3)
    vx_dist  = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=vx)
    vy_dist  = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=vy)

    fig1, ax1 = plt.subplots()
    ax1.pcolor(v_3_dist)
    fig1.colorbar()
    fig1.savefig('v3.png', dpi=600, bbox_inches="tight")
    fig2, ax2 = plt.subplots()
    ax2.pcolor(vx_dist)
    fig2.colorbar()
    fig2.savefig('vx.png', dpi=600, bbox_inches="tight")
    fig3, ax3 = plt.subplots()
    ax3.pcolor(vy_dist)
    fig3.colorbar()
    fig3.savefig('vy.png', dpi=600, bbox_inches="tight")

    term1 = np.multiply(vx_dist, vx_dist)
    term1 = np.multiply(term1, v_3_dist)
    fig4, ax4 = plt.subplots()
    ax4.pcolor(term1)
    fig4.colorbar()
    fig4.savefig('term1_cartx_before_grad.png', dpi=600, bbox_inches="tight")

    term1 = np.gradient(term1, dx, axis=0)
    fig5, ax5 = plt.subplots()
    ax5.pcolor(term1)
    fig5.colorbar()
    fig5.savefig('term1_cartx.png', dpi=600, bbox_inches="tight")

    term2 = np.multiply(vx_dist, vy_dist)
    term2 = np.multiply(term2, v_3_dist)
    fig6, ax6 = plt.subplots()
    ax6.pcolor(term2)
    fig6.colorbar()
    fig6.savefig('term2_cartx_before_grad.png', dpi=600, bbox_inches="tight")

    term2 = np.gradient(term2, dy, axis=1)
    fig7, ax7 = plt.subplots()
    ax7.pcolor(term2)
    fig7.colorbar()
    fig7.savefig('term2_cartx.png', dpi=600, bbox_inches="tight")

    div_num = term1 + term2

    const = -sc.m_e / (2 * sc.e)
    final = const * np.divide(div_num, v_3_dist, where=v_3_dist!=0.0)
    fig8, ax8 = plt.subplots()
    ax8.pcolor(final)
    fig8.colorbar()
    fig8.savefig('final_cartx.png', dpi=600, bbox_inches="tight")

    return final


def cart_higher_order_y(sdfdata, species=None):
    # grid = sdfdata.__dict__["Grid_Grid_mid"].data
    grid = sdfdata.__dict__["Grid_Grid"].data

    x = grid[0]
    y = grid[1]

    dx = grid[0][1] - grid[0][0]
    dy = grid[1][1] - grid[1][0]

    p_pos = sdfdata.__dict__[get_varname("Grid_Particles", species)].data
    vx = sdfdata.__dict__[get_varname("Particles_Vx", species)].data
    vy = sdfdata.__dict__[get_varname("Particles_Vy", species)].data
    vz = sdfdata.__dict__[get_varname("Particles_Vz", species)].data
    w = sdfdata.__dict__[get_varname("Particles_Weight", species)].data

    v_3 = (np.sqrt(vx**2 + vy**2 + vz**2))**3

    p_list = list(zip(p_pos[0], p_pos[1]))


    v_3_dist = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=v_3)
    vx_dist  = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=vx)
    vy_dist  = first_order_weight_2d(x, y, dx, dy, p_list, weight=w, values=vy)

    term1 = np.multiply(vy_dist, vx_dist)
    term1 = np.multiply(term1, v_3_dist)
    fig4, ax4 = plt.subplots()
    ax4.pcolor(term1)
    fig4.colorbar()
    fig4.savefig('term1_carty_before_grad.png', dpi=600, bbox_inches="tight")

    term1 = np.gradient(term1, dx, axis=0)
    fig5, ax5 = plt.subplots()
    ax5.pcolor(term1)
    fig5.colorbar()
    fig5.savefig('term1_carty.png', dpi=600, bbox_inches="tight")

    term2 = np.multiply(vy_dist, vy_dist)
    term2 = np.multiply(term2, v_3_dist)
    fig6, ax6 = plt.subplots()
    ax6.pcolor(term2)
    fig6.colorbar()
    fig6.savefig('term2_carty_before_grad.png', dpi=600, bbox_inches="tight")

    term2 = np.gradient(term2, dy, axis=1)
    fig7, ax7 = plt.subplots()
    ax7.pcolor(term2)
    fig7.colorbar()
    fig7.savefig('term2_carty.png', dpi=600, bbox_inches="tight")

    div_num = term1 + term2

    const = -sc.m_e / (2 * sc.e)
    final = const * np.divide(div_num, v_3_dist, where=v_3_dist!=0.0)
    fig8, ax8 = plt.subplots()
    ax8.pcolor(final)
    fig8.colorbar()
    fig8.savefig('final_carty.png', dpi=600, bbox_inches="tight")

    return final


def main():
    species = "Electron"
    path = "/scratch/lsa_flux/diiorios/2d_run/"
    fnums = ["0100", "0200", "0400"]
    fname = []
    for n in fnums:
        fname.append(path + n + ".sdf")

    x_axis_num = 0
    y_axis_num = 1
    # z_axis_num = 2

    fig, axarr = plt.subplots(len(fname), sharex=True)#, sharey=True)
    fig.set_facecolor("w")
    # only have a single file to plot, so we so this little hack since
    # axarr does not come out in an array in this case
    if not isinstance(axarr, np.ndarray):
        axarr = np.array([axarr])

    # axarr[0].set_title(species + " files " + str(fnums))
    axarr[0].set_title('Contribution to E Field')

    limit = 7E9

    for i in range(len(fname)):
        sdfdata = sdf.read(fname[i])

        e_var_x = sdfdata.Electric_Field_Ex
        e_var_y = sdfdata.Electric_Field_Ey
        x_grid = e_var_x.grid_mid.data[0]
        y_grid = e_var_y.grid_mid.data[1]
        x_grid_cent = e_var_x.grid.data[0]
        y_grid_cent = e_var_y.grid.data[1]

        e_avg, dummy_e_r, dummy_e_t = __radial_average(x_grid, y_grid, e_var_x.data, e_var_y.data)

        dummy_temp_grid_x, temp_data_x = Temp_Field(sdfdata,
                                                    x_axis_num,
                                                    species='Electron')
        dummy_temp_grid_y, temp_data_y = Temp_Field(sdfdata,
                                                    y_axis_num,
                                                    species='Electron')
        temp_avg, dummy_temp_r, dummy_temp_t = __radial_average(x_grid, y_grid,
                                                                temp_data_x,
                                                                temp_data_y)

        dummy_res_mhd_grid_x, res_mhd_data_x = Resistive_MHD_Field(sdfdata,
                                                                   x_axis_num)
        dummy_res_mhd_grid_y, res_mhd_data_y = Resistive_MHD_Field(sdfdata,
                                                                   y_axis_num)
        rmhd_avg, dummy_rmhd_r, dummy_rmhd_t = __radial_average(x_grid, y_grid,
                                                                res_mhd_data_x,
                                                                res_mhd_data_y)

        dummy_hall_grid_x, hall_data_x = Hall_Field(sdfdata, x_axis_num)
        dummy_hall_grid_y, hall_data_y = Hall_Field(sdfdata, y_axis_num)
        hall_avg, dummy_hall_r, dummy_hall_t = __radial_average(x_grid, y_grid, hall_data_x, hall_data_y)

        gen_data_x = Generalized_Ohm(sdfdata,
                                     x_axis_num,
                                     species='Electron')
        gen_data_y = Generalized_Ohm(sdfdata,
                                     y_axis_num,
                                     species='Electron')
        gen_avg, dummy_gen_r, dummy_gen_t = __radial_average(x_grid, y_grid, gen_data_x, gen_data_y)

        gen_avg2 = Generalized_Ohm_radavg(sdfdata,
                                          species='Electron')

        high_order = higher_order(sdfdata, species='Electron')

        cart_high_x = cart_higher_order_x(sdfdata, specied='Electron')
        cart_high_y = cart_higher_order_y(sdfdata, species='Electron')
        cart_high_avg, dummy_high_r, dummy_high_t = __radial_average(x_grid_cent, y_grid_cent, cart_high_x, cart_high_y)

        # ideal_mhd_grid_x, ideal_mhd_data_x = Ideal_MHD_Field(sdfdata,
        #                                                      x_axis_num,
        #                                                      species='Electron')
        # ideal_mhd_grid_y, ideal_mhd_data_y = Ideal_MHD_Field(sdfdata,
        #                                                      y_axis_num,
        #                                                      species='Electron')
        # imhd_avg, imhd_r, imhd_t = __radial_average(x_grid, y_grid, ideal_mhd_data_x, ideal_mhd_data_y)


        l1, = axarr[i].plot(#e_r,
                            e_avg,
                            'k-',
                            label='Simulation')
        l2, = axarr[i].plot(#temp_r,
                            np.clip(temp_avg, -limit, limit),
                            'r--',
                            label='Thermal')
        # l3, = axarr[i].plot(#rmhd_r,
        #                     np.clip(rmhd_avg, -limit, limit),
        #                     'b-.',
        #                     label='Resistive MHD')
        # l4, = axarr[i].plot(#hall_r,
        #                     np.clip(hall_avg, -limit, limit),
        #                     'g:',
        #                     label='Hall Term')
        # l5, = axarr[i].plot(#imhd_r,
        #                    np.clip(imhd_avg,
        #                            -limit, limit),
        #                    'm-.',
        #                    label='Ideal MHD Term')
        l6, = axarr[i].plot(np.clip(gen_avg, -limit, limit),
                            'm-.',
                            label='Generalized Ohm',
                            alpha=0.2)
        l7, = axarr[i].plot(np.clip(gen_avg2, -limit, limit),
                            'c-.',
                            label='Gen Ohm Avg 1st')
        l8, = axarr[i].plot(np.clip(high_order, -limit, limit),
                            'g:',
                            label='Higher Order')
        l9, = axarr[i].plot(np.clip(cart_high_avg, -limit, limit),
                            'b:',
                            label='Cart Higher Order')

    ls = [l1, l2, l6, l7, l8, l9]  #l3, l4, l5]
    labels = [l.get_label() for l in ls]
    lgd = fig.legend(ls, labels, bbox_to_anchor=(1.05, 1.0), loc=1)
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    # plt.xlim(-6.5e-6, 6.5e-6)
    # plt.ylim(-limit, limit)

    # xmult, xsym = get_si_prefix(np.max(e_r) - np.min(e_r))
    ymult, ysym = get_si_prefix(limit - (-limit))
    # axarr[0].xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * xmult)))
    axarr[0].yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult)))
    axarr[1].yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult)))
    axarr[2].yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult)))

    axarr[0].set_ylim([-limit, limit])
    axarr[1].set_ylim([-limit, limit])
    axarr[2].set_ylim([-limit, limit])

    # plt.xlabel('x' + ' $(' + xsym + 'm)$')
    plt.xlabel('x')
    axarr[1].set_ylabel('Radial Electric Field' + ' $(' + ysym + 'V/m)$')
    # plt.show()

    plt.savefig('ohm_rad.png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches = "tight", pad_inches=0.2)


if __name__ == "__main__":
    main()
