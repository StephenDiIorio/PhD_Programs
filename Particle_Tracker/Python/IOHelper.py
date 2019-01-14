import sys
sys.path.insert(0, '../../Utilities/')

from matplotlib import use
use('Agg')
import numpy as np
import PlottingTools as pt
from scipy.interpolate import griddata
# import seaborn as sns

try:
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    try:
        print("Trying mplot3d import workaround")
        # Workaround for broken macOS installation
        import os
        import sys
        import matplotlib
        sys.path.append(os.path.join(matplotlib.__path__[0],
                                     '..', 'mpl_toolkits'))
        from mplot3d import Axes3D
    except ImportError as e:
        print(e)
        pass
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=15)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize


def print_field_to_file(field, save_x=True, save_y=True, save_z=True):
    # print field to file to view in mathematica
    x_dim = 0
    y_dim = 1
    z_dim = 2

    if save_x:
        with open('x_field.dat', 'wb') as f:
            for index, val in np.ndenumerate(field.get_field()[x_dim]):
                to_write = [index[0], index[1], index[2], val]
                f.write(' '.join([str(float(num)) for num in to_write]))
                f.write('\n')
    if save_y:
        with open('y_field.dat', 'wb') as f:
            for index, val in np.ndenumerate(field.get_field()[y_dim]):
                to_write = [index[0], index[1], index[2], val]
                f.write(' '.join([str(float(num)) for num in to_write]))
                f.write('\n')
    if save_z:
        with open('z_field.dat', 'wb') as f:
            for index, val in np.ndenumerate(field.get_field()[z_dim]):
                to_write = [index[0], index[1], index[2], val]
                f.write(' '.join([str(float(num)) for num in to_write]))
                f.write('\n')
    return


def save_data_array(fname, *args):
    np.savez_compressed(fname, *args)
    return


def close_figures():
    plt.close('all')
    return


def plot_trajectories(data, fname=None):
    traj_fig = plt.figure()
    traj_ax = traj_fig.gca(projection='3d')
    traj_fig.suptitle('Particle Tracking Simulation', fontsize=18)
    for i in range(len(data)):
        x = []
        y = []
        z = []
        for j in range(len(data[i])):
            x.append(data[i][j][0])
            y.append(data[i][j][1])
            z.append(data[i][j][2])
        # traj_ax.plot(data[i, :, 0], data[i, :, 1], data[i, :, 2])
        traj_ax.plot(x, y, z)
    traj_ax.set_xlabel('X', fontsize=18)
    traj_ax.set_ylabel('Y', fontsize=18)
    traj_ax.set_zlabel('Z', fontsize=18)
    # traj_ax.set_xlim(0.0, 8.0)
    # traj_ax.set_ylim(0, ny * step[1])
    # traj_ax.set_zlim(4.8, 5.2)
    # # Hide tick numbers in plot
    # traj_ax.set_xticklabels([])
    # traj_ax.set_yticklabels([])
    # traj_ax.set_zticklabels([])

    if fname is not None:
        traj_fig.savefig(fname + '.png', dpi=600, bbox_inches="tight")
    return


def plot_heatmap(x1, x2, r, fname=None):
    dens_fig1 = plt.figure()
    dens_ax1 = dens_fig1.add_subplot(111)
    dens_fig2 = plt.figure()
    dens_ax2 = dens_fig2.add_subplot(111)

    # sns.jointplot(x=x1, y=x2, kind='hex')

    heatmap, yedges, zedges = np.histogram2d(x1, x2, bins=100)
    grid_x1, grid_x2 = np.meshgrid(yedges, zedges)
    dens_ax1.pcolormesh(grid_x1, grid_x2, heatmap.T)

    heatmap_r, yedges_r, zedges_r = np.histogram2d(x1, x2, range=r, bins=100)
    grid_x1_r, grid_x2_r = np.meshgrid(yedges_r, zedges_r)
    dens_ax2.pcolormesh(grid_x1_r, grid_x2_r, heatmap_r.T)

    # ymin, ymax = get_var_range(y)
    # ymult, ysym = get_si_prefix(ymax - ymin)
    # dens_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult)))

    # zmin, zmax = get_var_range(z)
    # zmult, zsym = get_si_prefix(zmax - zmin)
    # dens_ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * zmult)))

    if fname is not None:
        dens_fig1.savefig(fname + '.png', dpi=600, bbox_inches="tight")
        save_data_array(fname + '_full', heatmap.T, yedges, zedges)
        dens_fig2.savefig('r_' + fname + '.png', dpi=600, bbox_inches="tight")
        save_data_array(fname + '_range', heatmap_r.T, yedges_r, zedges_r)
    return heatmap, heatmap_r


def plot_rel_diff(radii, time, particle_count, base_count, fname=None):
    rel_fig = plt.figure()
    rel_ax = rel_fig.add_subplot(111)

    particle_count = np.array(particle_count)
    base_count = np.array(base_count)

    l1, = rel_ax.plot(time, particle_count, 'k-', label='Particle Count')
    l2, = rel_ax.plot(time, particle_count - base_count, 'r--', label='Difference')

    ls = [l1, l2]
    labels = [l.get_label() for l in ls]
    lgd = rel_ax.legend(ls, labels, loc=2)
    # lgd.get_frame().set_alpha(1)

    rel_ax.set_xlim(0., 100.)
    rel_ax.set_xlabel('Time $(ps)$')
    rel_ax.set_ylabel('Counts (arb. units)')
    rel_ax.set_title('Particle Tracking Counts')

    if fname is not None:
        plt.savefig(fname + '.png', dpi=600, bbox_extra_artists=(lgd,), bbox_inches="tight", pad_inches=0.2)
        save_data_array(fname, particle_count, base_count, particle_count - base_count)
    return


def plot_count_diff(radii, time, particle_count, fname=None):
    count_fig = plt.figure()
    count_ax = count_fig.add_subplot(111)

    particle_count = np.array(particle_count)

    l1, = count_ax.plot(time, particle_count, 'k-')

    count_ax.set_xlim(0., 100.)
    count_ax.set_xlabel('Time $(ps)$')
    count_ax.set_ylabel('Counts (arb. units)')
    count_ax.set_title('Particle Tracking Counts')

    if fname is not None:
        plt.savefig(fname + '.png', dpi=600, bbox_inches="tight", pad_inches=0.2)
        save_data_array(fname, particle_count)
    return


def plot_grid_heatmap(data, x1, x2, fname=None):
    near_fig = plt.figure()
    near_ax = near_fig.add_subplot(111)
    lin_fig = plt.figure()
    lin_ax = lin_fig.add_subplot(111)

    grid_x1, grid_x2 = np.mgrid[np.amin(x1):np.amax(x1):1000j, np.amin(x2):np.amax(x2):1000j]  # use imaginary numbers to specify how many steps, evenly spaced, to  use

    values = np.ones_like(x1)

    grid_z0 = griddata(data, values, (grid_x1, grid_x2), method='nearest')
    grid_z1 = griddata(data, values, (grid_x1, grid_x2), method='linear')

    shape = grid_z0.T.shape
    extent = [np.amin(x1), np.amax(x1), np.amin(x2), np.amax(x2)]

    near_ax.imshow(grid_z0.T, extent=extent, aspect=pt.calculate_aspect(shape, extent), origin='lower')
    near_ax.plot(x1, x2, 'k.')
    near_ax.set_title('Nearest')

    lin_ax.imshow(grid_z1.T, extent=extent, aspect=pt.calculate_aspect(shape, extent), origin='lower')
    lin_ax.plot(x1, x2, 'k.')
    lin_ax.set_title('Linear')

    if fname is not None:
        near_fig.savefig('near_' + fname + '.png', dpi=600, bbox_inches="tight")
        lin_fig.savefig('lin_' + fname + '.png', dpi=600, bbox_inches="tight")
    return


def show_plots():
    plt.show()
    return
