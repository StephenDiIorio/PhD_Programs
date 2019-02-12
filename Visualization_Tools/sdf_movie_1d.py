from matplotlib import use
use('Agg')

import glob
import os.path
import sys

package_directory = os.path.dirname(os.path.abspath(__file__))  # Get path to current file
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))  # Trace path back to Utilities folder to import modules

import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.ticker import FuncFormatter

import sdf
from PlottingTools import get_si_prefix, get_var_range_from_sdf_files


plt.rc('font', size=20)        # controls default text sizes
plt.rc('axes', titlesize=18)   # fontsize of the axes title
plt.rc('axes', labelsize=15)   # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)  # fontsize of the tick labels
plt.rc('ytick', labelsize=12)  # fontsize of the tick labels
plt.rc('legend', fontsize=16)  # legend fontsize


def get_files(wkdir='Data', base=None):
    """
    Get a list of SDF filenames belonging to the same run
    """

    import os.path

    if base:
        # take first element of base because parser gives args as a list
        wkdir = os.path.dirname(base[0])
        flist = glob.glob(wkdir + "/*.sdf")
        flist.remove(base[0])
        flist = base + sorted(flist)
    else:
        flist = glob.glob(wkdir + "/[0-9]*.sdf")
        flist = sorted(flist)

    # Find the job id
    # for f in flist:
    #     try:
    #         data = sdf.read(f)
    #         job_id = data.Header['jobid1']
    #         break
    #     except:
    #         pass

    # file_list = []

    # # Add all files matching the job id
    # for f in sorted(flist):
    #     try:
    #         data = sdf.read(f)
    #         file_job_id = data.Header['jobid1']
    #         # if file_job_id == job_id:
    #         file_list.append(f)
    #     except:
    #         pass

    return flist  # file_list


def clean_file_list(file_list, varname):
    new_file_list = []

    for f in file_list:
        try:
            data = sdf.read(f)
            dummy = data.__dict__[varname]
            new_file_list.append(f)
        except KeyError:
            pass
    return new_file_list


def plot_figure(filename, varname, vmin=None, vmax=None):
    """
    Plot the given variable for each file from a simulation
    """

    global verbose

    if verbose > 1:
        print('Plotting {} from file {}'.format(varname, filename))

    fig = plt.figure()

    if verbose > 1:
        print('Reading data')
    print(filename)
    data = sdf.read(filename)
    var = data.__dict__[varname]
    vdata = var.data
    grid = var.grid_mid
    x = grid.data[0]

    if verbose > 1:
        print('Plotting data')
    im, = plt.plot(x, vdata, 'k')
    plt.title('step={}, time={}'.format(data.Header['step'],
                                        data.Header['time']))

    plt.axis('tight')

    # Get all labels for axes and data
    xmult, xsym = get_si_prefix(np.max(x) - np.min(x))  # x axis
    ymult, ysym = get_si_prefix(np.max(vdata) - np.min(vdata))  # y axis

    if verbose > 1:
        print('Plotting axes')
    if verbose > 1:
        print('Scale axis by {} ({}, {})'.format(xmult, np.min(x), np.max(x)))
        print('Scale axis by {} ({}, {})'.format(ymult, np.min(vdata), np.max(vdata)))

    ax = plt.gca()
    ticker.Locator.MAXTICKS = 2000
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5e-6))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5e-6))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1e9))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1e9))

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * xmult)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult)))

    plt.xlabel(grid.labels[0] + ' $(' + xsym + grid.units[0] + ')$')
    plt.ylabel(var.name + ' $(' + ysym + var.units + ')$')

    # remove spines on right and top of plot
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # only show ticks on left and bottom side for these axes
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    # point ticks outside of plot
    ax.yaxis.set_tick_params(direction='out')
    ax.xaxis.set_tick_params(direction='out')

    if verbose > 1:
        print('Done plotting {} from file {}'.format(varname, filename))

    return im, fig


def plot_first_figure(file_list, varname, vmin=None, vmax=None):
    global verbose

    if verbose > 1:
        print('Getting data range')
    if vmin is None and vmax is None:
        vmin, vmax = get_var_range_from_sdf_files(file_list, varname)
    elif vmin is None:
        vmin = get_var_range_from_sdf_files(file_list, varname)[0]
    elif vmax is None:
        vmax = get_var_range_from_sdf_files(file_list, varname)[1]

    if verbose > 1:
        print('Found data range ({}, {})'.format(vmin, vmax))

    filename = file_list[0]
    im, fig = plot_figure(filename, varname, vmin, vmax)

    return im, fig


# def plot_figures(varname, scale=False, base=None):
# EDIT: ADD DIRECTORY AS OPTIONAL ARGUMENT
def plot_figures(varname, vmin=None, vmax=None, directory='Data', base=None):
    """
    Plot the given variable for each file from a simulation
    """

    global verbose, dpi

    # file_list = get_files(base=base)
    # EDIT: CALL GET_FILES WITH DIRECTORY ARGUMENT
    file_list = get_files(wkdir=directory, base=base)
    file_list = clean_file_list(file_list, varname)

    if verbose > 0:
        print('Found {} files to plot'.format(len(file_list)))

    im, fig = plot_first_figure(file_list, varname, vmin, vmax)
    # fig.savefig(varname + '.png', bbox_inches='tight', dpi=200)

    # Draw plots
    def init():
        data = sdf.read(file_list[0])
        plt.title('step={}, time={}'.format(data.Header['step'],
                                            data.Header['time']))
        return im,

    # Draw plots
    def update(filename):
        if verbose > 0:
            print('Generating frame for file {}'.format(filename))
        data = sdf.read(filename)
        var = data.__dict__[varname]
        vdata = var.data

        im.set_ydata(vdata)

        ax = plt.gca()
        ticker.Locator.MAXTICKS = 2000
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5e-6))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5e-6))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1e9))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1e9))
        ax.set_ylim(np.min(vdata), np.max(vdata))
        ymult, ysym = get_si_prefix(np.max(vdata) - np.min(vdata))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * ymult)))
        plt.ylabel(var.name + ' $(' + ysym + var.units + ')$')

        plt.title('step={}, time={}'.format(data.Header['step'],
                                            data.Header['time']))
        return im,

    ani = animation.FuncAnimation(fig, update, init_func=init,
                                  frames=file_list[::1],
                                  blit=False, interval=1)
    fps = len(file_list) / 8.
    writer = animation.writers['ffmpeg'](fps=fps,
                                         bitrate=5e5,
                                         codec='libx264',
                                         extra_args=['-pix_fmt', 'yuv420p'])
    ani.save(varname + '.mp4', writer=writer, dpi=dpi)
    return ani


def main():
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
    parser.add_argument('-f', '--filename', type=str, nargs=1,
                        help="Filename for one of the simulation dumps")
    # EDIT: OPTION TO SEND THE DIRECTORY NAME WITH DATA FILES
    parser.add_argument('-d', '--directory', type=str, default='Data', nargs=1,
                        help="Directory for one of the simulation dumps")
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
        vmin = args.vmin[0]
    else:
        vmin = args.vmin
    if isinstance(args.vmax, list):
        vmax = args.vmax[0]
    else:
        vmax = args.vmax

    # plot_figures(args.variable, args.scale, base=args.filename)
    # EDIT: CALL PLOT_FIGURES WITH DIRECTORY ARGUMENT
    # take first element of directory because parser gives args in a list
    plot_figures(var, vmin=vmin, vmax=vmax, directory=args.directory[0],
                 base=args.filename)


if __name__ == "__main__":
    main()
