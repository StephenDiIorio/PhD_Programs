# TODO: specifying only vmin or vmax might leave one larger (or smaller)
# than the other (vmin might be greater than vmax for example)
# TODO: allow negitve number be passed in as command line argument for
# vmin and vmax
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
from matplotlib.ticker import FuncFormatter
import numpy as np

import sdf
from PlottingTools import get_si_prefix, get_var_range_from_sdf_files

try:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except:
    try:
        # Workaround for broken macOS installation
        import os
        import sys
        import matplotlib
        sys.path.append(os.path.join(matplotlib.__path__[0],
                                     '..', 'mpl_toolkits'))
        from axes_grid1 import make_axes_locatable
    except:
        pass


dpi = 300
verbose = 0


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
    #         # file_job_id = data.Header['jobid1']
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


def find_renderer(fig):
    if hasattr(fig.canvas, "get_renderer"):
        # Some backends, such as TkAgg, have the get_renderer method, which
        # makes this easy.
        renderer = fig.canvas.get_renderer()
    else:
        # Other backends do not have the get_renderer method, so we have a work
        # around to find the renderer.  Print the figure to a temporary file
        # object, and then grab the renderer that was used.
        # (I stole this trick from the matplotlib backend_bases.py
        # print_figure() method.)
        import io
        fig.canvas.print_png(io.BytesIO())
        renderer = fig._cachedRenderer
    return(renderer)


def plot_figure(filename, varname, vmin=None, vmax=None, figsize=None, ax=None):
    """
    Plot the given variable for each file from a simulation
    """

    global verbose, title_info

    if verbose > 1:
        print('Plotting {} from file {}'.format(varname, filename))
    axis_type = 1  # 0 - off, 1 - left+bottom, 2 - all

    if ax is not None:
        fig = ax.get_figure()
    elif figsize is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    else:
        fig, ax = plt.subplots(1, 1, figsize=figsize)

    if verbose > 1:
        print('Adjusting subplot')

    if axis_type == 0:
        plt.subplots_adjust(left=0, right=0.9, bottom=0.05, top=0.92,
                            wspace=0, hspace=0)
    else:
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.92,
                            wspace=0, hspace=0)
    plt.set_cmap(cm.coolwarm)

    if verbose > 1:
        print('Reading data')

    # Draw initial plot
    f = filename
    data = sdf.read(f)
    var = data.__dict__[varname]
    grid = var.grid
    vdata = var.data

    i0 = 1
    i1 = 0
    if len(grid.dims) == 2:
        x = grid.data[i0]
        y = grid.data[i1]
    else:
        if (grid.dims[0] == 1):
            orig = var.grid.name.replace('/', '_') + '_Orig'
            grid = data.__dict__[orig]
            i0 = 2
            i1 = 1
            x = grid.data[i0]
            y = grid.data[i1]
        elif (grid.dims[1] == 1):
            i0 = 0
            i1 = 1
            x = grid.data[0][:, 0, :]
            y = grid.data[1][:, 0, :]
        else:
            i0 = 0
            i1 = 2
            x = grid.data[0][:, :, 0]
            y = grid.data[1][:, :, 0]
            z = grid.data[2][:, :, 0]
            x = np.sqrt(x * x + y * y)
            y = z

    if verbose > 1:
        print('Scaling data')

    vdata = vdata

    if vmin is None:
        vmin, vmax = get_var_range_from_sdf_files((filename,), varname)

    if verbose > 1:
        print('Plotting data')
    im = ax.pcolormesh(x, y, vdata, vmin=vmin, vmax=vmax)

    if verbose > 1:
        print('Plotting axes')

    xmult, xsym = get_si_prefix(np.max(x) - np.min(x))
    ymult, ysym = get_si_prefix(np.max(y) - np.min(y))
    if verbose > 1:
        print('Scale axis by {} ({}, {})'.format(xmult, np.min(x), np.max(x)))

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: '{:.2f}'.format(x * xmult)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: '{:.2f}'.format(x * ymult)))

    ax.set_xlabel(grid.labels[i0] + ' $(' + xsym + grid.units[i0] + ')$')
    ax.set_ylabel(grid.labels[i1] + ' $(' + ysym + grid.units[i1] + ')$')
    ax.axis('image')

    if axis_type == 0:
        ax.axis('off')
    elif axis_type == 1:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

    # Add colorbar
    mult, sym = get_si_prefix(vmax - vmin)
    if verbose > 0 and plot_figure.first:
        plot_figure.first = False
        print('Scale colorbar by {} ({}, {})'.format(mult, vmin, vmax))

    data_label = var.name + ' $(' + sym + var.units + ')$'
    if title_info:
        legend_label = ''
    else:
        legend_label = data_label

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="5%")
    cbar = fig.colorbar(im, cax=cax, ax=ax, label=legend_label,
                        format=FuncFormatter(lambda x, y: '{:.2f}'.format(x * mult)))

    if title_info:
        title_label = data_label + ', '
    else:
        title_label = ''

    fig.fs = 'large'
    fig.suptitle(title_label + 'step={}, time={}'.format(data.Header['step'],
                 data.Header['time']), fontsize=fig.fs, y=0.98)
    fig.sca(ax)
    fig.canvas.draw()

    if verbose > 1:
        print('Done plotting {} from file {}'.format(varname, filename))

    return im, ax


plot_figure.first = True


def plot_first_figure(file_list, varname, vmin=None, vmax=None):
    global verbose, title_info

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
    im, ax = plot_figure(filename, varname, vmin, vmax)
    fig = ax.get_figure()
    renderer = find_renderer(fig)
    box_orig = fig.bbox_inches
    box_raw = fig.get_tightbbox(renderer)
    pp = box_raw.get_points()
    figsize = (pp[1, 0] - pp[0, 0], pp[1, 1] - pp[0, 1])

    ratio = figsize[0] / figsize[1]
    if ratio > 1:
        w = box_orig.width
        h = w / ratio
        h = round(10 * h) / 10.
        figsize = (w, h)
    else:
        h = box_orig.height
        w = h * ratio
        w = round(10 * w) / 10.
        figsize = (w, h)

    im, ax = plot_figure(filename, varname, vmin, vmax, figsize)
    fig = ax.get_figure()

    # Get positions of the title's step and time text fields so that they
    # can be updated when animating

    data = sdf.read(file_list[-1])
    var = data.__dict__[varname]
    mult, sym = get_si_prefix(vmax - vmin)

    data_label = var.name + ' $(' + sym + var.units + ')$, '

    if title_info:
        title_label = data_label + ', '
        title_label2 = data_label + ',_'
    else:
        title_label = ''
        title_label2 = ''

    fig._suptitle.set_text(title_label +
                           'step={}, time={}'.format(data.Header['step'],
                                                     data.Header['time']))
    bbox = fig._suptitle.get_window_extent(renderer)
    inv = fig.transFigure.inverted()

    x0 = inv.transform(bbox)[0][0]
    y = inv.transform(bbox)[0][1]
    fig.y = y

    tt = fig.text(x0, y, title_label2,
                  fontsize=fig.fs, ha='left', va='bottom')
    bbox = tt.get_window_extent(renderer)
    fig.x1 = inv.transform(bbox)[1][0]
    fig.y = inv.transform(bbox)[0][1]
    tt.set_text(title_label)

    tt = fig.text(fig.x1, y, 'step={},_'.format(data.Header['step']),
                  fontsize=fig.fs, ha='left', va='bottom')
    bbox = tt.get_window_extent(renderer)
    fig.x2 = inv.transform(bbox)[1][0]
    tt.set_visible(False)

    fig._suptitle.set_visible(False)

    return im, fig


# EDIT: ADD DIRECTORY AS OPTIONAL ARGUMENT
def plot_figures(varname, vmin=None, vmax=None, directory='Data', base=None):
    """
    Plot the given variable for each file from a simulation
    """

    global verbose

    # file_list = get_files(base=base)
    # EDIT: CALL GET_FILES WITH DIRECTORY ARGUMENT
    file_list = get_files(wkdir=directory, base=base)
    file_list = clean_file_list(file_list, varname)
    if verbose > 0:
        print('Found {} files to plot'.format(len(file_list)))
    im, fig = plot_first_figure(file_list, varname, vmin, vmax)
    fig.savefig(varname + '.png', bbox_inches='tight', dpi=200)

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

        grid = var.grid
        i1 = 0
        y = grid.data[i1]

        ax = plt.gca()

        im.set_array(vdata.ravel())

        # Hack to change y axis for moving window. Change tick labels
        # instead of changing the limits of axes. This plays better
        # with pcolormesh.
        ymult, ysym = get_si_prefix(np.max(y) - np.min(y))
        curr_ticks = ax.get_yticks()
        new_ticks = np.linspace(np.min(y), np.max(y), num=np.size(curr_ticks))
        new_ticks = np.array(['{:.2f}'.format(x * ymult) for x in new_ticks])
        ax.set_yticklabels(new_ticks)
        ax.set_ylabel(grid.labels[i1] + ' $(' + ysym + grid.units[i1] + ')$')
        plt.title('step={}, time={}'.format(data.Header['step'],
                                            data.Header['time']))

        return im,

    ani = animation.FuncAnimation(fig, update, init_func=init,
                                  frames=file_list[::1],
                                  blit=False, interval=1)
    fps = len(file_list) / 8.
    writer = animation.writers['ffmpeg'](fps=fps, bitrate=5e5, codec='libx264',
                                         extra_args=['-pix_fmt', 'yuv420p'])
    ani.save(varname + '.mp4', writer=writer, dpi=dpi)
    return ani


def list_variables(data):
    import numpy as np
    dct = data.__dict__
    for key in sorted(dct):
        try:
            val = dct[key]
            if type(val) == sdf.BlockPlainVariable and len(val.dims) == 2 \
                    and val.dims[0] > 1 and val.dims[1] > 1:
                print('{} {}'.format(key,
                      np.array2string(np.array(val.dims), separator=', ')))
        except:
            pass


def list_file_variables(base=None):
    file_list = get_files(base=base)
    f = file_list[0]
    data = sdf.read(f)
    list_variables(data)


def main():
    import argparse
    global verbose, title_info

    varname = 'Electric_Field_Ey'
    vmin = None
    vmax = None
    # EDIT: CHANGE DEFAULT TO A LIST OF VARIABLES TO PLOT
    # varname = ['Electric_Field_Ey']

    parser = argparse.ArgumentParser(description='''
    This generates a movie from a series of SDF files
    ''')
    # parser.add_argument('variable', type=str, default=varname, nargs='?',
    #                     help="Variable to plot")
    # EDIT: TAKE IN MULTIPLE VARIABLES TO PLOT
    parser.add_argument('variable', type=str, default=varname, nargs='*',
                        help="Variable to plot")
    parser.add_argument('-l', '--list', action="store_true",
                        help="List variables in SDF file")
    parser.add_argument('-t', '--title-info', action="store_true",
                        help="Plot variable info in the title")
    parser.add_argument('-v', '--verbose', action="count", default=0,
                        help="Increase verbosity")
    parser.add_argument('-f', '--filename', type=str, nargs=1,
                        help="Filename for one of the simulation dumps")
    # EDIT: OPTION TO SEND THE DIRECTORY NAME WITH DATA FILES
    parser.add_argument('-d', '--directory', type=str, default='Data', nargs=1,
                        help="Directory for one of the simulation dumps")
    parser.add_argument('-n', '--noneg', action="count", default=0,
                        help="Makes vmin positive")
    parser.add_argument('-i', '--vmin', type=float, default=vmin, nargs=1,
                        help="Vmin for data display (changes min of colorbar)")
    parser.add_argument('-a', '--vmax', type=float, default=vmax, nargs=1,
                        help="Vmax for data display (changes max of colorbar)")
    args = parser.parse_args()

    verbose = args.verbose
    title_info = args.title_info

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

    if args.list:
        list_file_variables(base=args.filename)
    else:
        # EDIT: CALL PLOT_FIGURES WITH DIRECTORY ARGUMENT
        # take first element of directory because parser gives args in a list
        plot_figures(var, vmin=vmin, vmax=vmax,
                     directory=args.directory[0], base=args.filename)


if __name__ == "__main__":
    main()
