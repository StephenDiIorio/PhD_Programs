#! /usr/bin/env python

# import os
# import glob
try:
    import matplotlib as mpl
    mpl.use('Agg')
    import sdf
    import sdf_helper as sh
    import scipy.constants as sc
    import numpy as np
    # import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as colors
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
        except ImportError as e:
            print(e)
            pass
    # import matplotlib.pyplot as plt
    # import matplotlib.animation as animation
except ImportError as e:
    print(e)
    exit()


def get_varname(varname, species=None):
    if species is None:
        return varname
    else:
        return varname + "_" + species


def set_data_constants(fname, species=None):
    global units, time, step
    sdfdata = sdf.read(fname)

    units = sdfdata.__dict__["Electric_Field_Ey"].units
    time = sdfdata.Header['time']
    step = sdfdata.Header['step']
    return


def get_scale(grid, data, scale=False):
    global verbose

    if verbose > 1:
        print('Scaling not implemented with this yet.')
        # print('Calculating scaling factor')

    r2 = 1
    # if scale:
    #     data = sdf.read(filename, mmap=0)
    #     var = data.__dict__[varname]
    #     r0 = data.Grid_Grid_Orig.data[0][0]**2
    #     if (var.grid.dims[0] == 1):
    #         x = var.grid.data[0][0, :, :]
    #         y = var.grid.data[1][0, :, :]
    #     elif (var.grid.dims[1] == 1):
    #         x = var.grid.data[0][:, 0, :]
    #         y = var.grid.data[1][:, 0, :]
    #     else:
    #         x = var.grid.data[0][:, :, 0]
    #         y = var.grid.data[1][:, :, 0]
    #         z = var.grid.data[2][:, :, 0]
    #         x = np.sqrt(x * x + y * y)
    #         y = z
    #     r2 = x**2 + y**2
    #     r2 = 0.5 * (r2[1:, :] + r2[:-1, :])
    #     r2 = 0.5 * (r2[:, 1:] + r2[:, :-1])
    #     r2 = r2 / r0

    return r2


def get_var_range(grid, data, scale):
    """Get a the data range for a given variable across an entire run"""

    vmin = float("inf")
    vmax = -float("inf")

    vscale = get_scale(grid, data, scale=scale)

    var = vscale * data
    var_min = var.min()
    var_max = var.max()
    if var_min < vmin:
        vmin = var_min
    if var_max > vmax:
        vmax = var_max

    return vmin, vmax


def get_si_prefix(scale, full_units=False):
    scale = abs(scale)
    mult = 1
    sym = ''

    if scale < 1e-24:
        full_units = True
    elif scale < 1e-21:
        # yocto
        mult = 1e24
        sym = 'y'
    elif scale < 1e-19:
        # zepto
        mult = 1e21
        sym = 'z'
    elif scale < 1e-16:
        # atto
        mult = 1e18
        sym = 'a'
    elif scale < 1e-13:
        # femto
        mult = 1e15
        sym = 'f'
    elif scale < 1e-10:
        # pico
        mult = 1e12
        sym = 'p'
    elif scale < 1e-7:
        # nano
        mult = 1e9
        sym = 'n'
    elif scale < 1e-4:
        # micro
        mult = 1e6
        sym = '{\mu}'
    elif scale < 1e-1:
        # milli
        mult = 1e3
        sym = 'm'
    elif scale >= 1e27:
        full_units = True
    elif scale >= 1e24:
        # yotta
        mult = 1e-24
        sym = 'Y'
    elif scale >= 1e21:
        # zetta
        mult = 1e-21
        sym = 'Z'
    elif scale >= 1e18:
        # exa
        mult = 1e-18
        sym = 'E'
    elif scale >= 1e15:
        # peta
        mult = 1e-15
        sym = 'P'
    elif scale >= 1e12:
        # tera
        mult = 1e-12
        sym = 'T'
    elif scale >= 1e9:
        # giga
        mult = 1e-9
        sym = 'G'
    elif scale >= 1e6:
        # mega
        mult = 1e-6
        sym = 'M'
    elif scale >= 1e3:
        # kilo
        mult = 1e-3
        sym = 'k'

    if full_units:
        scale = scale * mult
        if scale <= 0:
            pwr = 0
        else:
            pwr = (-np.floor(np.log10(scale)))
        mult = mult * np.power(10.0, pwr)
        if np.rint(pwr) != 0:
            sym = "(10^{%.0f})" % (-pwr) + sym

    return mult, sym


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


def plot_figure(grid, vdata, vmin=None, vmax=None, xscale=0, yscale=0, scale=False, figsize=None, ax=None):
    """Plot the given variable for each file from a simulation"""

    from matplotlib.ticker import FuncFormatter
    global verbose, title_info, time, step

    if verbose > 1:
        # print('Plotting {} from file {}'.format(varname, filename))
        print('Plotting Electric Field from temp and num dens')
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

    vscale = get_scale(grid, vdata, scale=scale)

    i0 = 1
    i1 = 0
    if len(grid.dims) == 2:
        x = grid.data[i0]
        y = grid.data[i1]
    else:
        print('CRITICAL ERROR: Not implemented yet')
        # if (grid.dims[0] == 1):
        #     orig = var.grid.name.replace('/', '_') + '_Orig'
        #     grid = data.__dict__[orig]
        #     i0 = 2
        #     i1 = 1
        #     x = grid.data[i0]
        #     y = grid.data[i1]
        # elif (grid.dims[1] == 1):
        #     i0 = 0
        #     i1 = 1
        #     x = grid.data[0][:, 0, :]
        #     y = grid.data[1][:, 0, :]
        # else:
        #     i0 = 0
        #     i1 = 2
        #     x = grid.data[0][:, :, 0]
        #     y = grid.data[1][:, :, 0]
        #     z = grid.data[2][:, :, 0]
        #     x = np.sqrt(x * x + y * y)
        #     y = z

    if verbose > 1:
        print('Scaling data')

    vdata = vdata * vscale

    if vmin is None:
        vmin, vmax = get_var_range(grid, vdata, scale=scale)

    if verbose > 1:
        print('Plotting data')
    print x
    print y
    print vdata
    print len(x)
    print len(y)
    print len(vdata)
    # im = ax.pcolormesh(x, y, vdata, vmin=vmin, vmax=vmax)
    # im = ax.pcolormesh(x, y, vdata, norm=colors.Normalize(vmin=vmin, vmax=vmax))
    # im = ax.pcolormesh(x, y, vdata, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    im = ax.pcolormesh(x, y, vdata, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=vmin, vmax=vmax))

    if verbose > 1:
        print('Plotting axes')

    def fmt(x, pos):
        return r'${}$'.format(x)

    # xmult, xsym = get_si_prefix(np.max(x) - np.min(x))
    # ymult, ysym = get_si_prefix(np.max(y) - np.min(y))

    if xscale == 0:
        xmult, xsym = get_si_prefix(np.max(x) - np.min(x))
    else:
        xmult, xsym = get_si_prefix(xscale)

    if yscale == 0:
        ymult, ysym = get_si_prefix(np.max(y) - np.min(y))
    else:
        ymult, ysym = get_si_prefix(yscale)

    if verbose > 1:
        print('Scale axis by {} ({}, {})'.format(xmult, np.min(x), np.max(x)))

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * xmult)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (y * ymult)))

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
    # if verbose > 0 and plot_figure.first:
    #     plot_figure.first = False
    #     print('Scale colorbar by {} ({}, {})'.format(mult, vmin, vmax))

    # data_label = var.name + ' $(' + sym + var.units + ')$'
    data_label = 'Electric Field' + ' $(' + sym + 'V/m)$'
    if title_info:
        legend_label = ''
    else:
        legend_label = data_label

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="15%")
    cbar = fig.colorbar(im, cax=cax, ax=ax, label=legend_label,
                        format=FuncFormatter(lambda x, y: x * mult))

    if title_info:
        title_label = data_label + ', '
    else:
        title_label = ''

    fig.fs = 'large'
    fig.suptitle(title_label + 'step={}, time={}'.format(step,
                 time), fontsize=fig.fs, y=0.98)
    fig.sca(ax)
    fig.canvas.draw()

    if verbose > 1:
        print('Done plotting Electric Field from file')

    return im, ax


def plot_first_figure(grid, data, vvals=None, xscale=0, yscale=0, scale=False):
    global verbose, title_info, units, step, time

    if verbose > 1:
        print('Getting data range')
    if vvals is None:
        vmin, vmax = get_var_range(grid, data, scale=scale)
    else:
        vmin = vvals[0]
        vmax = vvals[1]
    if verbose > 1:
        print('Found data range ({}, {})'.format(vmin, vmax))

    im, ax = plot_figure(grid, data, vmin=vmin, vmax=vmax, xscale=xscale, yscale=yscale, scale=scale)
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

    im, ax = plot_figure(grid, data, vmin=vmin, vmax=vmax, scale=scale, figsize=figsize)
    fig = ax.get_figure()

    # Get positions of the title's step and time text fields so that they
    # can be updated when animating

    # data = sdf.read(file_list[-1], mmap=0)
    # var = data.__dict__[varname]
    mult, sym = get_si_prefix(vmax - vmin)

    data_label = 'Electric Field' + ' $(' + sym + units + ')$, '

    if title_info:
        title_label = data_label + ', '
        title_label2 = data_label + ',_'
    else:
        title_label = ''
        title_label2 = ''

    fig._suptitle.set_text(title_label +
                           'step={}, time={}'.format(step, time))
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

    tt = fig.text(fig.x1, y, 'step={},_'.format(step),
                  fontsize=fig.fs, ha='left', va='bottom')
    bbox = tt.get_window_extent(renderer)
    fig.x2 = inv.transform(bbox)[1][0]
    tt.set_visible(False)

    fig._suptitle.set_visible(False)

    return im, fig


def main():
    global verbose, title_info
    verbose = 2
    title_info = False

    # species = "Hydrogen"
    species = None
    distribution = "Electric_Field_Ey"
    fname = "0010.sdf"
    set_data_constants(fname, species=species)

    dist_key = get_varname(distribution, species=species)
    print dist_key
    # sdfdata = sdf.read(fname)
    sdfdata = sh.getdata(fname)
    grid = sdfdata.__dict__[dist_key].grid
    data = sdfdata.__dict__[dist_key].data
    print data
    # print sdfdata.__dict__['dist_fn_x_en_Hydrogen'].data

    vmin = 1e-10
    vmax = 1e10
    xscale = 1e-6
    yscale = 0
    scale = False
    im, fig = plot_first_figure(grid, data, vvals=[vmin, vmax], scale=scale)
    fig.savefig('out.png', bbox_inches='tight', dpi=200)
    return


if __name__ == "__main__":
    main()
