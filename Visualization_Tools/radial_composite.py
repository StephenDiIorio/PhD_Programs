import sdf
import glob
import numpy as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FuncFormatter
# plt.ion()

from scipy.ndimage import map_coordinates
from scipy.ndimage.interpolation import shift
from scipy.optimize import curve_fit, minimize

plt.rc('font', size=20)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=15)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize


def reproject_image_into_polar(data, origin=None, Jacobian=False, dr=1, dt=None):
    """
    Reprojects a 2D numpy array (``data``) into a polar coordinate system.
    "origin" is a tuple of (x0, y0) relative to the bottom-left image corner,
    and defaults to the center of the image.
    Parameters
    ----------
    data : 2D np.array
    origin : tuple
        The coordinate of the image center, relative to bottom-left
    Jacobian : boolean
        Include ``r`` intensity scaling in the coordinate transform.
        This should be included to account for the changing pixel size that
        occurs during the transform.
    dr : float
        Radial coordinate spacing for the grid interpolation
        tests show that there is not much point in going below 0.5
    dt : float
        Angular coordinate spacing (in radians)
        if ``dt=None``, dt will be set such that the number of theta values
        is equal to the maximum value between the height or the width of
        the image.
    Returns
    -------
    output : 2D np.array
        The polar image (r, theta)
    r_grid : 2D np.array
        meshgrid of radial coordinates
    theta_grid : 2D np.array
        meshgrid of theta coordinates

    Notes
    -----
    Adapted from:
    http://stackoverflow.com/questions/3798333/image-information-along-a-polar-coordinate-system
    """

    # bottom-left coordinate system requires numpy image to be np.flipud
    data = np.flipud(data)

    ny, nx = data.shape[:2]
    if origin is None:
        origin = (nx//2, ny//2)

    # Determine that the min and max r and theta coords will be...
    x, y = index_coords(data, origin=origin)  # (x,y) coordinates of each pixel
    r, theta = cart2polar(x, y)  # convert (x,y) -> (r,θ), note θ=0 is vertical

    nr = np.int(np.ceil((r.max()-r.min())/dr))

    if dt is None:
        nt = max(nx, ny)
    else:
        # dt in radians
        nt = np.int(np.ceil((theta.max()-theta.min())/dt))

    # Make a regular (in polar space) grid based on the min and max r & theta
    r_i = np.linspace(r.min(), r.max(), nr, endpoint=False)
    theta_i = np.linspace(theta.min(), theta.max(), nt, endpoint=False)
    theta_grid, r_grid = np.meshgrid(theta_i, r_i)

    # Project the r and theta grid back into pixel coordinates
    X, Y = polar2cart(r_grid, theta_grid)

    X += origin[0]  # We need to shift the origin
    Y += origin[1]  # back to the bottom-left corner...
    xi, yi = X.flatten(), Y.flatten()
    coords = np.vstack((yi, xi))  # (map_coordinates requires a 2xn array)

    zi = map_coordinates(data, coords)
    output = zi.reshape((nr, nt))

    if Jacobian:
        output = output*r_i[:, np.newaxis]

    return output, r_grid, theta_grid


def index_coords(data, origin=None):
    """
    Creates x & y coords for the indicies in a numpy array

    Parameters
    ----------
    data : numpy array
        2D data
    origin : (x,y) tuple
        defaults to the center of the image. Specify origin=(0,0)
        to set the origin to the *bottom-left* corner of the image.

    Returns
    -------
        x, y : arrays
    """

    ny, nx = data.shape[:2]
    if origin is None:
        origin_x, origin_y = nx//2, ny//2
    else:
        origin_x, origin_y = origin

    x, y = np.meshgrid(np.arange(float(nx)), np.arange(float(ny)))

    x -= origin_x
    y -= origin_y
    return x, y


def cart2polar(x, y):
    """
    Transform Cartesian coordinates to polar

    Parameters
    ----------
    x, y : floats or arrays
        Cartesian coordinates

    Returns
    -------
    r, theta : floats or arrays
        Polar coordinates

    """

    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(x, y)  # θ referenced to vertical
    return r, theta


def polar2cart(r, theta):
    """
    Transform polar coordinates to Cartesian

    Parameters
    -------
    r, theta : floats or arrays
        Polar coordinates

    Returns
    ----------
    x, y : floats or arrays
        Cartesian coordinates
    """

    y = r * np.cos(theta)   # θ referenced to vertical
    x = r * np.sin(theta)
    return x, y


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


def get_var_range(file_list, varname):
    """Get a the data range for a given variable across an entire run"""

    vmin = float("inf")
    vmax = -float("inf")

    for f in file_list:
        try:
            data = sdf.read(f, mmap=0)
            var = data.__dict__[varname].data
            var_min = var.min()
            var_max = var.max()
            if var_min < vmin:
                vmin = var_min
            if var_max > vmax:
                vmax = var_max
        except:
            pass

    if verbose > 0:
        print(vmin, vmax)

    return vmin, vmax


def get_files(wkdir='Data'):
    """Get a list of SDF filenames belonging to the same run"""
    flist = glob.glob(wkdir + "/[0-9]*.sdf")
    flist = sorted(flist)

    return flist


def clean_file_list(file_list, varname):
    new_file_list = []

    for f in file_list:
        try:
            data = sdf.read(f, mmap=0)
            data.__dict__[varname]
            new_file_list.append(f)
        except KeyError:
            pass
    return new_file_list


def calculate_aspect(shape, extent):
    dx = (extent[1] - extent[0]) / float(shape[1])
    dy = (extent[3] - extent[2]) / float(shape[0])
    if dx / dy > 0:
        return dx / dy
    else:
        return -dx / dy


def electric_radial_average(sdf_data):
    ex = sdf_data.__dict__['Electric_Field_Ex']
    ey = sdf_data.__dict__['Electric_Field_Ey']

    x_grid = ex.grid_mid.data[0]
    y_grid = ey.grid_mid.data[1]
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


def composite_field_plot(varname, vmin=None, vmax=None, directory='Data'):
    global verbose, dpi

    file_list = get_files(wkdir=directory)
    file_list = clean_file_list(file_list, varname)

    file_list.remove(directory + '0000.sdf')
    file_list.remove(directory + '0002.sdf')
    file_list.remove(directory + '0003.sdf')
    file_list.remove(directory + '0004.sdf')
    file_list.remove(directory + '0005.sdf')
    file_list.remove(directory + '0006.sdf')
    file_list.remove(directory + '0007.sdf')
    file_list.remove(directory + '0008.sdf')
    file_list.remove(directory + '0009.sdf')
    file_list.remove(directory + '0011.sdf')
    file_list.remove(directory + '0012.sdf')
    file_list.remove(directory + '0013.sdf')
    file_list.remove(directory + '0014.sdf')
    file_list.remove(directory + '0015.sdf')
    file_list.remove(directory + '0016.sdf')
    file_list.remove(directory + '0017.sdf')
    file_list.remove(directory + '0018.sdf')
    file_list.remove(directory + '0019.sdf')
    file_list.remove(directory + '0021.sdf')
    file_list.remove(directory + '0022.sdf')
    file_list.remove(directory + '0023.sdf')
    file_list.remove(directory + '0024.sdf')
    file_list.remove(directory + '0025.sdf')
    file_list.remove(directory + '0026.sdf')
    file_list.remove(directory + '0027.sdf')
    file_list.remove(directory + '0028.sdf')
    file_list.remove(directory + '0029.sdf')
    file_list.remove(directory + '0031.sdf')
    file_list.remove(directory + '0032.sdf')
    file_list.remove(directory + '0033.sdf')
    file_list.remove(directory + '0034.sdf')
    file_list.remove(directory + '0035.sdf')
    file_list.remove(directory + '0036.sdf')
    file_list.remove(directory + '0037.sdf')
    file_list.remove(directory + '0038.sdf')
    file_list.remove(directory + '0039.sdf')
    file_list.remove(directory + '0041.sdf')
    file_list.remove(directory + '0042.sdf')
    file_list.remove(directory + '0043.sdf')
    file_list.remove(directory + '0044.sdf')
    file_list.remove(directory + '0045.sdf')
    file_list.remove(directory + '0046.sdf')
    file_list.remove(directory + '0047.sdf')
    file_list.remove(directory + '0048.sdf')
    file_list.remove(directory + '0049.sdf')
    file_list.remove(directory + '0051.sdf')
    file_list.remove(directory + '0052.sdf')
    file_list.remove(directory + '0053.sdf')
    file_list.remove(directory + '0054.sdf')
    file_list.remove(directory + '0055.sdf')
    file_list.remove(directory + '0056.sdf')
    file_list.remove(directory + '0057.sdf')
    file_list.remove(directory + '0058.sdf')
    file_list.remove(directory + '0059.sdf')
    file_list.remove(directory + '0061.sdf')
    file_list.remove(directory + '0062.sdf')
    file_list.remove(directory + '0063.sdf')
    file_list.remove(directory + '0064.sdf')
    file_list.remove(directory + '0065.sdf')
    file_list.remove(directory + '0066.sdf')
    file_list.remove(directory + '0067.sdf')
    file_list.remove(directory + '0068.sdf')
    file_list.remove(directory + '0069.sdf')
    file_list.remove(directory + '0071.sdf')
    file_list.remove(directory + '0072.sdf')
    file_list.remove(directory + '0073.sdf')
    file_list.remove(directory + '0074.sdf')
    file_list.remove(directory + '0075.sdf')
    file_list.remove(directory + '0076.sdf')
    file_list.remove(directory + '0077.sdf')
    file_list.remove(directory + '0078.sdf')
    file_list.remove(directory + '0079.sdf')
    file_list.remove(directory + '0081.sdf')
    file_list.remove(directory + '0082.sdf')
    file_list.remove(directory + '0083.sdf')
    file_list.remove(directory + '0084.sdf')
    file_list.remove(directory + '0085.sdf')
    file_list.remove(directory + '0086.sdf')
    file_list.remove(directory + '0087.sdf')
    file_list.remove(directory + '0088.sdf')
    file_list.remove(directory + '0089.sdf')
    file_list.remove(directory + '0091.sdf')
    file_list.remove(directory + '0092.sdf')
    file_list.remove(directory + '0093.sdf')
    file_list.remove(directory + '0094.sdf')
    file_list.remove(directory + '0095.sdf')
    file_list.remove(directory + '0096.sdf')
    file_list.remove(directory + '0097.sdf')
    file_list.remove(directory + '0098.sdf')
    file_list.remove(directory + '0099.sdf')

    if verbose > 0:
        print('Found {} files to plot'.format(len(file_list)))

    data = []
    for f in file_list:
        d = sdf.read(f, mmap=0)
        avg, r, t = electric_radial_average(d)
        data.append(avg)
    var = d.__dict__[varname]  # for units later
    data = np.asarray(data)
    data = data.T

    tmin = sdf.read(file_list[0], mmap=0).Header['time']
    tmax = sdf.read(file_list[-1], mmap=0).Header['time']
    rmin = np.min(r)
    rmax = np.max(r)

    shape = data.shape
    extent = [tmin, tmax, rmax, rmin]

    rmult, rsym = get_si_prefix(rmax - rmin)  # y axis
    tmult, tsym = get_si_prefix(tmax - tmin)  # x axis

    if vmin is None and vmax is None:
        vmin = np.min(data)
        vmax = np.max(data)
        # vmin, vmax = get_var_range(file_list, varname)
    elif vmin is None:
        vmin = np.min(data)
        # vmin = get_var_range(file_list, varname)[0]
    elif vmax is None:
        vmax = np.max(data)
        # vmax = get_var_range(file_list, varname)[1]
    mult, sym = get_si_prefix(vmax - vmin)

    fig, ax = plt.subplots()
    im = ax.imshow(data, extent=extent, aspect=calculate_aspect(shape, extent), interpolation='none', cmap=cm.coolwarm, vmin=vmin, vmax=vmax, origin='lower')

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * tmult)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * rmult)))
    plt.xlabel('t $(' + tsym + 's)$')
    plt.ylabel('r $(' + rsym + 'm)$')
    # data_label = var.name + ' $(' + sym + var.units + ')$'
    data_label = 'Radial Electric Field $(' + sym + var.units + ')$'
    plt.title('Radial Electric Field Evolution')

    cbar = fig.colorbar(im, label=data_label, format=FuncFormatter(lambda x, y: x * mult))
    plt.tight_layout()
    plt.savefig('rad_test.png', dpi=600, bbox_inches="tight")
    # plt.show()


if __name__ == "__main__":
    import argparse
    global verbose, dpi

    verbose = 0
    dpi = 300

    varname = 'Electric_Field_Ex'
    vmin = None
    vmax = None

    parser = argparse.ArgumentParser(description='''
    This composite plot of variable evolution
    ''')
    # parser.add_argument('variable', type=str, default=varname, nargs='?',
    #                     help="Variable to plot")
    # EDIT: TAKE IN MULTIPLE VARIABLES TO PLOT
    parser.add_argument('variable', type=str, default=varname, nargs='*',
                        help="Variable(s) to plot")
    parser.add_argument('-v', '--verbose', action="count", default=0,
                        help="Increase verbosity")
    # EDIT: OPTION TO SEND THE DIRECTORY NAME WITH DATA FILES
    parser.add_argument('-d', '--directory', type=str, default='Data', nargs=1,
                        help="Directory for one of the simulation dumps")
    parser.add_argument('-n', '--noneg', action="count", default=0,
                        help="Makes vmin positive")
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

    # EDIT: CALL PLOT_FIGURES WITH DIRECTORY ARGUMENT
    # take first element of directory because parser gives args in a list
    composite_field_plot(var, vmin=vmin, vmax=vmax, directory=args.directory[0])
