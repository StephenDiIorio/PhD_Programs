import sys
import numpy as np
try:
    import sdf
    sdf_mod_name = 'sdf'
except ImportError:
    pass


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


def get_var_range(data):
    """Get a the data range for a given variable across an entire run."""

    vmin = float("inf")
    vmax = -float("inf")

    var_min = data.min()
    var_max = data.max()
    if var_min < vmin:
        vmin = var_min
    if var_max > vmax:
        vmax = var_max

    return vmin, vmax


def get_var_range_from_sdf_files(file_list, varname):
    """Get a the data range for a given variable across an entire run"""
    if sdf_mod_name not in sys.modules:
        raise ImportError("SDF module from EPOCH is not properly installed. Cannot use this function.")

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

    return vmin, vmax


def calculate_aspect(shape, extent):
    dx = (extent[1] - extent[0]) / float(shape[1])
    dy = (extent[3] - extent[2]) / float(shape[0])
    if dx / dy > 0:
        return dx / dy
    else:
        return -dx / dy
