import numpy as np

import sdf


def first_order_weight_2d(x, y, dx, dy, p_values, values=None, weight=None):
    """
    Perform first-order (linear or area) weighting for mapping
    particles onto a grid

    Parameters
    ----------
    x, y : arrays containing the physical x and y values for the
        original grid
    dx, dy : floats providing the physical spacing between cells
    p_values : array containing the physical particle positions;
        elements should be of form (x_pos, y_pos)
    values (opt) : an array containing some physical quantity
        of the particle (velocity, temp, etc.) to weight onto
        the grid
    weight (opt) : an array containing the particle weightings
        if they have unequal weighting

    Returns
    -------
    weighted_grid : the grid with the weighted particle
        values mapped onto it; has the same shape as
        (size(x) x size(y)); the origin of this matrix
        is the upper left corner
    """

    weighted_grid = np.zeros((np.size(x), np.size(y)), dtype=np.float64)

    # x_grid, y_grid = np.meshgrid(x, y)

    x_min = np.min(x)
    y_min = np.min(y)

    for p_index, p in enumerate(p_values):
        x_pos = p[0]
        y_pos = p[1]

        fi = (x_pos - x_min) / dx
        i = np.floor(fi, dtype=np.int64)
        hx = fi - i

        fj = (y_pos - y_min) / dy
        j = np.floor(fj, dtype=np.int64)
        hy = fj - j

        if values is not None:
            v = values[p_index]
        else:
            v = 1.0
        if weight is not None:
            w = weight[p_index]
        else:
            w = 1.0

        weighted_grid[i, j]     += (1 - hx) * (1 - hy) * v * w
        weighted_grid[i+1, j]   += (1 - hx) * hy       * v * w
        weighted_grid[i, j+1]   += hx       * (1 - hy) * v * w
        weighted_grid[i+1, j+1] += hx       * hy       * v * w

    return weighted_grid
