from matplotlib import use
use('Agg')

import glob

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, rc
from matplotlib.ticker import FuncFormatter

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


plt.rc('font', size=20)        # controls default text sizes
plt.rc('axes', titlesize=18)   # fontsize of the axes title
plt.rc('axes', labelsize=15)   # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)  # fontsize of the tick labels
plt.rc('ytick', labelsize=12)  # fontsize of the tick labels
plt.rc('legend', fontsize=16)  # legend fontsize


def main():
    filepath = "/Users/stephendiiorio/Desktop/asdf/"
    files = sorted(glob.glob(filepath + "*_heat_range.npz"))
    # print files
    print(len(files))

    heatmap = []

    for f in files:
        # print f
        data = np.load(f)
        r_heatmap = data['arr_0']
        r_heatmap = np.sum(r_heatmap, axis=0)
        heatmap.append(r_heatmap)
        data.close()

    heatmap = np.array(heatmap).T

    shape = heatmap.shape
    print(shape)
    extent = [0., 100., 5., -5.]

    fig, ax = plt.subplots()
    im = ax.imshow(heatmap, extent=extent, interpolation='none')
    plt.xlabel('t $(ps)$')
    plt.ylabel('x')
    plt.locator_params(axis='y', nbins=2)
    plt.title('Detected Electron Density')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, ticks=[np.amin(heatmap), np.amax(heatmap)])
    cbar.ax.set_yticklabels([str(np.amin(heatmap)), str(np.amax(heatmap))])

    plt.tight_layout()
    plt.savefig('asdf.png', dpi=600, bbox_inches="tight")

    return


if __name__ == "__main__":
    main()
