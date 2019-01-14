from scipy.constants import physical_constants
import sdf_helper as sh
import numpy as np
from matplotlib import use, rc
from matplotlib.ticker import FuncFormatter
use('Agg')
rc('mathtext', default='regular')
import matplotlib.pyplot as plt
# plt.ion()


def get_si_prefix(scale):
    scale = abs(scale)
    mult = 1
    sym = ''

    scale = scale * mult
    if scale <= 0:
        pwr = 0
    else:
        pwr = (-np.floor(np.log10(scale)))
    mult = mult * np.power(10.0, pwr)
    if np.rint(pwr) != 0:
        sym = "x10^{%.0f}" % (-pwr) + sym

    return mult, sym


if __name__ == "__main__":
    path = "/Users/stephendiiorio/Desktop/"
    f_num = "0015"
    f_name = path + f_num + ".sdf"

    sdfdata = sh.getdata(f_name)
    electron_en = sdfdata.Particles_Ek_Electron
    electron_en_data = np.copy(electron_en.data)
    electron_en_data *= physical_constants["joule-electron volt relationship"][0]

    xmult, xsym = get_si_prefix(np.max(electron_en_data) - np.min(electron_en_data))

    fig = plt.figure()
    ax = fig.gca()
    ax.hist(electron_en_data[electron_en_data != 0], bins=100, normed=True)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * xmult)))
    ax.set_xlabel(r'Electron Energy $(' + xsym + 'eV)$')

    # plt.show()
    plt.savefig('spectrum.png')
