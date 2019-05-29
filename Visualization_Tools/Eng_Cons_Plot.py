from matplotlib import use
use('Agg')

import glob
import os.path
import sys

package_directory = os.path.dirname(
    os.path.abspath(__file__))  # Get path to current file
# Trace path back to Utilities folder to import modules
sys.path.insert(0, os.path.join(package_directory, os.pardir, 'Utilities'))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

import sdf
from PlottingTools import get_si_prefix


def get_files(wkdir='Data', base=None):
    """Get a list of SDF filenames belonging to the same run"""
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

    return flist


def clean_file_list(file_list, varname):
    new_file_list = []

    for f in file_list:
        try:
            data = sdf.read(f, mmap=0)
            dummy = data.__dict__[varname]
            new_file_list.append(f)
        except KeyError:
            pass
    return new_file_list


def main():
    absorb_key = 'Absorption_Fraction_of_Laser_Energy_Absorbed____'
    inject_key = 'Absorption_Total_Laser_Energy_Injected__J_'
    field_en_key = 'Total_Field_Energy_in_Simulation__J_'
    part_en_key = 'Total_Particle_Energy_in_Simulation__J_'

    filepath = "/scratch/lsa_flux/diiorios/Exp_81418/energy_test"
    flist = get_files(wkdir=filepath)
    flist = clean_file_list(flist, absorb_key)
    flist = clean_file_list(flist, inject_key)

    time = []
    absorbed = []
    injected = []
    field = []
    part = []

    for f in flist:
        data = sdf.read(f, mmap=0)
        time.append(float(data.Header['time']))
        absorbed.append(data.__dict__[absorb_key].data)
        injected.append(data.__dict__[inject_key].data)
        field.append(data.__dict__[field_en_key].data)
        part.append(data.__dict__[part_en_key].data)

    tmult, tsym = get_si_prefix(np.max(time) - np.min(time))
    amult, asym = get_si_prefix(np.max(absorbed) - np.min(absorbed))
    imult, isym = get_si_prefix(np.max(injected) - np.min(injected))
    fmult, fsym = get_si_prefix(np.max(field) - np.min(field))
    pmult, psym = get_si_prefix(np.max(part) - np.min(part))

    fig = plt.figure()

    ax = fig.add_subplot(111)
    l1, = ax.plot(time, absorbed, 'k', label='Laser Energy Absorbed')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * tmult)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * amult)))
    ax.set_xlabel('Time $(' + tsym + 's)$')
    ax.set_ylabel('Absorption/Abs_frac')

    ax2 = ax.twinx()
    l2, = ax2.plot(time, injected, 'r', label='Laser Energy Injected')
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * imult)))
    ax2.set_ylabel('Absorption/Laser_enTotal $(' + isym + 'J)$')

    ax3 = ax.twinx()
    ax3.spines['right'].set_position(('outward', 60))
    l3, = ax3.plot(time, field, 'b', label='Total Field Energy')
    ax3.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * fmult)))
    ax3.set_ylabel('Total Field Energy $(' + fsym + 'J)$')

    ax4 = ax.twinx()
    ax4.spines['right'].set_position(('outward', 120))
    l4, = ax4.plot(time, part, 'g', label='Total Particle Energy')
    ax4.yaxis.set_major_formatter(FuncFormatter(lambda x, y: (x * pmult)))
    ax4.set_ylabel('Total Particle Energy $(' + psym + 'J)$')

    ax.yaxis.label.set_color(l1.get_color())
    ax2.yaxis.label.set_color(l2.get_color())
    ax3.yaxis.label.set_color(l3.get_color())
    ax4.yaxis.label.set_color(l4.get_color())

    ls = [l1, l2, l3, l4]
    labels = [l.get_label() for l in ls]
    ax.legend(ls, labels, loc='best')

    plt.tight_layout()
    plt.savefig('Eng_Cons.png')

    return


if __name__ == "__main__":
    main()
