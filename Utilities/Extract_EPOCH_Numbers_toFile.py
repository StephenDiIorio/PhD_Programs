import numpy as np

import sdf


def main(file_name, varname, data_index, output_base):
    data = sdf.read(file_name)

    with open(output_base + '_numdens.dat', mode='wb') as f:
        data.__dict__[varname].data[data_index].astype('float64', order='F').tofile(f)
    with open(output_base + '_temp.dat', mode='wb') as f:
        data.__dict__[varname].data[data_index].astype('float64', order='F').tofile(f)

    return
