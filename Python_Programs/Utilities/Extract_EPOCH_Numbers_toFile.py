import numpy as np

import sdf


def main():
    directory = '/scratch/lsa_flux/diiorios/2d_run/'
    file_name = '0149.sdf'
    ex_varname = 'Electric_Field_Ex'
    ey_varname = 'Electric_Field_Ey'
    ez_varname = 'Electric_Field_Ez'

    data = sdf.read(directory + file_name)

    with open('xfield.dat', mode='w') as f:
        data.__dict__[ex_varname].data.astype('float64', order='F').tofile(f, sep='\t')
    with open('yfield.dat', mode='w') as f:
        data.__dict__[ey_varname].data.astype('float64', order='F').tofile(f, sep='\t')
    with open('zfield.dat', mode='w') as f:
        data.__dict__[ez_varname].data.astype('float64', order='F').tofile(f, sep='\t')

    return


if __name__ == "__main__":
    main()
