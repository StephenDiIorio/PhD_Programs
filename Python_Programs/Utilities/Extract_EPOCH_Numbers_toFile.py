import glob

import numpy as np

import sdf


def get_files(wkdir='Data'):
    """Get a list of SDF filenames belonging to the same run"""
    flist = glob.glob(wkdir + "/[0-9]*.sdf")
    flist = sorted(flist)

    return flist


def clean_file_list(file_list, varname):
    new_file_list = []

    for f in file_list:
        try:
            data = sdf.read(f)
            data.__dict__[varname]
            new_file_list.append(f)
        except KeyError:
            pass
    return new_file_list


def main():
    directory = '/scratch/lsa_flux/diiorios/2d_run/'
    ex_varname = 'Electric_Field_Ex'
    ey_varname = 'Electric_Field_Ey'
    ez_varname = 'Electric_Field_Ez'

    file_list = get_files(wkdir=directory)
    e_file_list = clean_file_list(file_list, ex_varname)

    e_file_list.remove(directory + '0000.sdf')
    e_file_list.remove(directory + '0002.sdf')
    e_file_list.remove(directory + '0003.sdf')
    e_file_list.remove(directory + '0004.sdf')
    e_file_list.remove(directory + '0005.sdf')
    e_file_list.remove(directory + '0006.sdf')
    e_file_list.remove(directory + '0007.sdf')
    e_file_list.remove(directory + '0008.sdf')
    e_file_list.remove(directory + '0009.sdf')
    e_file_list.remove(directory + '0011.sdf')
    e_file_list.remove(directory + '0012.sdf')
    e_file_list.remove(directory + '0013.sdf')
    e_file_list.remove(directory + '0014.sdf')
    e_file_list.remove(directory + '0015.sdf')
    e_file_list.remove(directory + '0016.sdf')
    e_file_list.remove(directory + '0017.sdf')
    e_file_list.remove(directory + '0018.sdf')
    e_file_list.remove(directory + '0019.sdf')
    e_file_list.remove(directory + '0021.sdf')
    e_file_list.remove(directory + '0022.sdf')
    e_file_list.remove(directory + '0023.sdf')
    e_file_list.remove(directory + '0024.sdf')
    e_file_list.remove(directory + '0025.sdf')
    e_file_list.remove(directory + '0026.sdf')
    e_file_list.remove(directory + '0027.sdf')
    e_file_list.remove(directory + '0028.sdf')
    e_file_list.remove(directory + '0029.sdf')
    e_file_list.remove(directory + '0031.sdf')
    e_file_list.remove(directory + '0032.sdf')
    e_file_list.remove(directory + '0033.sdf')
    e_file_list.remove(directory + '0034.sdf')
    e_file_list.remove(directory + '0035.sdf')
    e_file_list.remove(directory + '0036.sdf')
    e_file_list.remove(directory + '0037.sdf')
    e_file_list.remove(directory + '0038.sdf')
    e_file_list.remove(directory + '0039.sdf')
    e_file_list.remove(directory + '0041.sdf')
    e_file_list.remove(directory + '0042.sdf')
    e_file_list.remove(directory + '0043.sdf')
    e_file_list.remove(directory + '0044.sdf')
    e_file_list.remove(directory + '0045.sdf')
    e_file_list.remove(directory + '0046.sdf')
    e_file_list.remove(directory + '0047.sdf')
    e_file_list.remove(directory + '0048.sdf')
    e_file_list.remove(directory + '0049.sdf')
    e_file_list.remove(directory + '0051.sdf')
    e_file_list.remove(directory + '0052.sdf')
    e_file_list.remove(directory + '0053.sdf')
    e_file_list.remove(directory + '0054.sdf')
    e_file_list.remove(directory + '0055.sdf')
    e_file_list.remove(directory + '0056.sdf')
    e_file_list.remove(directory + '0057.sdf')
    e_file_list.remove(directory + '0058.sdf')
    e_file_list.remove(directory + '0059.sdf')
    e_file_list.remove(directory + '0061.sdf')
    e_file_list.remove(directory + '0062.sdf')
    e_file_list.remove(directory + '0063.sdf')
    e_file_list.remove(directory + '0064.sdf')
    e_file_list.remove(directory + '0065.sdf')
    e_file_list.remove(directory + '0066.sdf')
    e_file_list.remove(directory + '0067.sdf')
    e_file_list.remove(directory + '0068.sdf')
    e_file_list.remove(directory + '0069.sdf')
    e_file_list.remove(directory + '0071.sdf')
    e_file_list.remove(directory + '0072.sdf')
    e_file_list.remove(directory + '0073.sdf')
    e_file_list.remove(directory + '0074.sdf')
    e_file_list.remove(directory + '0075.sdf')
    e_file_list.remove(directory + '0076.sdf')
    e_file_list.remove(directory + '0077.sdf')
    e_file_list.remove(directory + '0078.sdf')
    e_file_list.remove(directory + '0079.sdf')
    e_file_list.remove(directory + '0081.sdf')
    e_file_list.remove(directory + '0082.sdf')
    e_file_list.remove(directory + '0083.sdf')
    e_file_list.remove(directory + '0084.sdf')
    e_file_list.remove(directory + '0085.sdf')
    e_file_list.remove(directory + '0086.sdf')
    e_file_list.remove(directory + '0087.sdf')
    e_file_list.remove(directory + '0088.sdf')
    e_file_list.remove(directory + '0089.sdf')
    e_file_list.remove(directory + '0091.sdf')
    e_file_list.remove(directory + '0092.sdf')
    e_file_list.remove(directory + '0093.sdf')
    e_file_list.remove(directory + '0094.sdf')
    e_file_list.remove(directory + '0095.sdf')
    e_file_list.remove(directory + '0096.sdf')
    e_file_list.remove(directory + '0097.sdf')
    e_file_list.remove(directory + '0098.sdf')
    e_file_list.remove(directory + '0099.sdf')

    print('Found {} files to plot'.format(len(file_list)))

    for i, file in enumerate(e_file_list):
        data = sdf.read(file)

        with open('x' + str(i) +'_field.dat', mode='w') as f:
            data.__dict__[ex_varname].data.astype('float64', order='F').tofile(f, sep='\t')
        with open('y' + str(i) + '_field.dat', mode='w') as f:
            data.__dict__[ey_varname].data.astype('float64', order='F').tofile(f, sep='\t')
        with open('z' + str(i) + '_field.dat', mode='w') as f:
            data.__dict__[ez_varname].data.astype('float64', order='F').tofile(f, sep='\t')

    return


if __name__ == "__main__":
    main()
