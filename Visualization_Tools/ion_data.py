import numpy as np
from scipy.integrate import simps

import sdf


def main():
    path = "/Users/stephendiiorio/Documents/Research/ion_test/"

    e_num = "0001"
    e_name = path + e_num + ".sdf"

    sdfdata = sdf.read(e_name)
    electron_data = sdfdata.Derived_Number_Density_Electron
    e_num = electron_data.data
    e_grid = electron_data.grid_mid.data[0]
    int_electron = simps(e_num, e_grid)
    print(int_electron)
    e_count = np.array([int_electron])

    ex_num = "0000"
    ex_name = path + ex_num + ".sdf"

    sdfdata = sdf.read(ex_name)
    ex_data = sdfdata.Electric_Field_Ex
    ex_amp = ex_data.data[0]
    print(ex_amp)
    e_field = np.array([ex_amp])

    np.savez(str(int(ex_amp)), e_count, e_field)

    return


if __name__ == "__main__":
    main()
