#########################################################
# MAIN FUNCTIONS TO CALL
#########################################################


def osirisEzSlice(dump, theta=0., e0=1.0, modes=[0, 1, 2], path='./', x0=1.0, use2r=False):
    """
    Returns the longitudinal component (z) of the electric field

    dump: output number to gather data from
    theta: angle to take slice of field from
    e0: electric field amplitude to normalize to
    modes: list of modenumbers to construct field from
    path: path to the MS folder
    x0: reference distance to scale length units by
    use2r: whether or not to return the grid and data ranging from rmin to rmax
           or from -rmax to rmax

    returns the r and z grids and Ez values, in that order
    """

    return __osirisLoadECylSlice(dump, '1', theta=theta, e0=e0, modes=modes, path=path, x0=x0, use2r=use2r)


def osirisErSlice(dump, theta=0., e0=1.0, modes=[0, 1, 2], path='./', x0=1.0, use2r=False):
    """
    Returns the tangential component (r) of the electric field

    dump: output number to gather data from
    theta: angle to take slice of field from
    e0: electric field amplitude to normalize to
    modes: list of modenumbers to construct field from
    path: path to the MS folder
    x0: reference distance to scale length units by
    use2r: whether or not to return the grid and data ranging from rmin to rmax
           or from -rmax to rmax

    returns the r and z grids and Er values, in that order
    """

    return __osirisLoadECylSlice(dump, '2', theta=theta, e0=e0, modes=modes, path=path, x0=x0, use2r=use2r)


def osirisExySlice(dump, theta=0., e0=1.0, modes=[0, 1, 2], path='./', x0=1.0, use2r=False):
    """
    Converts the cylindrical electric field into cartesian coordinates and returns the x and y components

    dump: output number to gather data from
    theta: angle to take slice of field from
    e0: electric field amplitude to normalize to
    modes: list of modenumbers to construct field from
    path: path to the MS folder
    x0: reference distance to scale length units by
    use2r: whether or not to return the grid and data ranging from xmin to xmax
           or from -xmax to xmax

    returns the x and y grids and the Ex and Ey values, in that order
    """

    Er_modes = __osirisLoadEField(dump, component='2', e0=e0, modes=modes, path=path)
    Et_modes = __osirisLoadEField(dump, component='3', e0=e0, modes=modes, path=path)
    r, z, t = __osirisLoadEFieldGrid(dump=dump, x0=x0, component='2', mode=modes[0], path=path)
    x, Ex, Ey = __EXYSlices(Er_modes, Et_modes, r, theta, use2r=use2r)

    return t, z, x, Ex, Ey


def osirisDensitySlice(dump, theta=0., n0=1., species='plasma', modes=[0, 1, 2], path='./', x0=1.0, use2r=False):
    import numpy as np

    """
    Returns a slice of the density profile at a given angle for a given species

    dump: output number to gather data from
    theta: angle to take slice of field from
    n0: density to normalize to
    species: name of species to gather data from
    modes: list of modenumbers to construct field from
    path: path to the MS folder
    x0: reference distance to scale length units by
    use2r: whether or not to return the grid and data ranging from rmin to rmax
           or from -rmax to rmax

    returns the r and z grids and Ez values, in that order
    """

    return __osirisLoadDensCylSlice(dump, theta=theta, n0=n0, species=species, modes=modes, path=path, x0=x0, use2r=use2r)


def osirisLoadGrid(dump, x0=1., species='plasma', path='./', use2r=False):
    import h5py
    import numpy as np

    """
    Returns the grid of the simulation box.

    dump: output number to gather data from
    x0: reference distance to scale length units by
    species: name of species to gather data from
    path: path to the MS folder
    use2r: whether or not to return the grid and data ranging from rmin to rmax
           or from -rmax to rmax

    returns the r and z grids, in that order
    """

    f = h5py.File(path + 'MS/DENSITY/%s/MODE-0-RE/charge_cyl_m/charge_cyl_m-%s-0-re-%.6d.h5' % (species, species, dump), 'r')

    Nr, Nz = f['charge_cyl_m'].shape

    zmin, zmax = f['AXIS']['AXIS1'][:]
    rmin, rmax = f['AXIS']['AXIS2'][:]

    # dz = (zmax - zmin) / Nz
    dr = (rmax - rmin) / Nr

    z = x0 * np.linspace(zmin, zmax, Nz)  # + 0.5 * dz
    if use2r:
        r = x0 * np.linspace(-rmax, rmax, 2 * Nr)
    else:
        r = x0 * np.linspace(rmin, rmax, Nr) + 0.5 * dr

    return r, z


def osirisParticles(dump, species='plasma', path='./', x0=1.0):
    import h5py
    import numpy as np

    """
    Obtain all the raw data from a single Osiris dump

    field: output number to gather data from
    species: name of species to gather data from
    path: path to the MS folder
    x0: reference distance to scale length units by

    returns a dictionary of all the raw variables stored in an Osiris dump
    """

    f = h5py.File(path + 'MS/RAW/%s/RAW-%s-%.6d.h5' % (species, species, dump), 'r')

    ene = f['ene'][()]
    p1 = f['p1'][()]
    p2 = f['p2'][()]
    p3 = f['p3'][()]
    wgt = f['q'][()]
    z = f['x1'][()] * x0
    r = f['x2'][()] * x0
    x = f['x3'][()] * x0
    y = f['x4'][()] * x0

    return {'z': z, 'r': r, 'x': x, 'y': y, 'p1': p1, 'p2': p2, 'p3': p3, 'ene': ene, 'wgt': wgt}


def osirisToVTK(theta, r, z, dump,
                diagnostics=['charge', 'e1', 'e2', 'e3'],
                modes=[0, 1, 2],
                n0=1.,
                downsampling=(slice(None, None, None), slice(None, None, None), slice(None, None, None)),
                path='./',
                output_filename='paraview'):

    import numpy as np
    from evtk.hl import gridToVTK
    from scipy import constants

    """
    Save osiris grid diagnostics to VTK array

    INPUTS:
    =======
    theta = linear array of azimuthal grid points, len(theta) = Ntheta
    r     = linear array of radial grid points   , len(r)     = Nr
    z     = linear array of axial grid points    , len(z)     = Nz

    dump  = integer simulation dump identifier

    path           : string, where to find the simulation outputs
    output_filename: string, where to save the VTK structured arrays


    diagnostics: list of strings indicating the field diagnostics to be saved
        example: diagnostics = [ 'charge' , 'e1' , 'e2' , 'e3' ]

    downsampling: triple of python slice objects for (theta,r,z)
    """

    # physical scale for the electric field
    eref = np.sqrt(constants.electron_mass * constants.c**2 / constants.epsilon_0) * 1000.
    e0 = eref * np.sqrt(n0)  # reference electric field strength in [V/m]
                             # n0 in 1/cm^3

    # turn linspaces into meshgrid
    thetagrid,rgrid,zgrid = np.meshgrid(theta[downsampling[0]], r[downsampling[1]], z[downsampling[2]], indexing='ij')

    xgrid = rgrid * np.cos(thetagrid)
    ygrid = rgrid * np.sin(thetagrid)

    # print(xgrid.shape)
    # print(ygrid.shape)

    # initialize an empty dictionary for the reconstructed field diagnostics
    reconstructed_diags = {}

    for diag in diagnostics:

        if diag == 'density':
            rho = __osirisLoadDensity(dump=dump, n0=n0, modes=modes, path=path)
            rec_rho = __basicFieldReconstructionFromNodes(rho[downsampling], thetagrid)

            reconstructed_diags['density (10^18 cm^-3)'] = np.nan_to_num(rec_rho) / 1.e18

        if diag == 'e1':
            e1 = __osirisLoadEField(dump=dump, component='1', e0=e0, modes=modes, path=path)
            rec_e1 = __basicFieldReconstructionFromNodes(e1[downsampling], thetagrid)
            reconstructed_diags['Ez (GV/m)'] = rec_e1 / 1e9
            del e1

        if diag == 'e2':
            e2 = __osirisLoadEField(dump=dump, component='2', e0=e0, modes=modes, path=path)
            rec_e2 = __basicFieldReconstructionFromNodes(e2[downsampling], thetagrid)
            reconstructed_diags['Er (TV/m)'] = rec_e2 / 1e12
            del e2

        if diag == 'e3':
            e3 = __osirisLoadEField(dump=dump, component='3', e0=e0, modes=modes, path=path)
            rec_e3 = __basicFieldReconstructionFromNodes(e3[downsampling], thetagrid)
            reconstructed_diags['Et (TV/m)'] = rec_e3 / 1e12
            del e3

    rec_ex = rec_e2 * np.cos(thetagrid) - rec_e3 * np.sin(thetagrid)
    rec_ey = rec_e2 * np.sin(thetagrid) + rec_e3 * np.cos(thetagrid)

    reconstructed_diags['Ex (TV/m)'] = rec_ex / 1e12
    reconstructed_diags['Ey (TV/m)'] = rec_ey / 1e12

    # print('output dump              :', dump)
    # print('save file to             :', output_filename + '%.6d' % dump)
    # print('full grid size           :', (len(theta), len(r), len(z)))
    # print('reduced grid size        :', xgrid.shape)
    # print('save structured grids for:', reconstructed_diags.keys())

    # print(xgrid.shape)
    # print(ygrid.shape)
    # print(zgrid.shape)
    # for k in reconstructed_diags.keys(): print(k, reconstructed_diags[k].shape)

    gridToVTK(output_filename + '%.6d' % dump,
              xgrid,
              ygrid,
              zgrid,
              pointData=reconstructed_diags)

    return


#########################################################
# PRIVATE FUNCTIONS
#########################################################


def __osirisLoadDensity(dump, n0=1., species='plasma', modes=[0, 1, 2], path='./'):
    import h5py
    import numpy as np

    loaded_modes = []

    for M in modes:
        ne_RE = h5py.File(path + 'MS/DENSITY/%s/MODE-%d-RE/charge_cyl_m/charge_cyl_m-%s-%d-re-%.6d.h5' % (species, M, species, M, dump), 'r')['charge_cyl_m']
        loaded_modes.append(-n0 * ne_RE[:])

        if M > 0:
            ne_IM = h5py.File(path + 'MS/DENSITY/%s/MODE-%d-IM/charge_cyl_m/charge_cyl_m-%s-%d-im-%.6d.h5' % (species, M, species, M, dump), 'r')['charge_cyl_m']
            loaded_modes.append(-n0 * ne_IM[:])

    return np.asarray(loaded_modes)


def __osirisLoadEField(dump, component='1', e0=1., modes=[0, 1, 2], path='./'):
    import h5py
    import numpy as np

    loaded_modes = []

    for M in modes:
        RE = h5py.File(path + 'MS/FLD/MODE-%d-RE/e%s_cyl_m/e%s_cyl_m-%d-re-%.6d.h5' % (M, component, component, M, dump), 'r')['e%s_cyl_m' % component]
        loaded_modes.append(e0 * RE[:])

        if M > 0:
            IM = h5py.File(path + 'MS/FLD/MODE-%d-IM/e%s_cyl_m/e%s_cyl_m-%d-im-%.6d.h5' % (M, component, component, M, dump), 'r')['e%s_cyl_m' % component]
            loaded_modes.append(e0 * IM[:])

    return np.asarray(loaded_modes)


def __basicFieldReconstructionFromNodes(field, theta):
    import numpy as np

    """
    Reconstruct 3D diagnostic on a cylindrical grid

    field: numpy array of the shape (2*Nmodes + 1, Nr, Nz)
    theta: scalar or numpy array of the shape (Ntheta, Nr, Nz) or (Ntheta, newaxis, newaxis)

    returns the sum over all modes for given theta
    shape of output array is [(Ntheta)] + (Nr, Nz)
    """

    Nmodes = field.shape[0] // 2 + 1

    modefunctions = [1 + 0. * theta]

    for m in range(1, Nmodes):
        modefunctions.append(np.cos(m * theta))
        modefunctions.append(np.sin(m * theta))

    return np.sum([mode * field[i] for i, mode in enumerate(modefunctions)], axis=0)


def __fieldSlice(field, r, theta, use2r=False):
    import numpy as np

    f1 = __basicFieldReconstructionFromNodes(field, theta)

    if use2r:
        f2 = __basicFieldReconstructionFromNodes(field, theta + np.pi)

        x = np.concatenate([-r[::-1], r], axis=0)
        f = np.concatenate([f2[::-1], f1], axis=0)
        return x, f
    else:
        return r, f1


def __osirisLoadECylSlice(dump, component, theta=0., e0=1.0, modes=[0, 1, 2], path='./', x0=1.0, use2r=False):
    E_modes = __osirisLoadEField(dump, component=component, e0=e0, modes=modes, path=path)
    r, z, t = __osirisLoadEFieldGrid(dump, x0=x0, component=component, mode=modes[0], path=path)
    r, E_slice = __fieldSlice(E_modes, r, theta=theta, use2r=use2r)

    return t, r, z, E_slice


def __osirisLoadDensCylSlice(dump, theta=0., n0=1., species='plasma', modes=[0, 1, 2], path='./', x0=1.0, use2r=False):
    density_modes = __osirisLoadDensity(dump=dump, n0=n0, species=species, modes=modes, path=path)  # in 1/cm^3
    r, z, t = __osirisLoadDensityGrid(dump, x0=x0, species=species, path=path)
    r, dens_slice = __fieldSlice(density_modes, r, theta=theta, use2r=use2r)

    return t, r, z, dens_slice


def __EXYSlices(E_r, E_t, r, theta, use2r=False):
    import numpy as np

    E_r_pos = __basicFieldReconstructionFromNodes(E_r, theta)
    E_t_pos = __basicFieldReconstructionFromNodes(E_t, theta)

    E_x_pos = E_r_pos * np.cos(theta) - E_t_pos * np.sin(theta)
    E_y_pos = E_r_pos * np.sin(theta) + E_t_pos * np.cos(theta)

    if use2r:
        E_r_neg = __basicFieldReconstructionFromNodes(E_r, theta + np.pi)
        E_t_neg = __basicFieldReconstructionFromNodes(E_t, theta + np.pi)

        E_x_neg = E_r_neg * np.cos(theta + np.pi) - E_t_neg * np.sin(theta + np.pi)
        E_y_neg = E_r_neg * np.sin(theta + np.pi) + E_t_neg * np.cos(theta + np.pi)

        x = np.concatenate([-r[::-1], r], axis=0)
        E_x = np.concatenate([E_x_neg[::-1], E_x_pos], axis=0)
        E_y = np.concatenate([E_y_neg[::-1], E_y_pos], axis=0)
    else:
        x = r
        E_x = np.concatenate([E_x_pos], axis=0)
        E_y = np.concatenate([E_y_pos], axis=0)

    return x, E_x, E_y


def __osirisLoadDensityGrid(dump, x0=1., species='plasma', path='./'):
    import h5py
    import numpy as np

    f = h5py.File(path + 'MS/DENSITY/%s/MODE-0-RE/charge_cyl_m/charge_cyl_m-%s-0-re-%.6d.h5' % (species, species, dump), 'r')

    t = f.attrs['TIME']

    Nr, Nz = f['charge_cyl_m'].shape

    zmin, zmax = f['AXIS']['AXIS1'][:]
    rmin, rmax = f['AXIS']['AXIS2'][:]

    # dz = (zmax - zmin) / Nz
    dr = (rmax - rmin) / Nr

    z = x0 * np.linspace(zmin, zmax, Nz)  # + 0.5 * dz
    r = x0 * np.linspace(rmin, rmax, Nr) + 0.5 * dr

    return r, z, t


def __osirisLoadEFieldGrid(dump, x0=1., component='1', mode=0, path='./'):
    import h5py
    import numpy as np

    f = h5py.File(path + 'MS/FLD/MODE-%d-RE/e%s_cyl_m/e%s_cyl_m-%d-re-%.6d.h5' % (mode, component, component, mode, dump), 'r')

    t = f.attrs['TIME']

    Nr, Nz = f['e%s_cyl_m' % component].shape

    zmin, zmax = f['AXIS']['AXIS1'][:]
    rmin, rmax = f['AXIS']['AXIS2'][:]

    # dz = (zmax - zmin) / Nz
    dr = (rmax - rmin) / Nr

    z = x0 * np.linspace(zmin, zmax, Nz)  # + 0.5 * dz
    r = x0 * np.linspace(rmin, rmax, Nr) + 0.5 * dr

    return r, z, t
