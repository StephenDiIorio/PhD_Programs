import numpy as np
from scipy.constants import c, epsilon_0, e, m_e


__wp_coeff = c * 100.  # speed of light in [cm / s]
__n0_coeff = c / np.sqrt(100 * np.power(e, 2) / m_e / epsilon_0)  # coefficient terms in the equation for plasma frequency after converting n_e to cm^-3 and epsilon_0 to cm as well


#################################
# POSITION CONVERSIONS
#################################
def pos_si_to_osiris_wp(pos_si, wp):
    """
    Convert position to OSIRIS units normalized with a known plasma frequency.

    Parameters
    ----------
    pos_si : float
        The position in SI units [m] to convert.
    wp : float
        The plasma frequency [rad / s] to normalize against.

    Returns
    -------
    pos_osiris : numpy.float64
        Position in OSIRIS units normalized to the plasma frequency.

    """
    pos_si = np.float64(pos_si) * np.float64(100.)  # convert to [cm]
    wp = np.float64(wp)
    return pos_si * wp / __wp_coeff


def pos_osiris_to_si_wp(pos_osiris, wp):
    """
    Convert position to SI units from OSIRIS units normalized to a known
    plasma frequency.

    Parameters
    ----------
    pos_osiris : float
        The position in OSIRIS units to convert.
    wp : float
        The plasma frequency [rad / s] the OSIRIS value was normalized
        against.

    Returns
    -------
    pos_si : numpy.float64
        Position in SI units [m].

    """
    pos_osiris = np.float64(pos_osiris) / np.float64(100.)  # convert to [m]
    wp = np.float64(wp)
    return pos_osiris * __wp_coeff / wp


def pos_si_to_osiris_n0(pos_si, n0):
    """
    Convert position to OSIRIS units normalized with a known number density.

    Parameters
    ----------
    pos_si : float
        The position in SI units [m] to convert.
    n0 : float
        The number density [m^-3] to normalize against.

    Returns
    -------
    pos_osiris : numpy.float64
        Position in OSIRIS units normalized to the number density.

    """
    pos_si = np.float64(pos_si) * np.float64(100.)  # convert to [cm]
    n0 = np.float64(n0) * np.float64(1.e-6)  # convert to [cm^3]
    return pos_si * np.sqrt(n0) / __n0_coeff


def pos_osiris_to_si_n0(pos_osiris, n0):
    """
    Convert position to SI units from OSIRIS units normalized to a known
    number density.

    Parameters
    ----------
    pos_osiris : float
        The position in OSIRIS units to convert.
    n0 : float
        The number density [m^-3] the OSIRIS value was normalized against.

    Returns
    -------
    pos_si : numpy.float64
        Position in SI units [m].

    """
    pos_osiris = np.float64(pos_osiris) / np.float64(100.)  # convert to [m]
    n0 = np.float64(n0) * np.float64(1.e-6)  # convert to [cm^3]
    return pos_osiris * __n0_coeff / np.sqrt(n0)


#################################
# MOMENTUM CONVERSIONS
#################################
def mom_si_to_osiris(mom_si, wp):
    """
    """
    return


def mom_osiris_to_si(mom_osiris, wp):
    """
    """
    return


#################################
# ELECTRIC FIELD CONVERSIONS
#################################
__ef_wp_coeff = m_e * c / e / 1.e11  # correct terms and conversion into GV / cm
__ef_n0_coeff = np.float64(9.613e-10)


def efield_si_to_osiris_wp(efield_si, wp):
    """
    Convert electric field amplitude to OSIRIS units normalized with a known
    plasma frequency.

    Parameters
    ----------
    efield_si : float
        The electric field amplitude in SI units [V / m] to convert.
    wp : float
        The plasma frequency [rad / s] to normalize against.

    Returns
    -------
    efield_osiris : numpy.float64
        Electric field amplitude in OSIRIS units normalized to the plasma
        frequency.

    """
    efield_si = np.float64(efield_si) / np.float64(1.e11)  # convert to GV/cm
    wp = np.float64(wp)
    return efield_si / __ef_wp_coeff / wp


def efield_osiris_to_si_wp(efield_osiris, wp):
    """
    Convert electric field amplitude to SI units from OSIRIS units normalized
    to a known plasma frequency.

    Parameters
    ----------
    efield_osiris : float
        The electric field amplitude in OSIRIS units to convert.
    wp : float
        The plasma frequency [rad / s] the OSIRIS value was normalized
        against.

    Returns
    -------
    efield_si : numpy.float64
        Electric field amplitude in SI units [V / m].

    """
    efield_osiris = np.float64(efield_osiris) * np.float64(1.e11)  # convert to V/m
    wp = np.float64(wp)
    return efield_osiris * __ef_wp_coeff * wp


def efield_si_to_osiris_n0(efield_si, n0):
    """
    Convert electric field amplitude to OSIRIS units normalized with a known
    number density.

    Parameters
    ----------
    efield_si : float
        The electric field amplitude in SI units [V / m] to convert.
    n0 : float
        The number density [m^-3] to normalize against.

    Returns
    -------
    efield_osiris : numpy.float64
        Electric field amplitude in OSIRIS units normalized to the number
        density.

    """
    efield_si = np.float64(efield_si) / np.float64(1.e11)  # convert to GV/cm
    n0 = np.float64(n0) * np.float64(1.e-6)  # convert to [cm^3]
    return efield_si / __ef_n0_coeff / np.sqrt(n0)


def efield_osiris_to_si_n0(efield_osiris, n0):
    """
    Convert electric field amplitude to SI units from OSIRIS units normalized
    to a known number density.

    Parameters
    ----------
    efield_osiris : float
        The electric field amplitude in OSIRIS units to convert.
    n0 : float
        The number density [m^-3] the OSIRIS value was normalized against.

    Returns
    -------
    efield_si : numpy.float64
        Electric field amplitude in SI units [V / m].

    """
    efield_osiris = np.float64(efield_osiris) * np.float64(1.e11)  # convert to V/m
    n0 = np.float64(n0) * np.float64(1.e-6)  # convert to [cm^3]
    return efield_osiris * __ef_n0_coeff * np.sqrt(n0)


#################################
# MAGNETIC FIELD CONVERSIONS
#################################
__bf_wp_coeff = np.float64(5.681e-8)
__bf_n0_coeff = np.float64(3.204e-3)


def bfield_si_to_osiris_wp(bfield_si, wp):
    """
    Convert magnetic field amplitude to OSIRIS units normalized with a known
    plasma frequency.

    Parameters
    ----------
    bfield_si : float
        The magnetic field amplitude in SI units [T] to convert.
    wp : float
        The plasma frequency [rad / s] to normalize against.

    Returns
    -------
    bfield_osiris : numpy.float64
        Magnetic field amplitude in OSIRIS units normalized to the plasma
        frequency.

    """
    bfield_si = np.float64(bfield_si) * np.float64(1.e4)  # convert to [G]
    wp = np.float64(wp)
    return bfield_si / __bf_wp_coeff / wp


def bfield_osiris_to_si_wp(bfield_osiris, wp):
    """
    Convert magnetic field amplitude to SI units from OSIRIS units normalized
    to a known plasma frequency.

    Parameters
    ----------
    bfield_osiris : float
        The magnetic field amplitude in OSIRIS units to convert.
    wp : float
        The plasma frequency [rad / s] the OSIRIS value was normalized
        against.

    Returns
    -------
    bfield_si : numpy.float64
        Magnetic field amplitude in SI units [T].

    """
    bfield_osiris = np.float64(bfield_osiris) / np.float64(1.e4)  # convert to [T]
    wp = np.float64(wp)
    return bfield_osiris * __bf_wp_coeff * wp


def bfield_si_to_osiris_n0(bfield_si, n0):
    """
    Convert magnetic field amplitude to OSIRIS units normalized with a known
    number density.

    Parameters
    ----------
    bfield_si : float
        The magnetic field amplitude in SI units [T] to convert.
    n0 : float
        The number density [m^-3] to normalize against.

    Returns
    -------
    bfield_osiris : numpy.float64
        Magnetic field amplitude in OSIRIS units normalized to the number
        density.

    """
    bfield_si = np.float64(bfield_si) * np.float64(1.e4)  # convert to [G]
    n0 = np.float64(n0) * np.float64(1.e-6)  # convert to [cm^3]
    return bfield_si / __bf_n0_coeff / np.sqrt(n0)


def bfield_osiris_to_si_n0(bfield_osiris, n0):
    """
    Convert magnetic field amplitude to SI units from OSIRIS units normalized
    to a known number density.

    Parameters
    ----------
    bfield_osiris : float
        The magnetic field amplitude in OSIRIS units to convert.
    n0 : float
        The number density [m^-3] the OSIRIS value was normalized against.

    Returns
    -------
    bfield_si : numpy.float64
        Magnetic field amplitude in SI units [T].

    """
    bfield_osiris = np.float64(bfield_osiris) / np.float64(1.e4)  # convert to [T]
    n0 = np.float64(n0) * np.float64(1.e-6)  # convert to [cm^3]
    return bfield_osiris * __bf_n0_coeff * np.sqrt(n0)
