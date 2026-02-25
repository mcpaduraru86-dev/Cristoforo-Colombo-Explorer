import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


def resolve_name(name: str):
    try:
        coord = SkyCoord.from_name(name)
        return coord.ra.deg, coord.dec.deg
    except Exception:
        return None


def build_cosmology():
    return FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)


def radec_z_to_xyz(ra_deg, dec_deg, z, cosmo=None):
    if cosmo is None:
        cosmo = build_cosmology()
    Dc = cosmo.comoving_distance(z).value
    ra_rad = np.deg2rad(ra_deg)
    dec_rad = np.deg2rad(dec_deg)
    x = Dc * np.cos(dec_rad) * np.cos(ra_rad)
    y = Dc * np.cos(dec_rad) * np.sin(ra_rad)
    z3d = Dc * np.sin(dec_rad)
    return x, y, z3d
