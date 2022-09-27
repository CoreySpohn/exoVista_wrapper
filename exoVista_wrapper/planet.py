import astropy.constants as const
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io.fits import getdata
from astropy.time import Time
from keplertools import fun as kt


class Planet:
    """
    Class for the planets in the exoVista systems
    Assumes e=0, which is the case for all downloaded systems
    """

    def __init__(self, infile, fits_ext, star):

        # Get the object's data from the fits file
        with open(infile, "rb") as f:
            obj_data, obj_header = getdata(f, ext=fits_ext, header=True, memmap=False)

        # Time data
        self._t = obj_data[:, 0] * u.yr
        self.t0 = Time(self._t[0].to(u.yr).value, format="decimalyear")

        # Position data
        self._x = obj_data[:, 9] * u.AU
        self._y = obj_data[:, 10] * u.AU
        self._z = obj_data[:, 11] * u.AU

        # Velocity data
        self._vx = obj_data[:, 12] * u.AU / u.yr
        self._vy = obj_data[:, 13] * u.AU / u.yr
        self._vz = obj_data[:, 14] * u.AU / u.yr

        # Assign the planet's keplerian orbital elements
        self.a = obj_header["A"] * u.AU
        self.e = obj_header["E"]
        self.i = (obj_header["I"] * u.deg).to(u.rad)
        self.W = (obj_header["LONGNODE"] * u.deg).to(u.rad)
        # self.w = (obj_header["ARGPERI"] * u.deg).to(u.rad)
        self.w = 0 * u.rad

        # Assign the planet's mass/radius information
        self.mass = obj_header["M"] * u.M_earth
        self.radius = obj_header["R"] * u.R_earth

        # Assign the planet's time-varying mean anomaly, argument of pericenter,
        # true anomaly, and contrast
        self.rep_w = (obj_data[:, 7] * u.deg + 180 * u.deg) % (2 * np.pi * u.rad)
        self.M = (obj_data[:, 8] * u.deg + 180 * u.deg) % (2 * np.pi * u.rad)
        self.nu = (self.rep_w + self.M) % (
            2 * np.pi * u.rad
        )  # true anomaly for circular orbits
        self.contrast = obj_data[:, 15]

        # Initial mean anomaly
        self.M0 = self.nu[0]

        # Gravitational parameter
        self.mu = (const.G * (self.mass + star.mass)).decompose()

        # Mean angular motion
        self.n = (np.sqrt(self.mu / self.a**3)).decompose()

        # Propagation table
        self.vectors = pd.DataFrame(
            {
                "t": [self._t[0].decompose().value],
                "x": [self._x[0].decompose().value],
                "y": [self._y[0].decompose().value],
                "z": [self._z[0].decompose().value],
                "vx": [self._vx[0].decompose().value],
                "vy": [self._vy[0].decompose().value],
                "vz": [self._vz[0].decompose().value],
            }
        )

        self.star = star

    def calc_vectors(self, t, return_r=True, return_v=False):
        """
        Given a time, calculate the planet's position and/or velocity vectors
        Args:
            t (Time):
                Time to calculate the position vectors at
        Returns:
            r(ndarray):
                3 x n stacked position vector
            v(ndarray):
                3 x n stacked velocity vector

        """
        # This will find the radial and velocity vectors at an epoch
        M = self.mean_anom(t)
        E = kt.eccanom(M.value, self.e)
        a, e, Omega, inc, w = self.a.decompose(), self.e, self.W, self.i, self.w
        if not np.isscalar(E):
            a = np.ones(len(E)) * a
            e = np.ones(len(E)) * e
            Omega = np.ones(len(E)) * Omega
            inc = np.ones(len(E)) * inc
            w = np.ones(len(E)) * w
        A = np.vstack(
            (
                a
                * (np.cos(Omega) * np.cos(w) - np.sin(Omega) * np.cos(inc) * np.sin(w)),
                a
                * (np.sin(Omega) * np.cos(w) + np.cos(Omega) * np.cos(inc) * np.sin(w)),
                a * np.sin(inc) * np.sin(w),
            )
        )

        B = np.vstack(
            (
                -a
                * np.sqrt(1 - e**2)
                * (np.cos(Omega) * np.sin(w) + np.sin(Omega) * np.cos(inc) * np.cos(w)),
                a
                * np.sqrt(1 - e**2)
                * (
                    -np.sin(Omega) * np.sin(w) + np.cos(Omega) * np.cos(inc) * np.cos(w)
                ),
                a * np.sqrt(1 - e**2) * np.sin(inc) * np.cos(w),
            )
        )

        # Calculate position vectors
        if return_r:
            if np.isscalar(self.mu) and not (np.isscalar(E)):
                r = np.matmul(A, np.array((np.cos(E) - e), ndmin=2)) + np.matmul(
                    B, np.array(np.sin(E), ndmin=2)
                )
            else:
                r = np.matmul(A, np.diag(np.cos(E) - e)) + np.matmul(
                    B, np.diag(np.sin(E))
                )

        # Calculate velocity vectors
        if return_v:
            if np.isscalar(self.mu) and not (np.isscalar(E)):
                v = (
                    np.matmul(-A, np.array(np.sin(E), ndmin=2))
                    + np.matmul(B, np.array(np.cos(E), ndmin=2))
                ) * np.tile(
                    np.sqrt(self.mu * a ** (-3.0)) / (1 - e * np.cos(E)), (3, 1)
                )
            else:
                v = np.matmul(
                    np.matmul(-A, np.diag(np.sin(E)))
                    + np.matmul(B, np.diag(np.cos(E))),
                    np.diag(np.sqrt(self.mu * a ** (-3.0)) / (1 - e * np.cos(E))),
                )

        # Rotate around x axis with midplane inclination
        rot1 = np.array(
            [
                [1, 0, 0],
                [
                    0,
                    np.cos(self.star.midplane_I.to(u.rad)),
                    -np.sin(self.star.midplane_I.to(u.rad)),
                ],
                [
                    0,
                    np.sin(self.star.midplane_I.to(u.rad)),
                    np.cos(self.star.midplane_I.to(u.rad)),
                ],
            ]
        )

        # Rotate around z axis with midplane position angle
        rot2 = np.array(
            [
                [
                    np.cos(self.star.midplane_PA.to(u.rad)),
                    np.sin(self.star.midplane_PA.to(u.rad)),
                    0,
                ],
                [
                    -np.sin(self.star.midplane_PA.to(u.rad)),
                    np.cos(self.star.midplane_PA.to(u.rad)),
                    0,
                ],
                [0, 0, 1],
            ]
        )
        if return_r:
            r = np.matmul(r.T, rot1)
            r = np.matmul(r, rot2)
        if return_v:
            v = np.matmul(v.T, rot1)
            v = np.matmul(v, rot2)

        if return_r and return_v:
            return r, v
        if return_r:
            return r
        if return_v:
            return v

    def calc_vs(self, t, return_nu=False):
        """
        Calculate the radial velocities

        Args:
            t (astropy Time):
                input time
            return_nu (Bool):
                Whether the true anomaly should be returned

        Returns:
            vs (Astropy quantity):
                The radial velocity at time t
            nu (astropy Quantity):
                True anomaly
        """
        M = self.mean_anom(t)
        E = kt.eccanom(M, self.e)
        nu = kt.trueanom(E, self.e) * u.rad
        vs = (
            np.sqrt(
                const.G / ((self.mass + self.star.mass) * self.a * (1 - self.e**2))
            )
            * self.mass
            * np.sin(self.i)
            * (np.cos(self.w + nu) + self.e * np.cos(self.w))
        )
        if return_nu:
            return vs.decompose(), nu
        else:
            return vs.decompose()

    def mean_anom(self, t):
        """
        Calculate the mean anomaly at a given time, assumes initial time is 0
        Args:
            t (Time):
                Time to calculate mean anomaly

        Returns:
            M (astropy Quantity):
                Planet's mean anomaly at t (radians)
        """
        # t is a Time quantity
        # n is in standard SI units so the times are converted to seconds
        M1 = (self.n * t).decompose() * u.rad
        M = ((M1 + self.M0).to(u.rad)) % (2 * np.pi * u.rad)
        return M
