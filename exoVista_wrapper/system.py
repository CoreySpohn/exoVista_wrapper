import astropy.constants as const
import astropy.units as u
import numpy as np
import pandas as pd
import rebound
from astropy.io.fits import getheader
from tqdm import tqdm

from exoVista_wrapper.planet import Planet
from exoVista_wrapper.star import Star


class System:
    """
    Class for the whole stellar system
    """

    def __init__(self, infile):
        self.file = infile

        # fits file extensions, exoVista hard codes these
        planet_ext = 4

        # Get the number of planets
        with open(infile, "rb") as f:
            # read header of first extension
            h = getheader(f, ext=0, memmap=False)
        n_ext = h["N_EXT"]  # get the largest extension
        nplanets = n_ext - 3

        # Create star object
        self.star = Star(infile)
        self.planets = []
        # loop over all planets
        for i in range(nplanets):
            self.planets.append(Planet(infile, planet_ext + i, self.star))

        # Set up rebound simulation
        self.sim = rebound.Simulation()
        self.sim.G = const.G.value
        self.sim.add(
            m=self.star.mass.decompose().value,
            x=self.star._x[0].decompose().value,
            y=self.star._y[0].decompose().value,
            z=self.star._z[0].decompose().value,
            vx=self.star._vx[0].decompose().value,
            vy=self.star._vy[0].decompose().value,
            vz=self.star._vz[0].decompose().value,
        )
        for planet in self.planets:
            self.sim.add(
                m=planet.mass.decompose().value,
                x=planet._x[0].decompose().value,
                y=planet._y[0].decompose().value,
                z=planet._z[0].decompose().value,
                vx=planet._vx[0].decompose().value,
                vy=planet._vy[0].decompose().value,
                vz=planet._vz[0].decompose().value,
            )
        self.sim.move_to_com()

    def get_rv(self, t):
        """
        Calculate the radial velocity of the system

        Args:
            t (Astropy Time array):
                Astropy time quantities that RV values are desired for

        Returns:
            rv (Astropy Quantity array):
                Star's radial velocity in z direction
        """
        rv = np.zeros(len(t)) * u.m / u.s
        for planet in self.planets:
            pv = planet.calc_vectors(t, return_r=True)
            x, y, z = pv[:, 0], pv[:, 1], pv[:, 2]
            vz = planet.calc_vs(t)
            rv += vz
        return rv

    def propagate_system(self, t):
        """
        Propage system with rebound
        """
        for time in tqdm(t, desc="System propagation"):
            self.sim.integrate(time)
            for j, p in enumerate(self.sim.particles):
                p_vectors = {
                    "t": [time],
                    "x": [p.x],
                    "y": [p.y],
                    "z": [p.z],
                    "vx": [p.vx],
                    "vy": [p.vy],
                    "vz": [p.vz],
                }
                if j == 0:
                    self.star.vectors = pd.concat(
                        [self.star.vectors, pd.DataFrame(p_vectors)]
                    )
                else:
                    planet = self.planets[j - 1]
                    planet.vectors = pd.concat(
                        [planet.vectors, pd.DataFrame(p_vectors)]
                    )
        self.star.vectors = self.star.vectors.sort_values("t").reset_index()
        for planet in self.planets:
            planet.vectors = planet.vectors.sort_values("t").reset_index()
