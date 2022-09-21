import numpy as np
from astropy.io.fits import getheader, info

from exoVista_wrapper.systemobject import SystemObject


class System:
    """
    Class for the whole stellar system
    """

    def __init__(self, infile):
        star_ext = 3
        planet_ext = 4
        system_info = info(infile, output=False)
        with open(infile, "rb") as f:
            # read header of first extension
            h = getheader(f, ext=0, memmap=False)
        n_ext = h["N_EXT"]  # get the largest extension
        nplanets = n_ext - 3
        if len(system_info) != n_ext + 1:
            # Catch error where part of the fits doesn't load
            return None

        self.star = SystemObject(infile, star_ext, 0, 0)
        self.planets = []
        # loop over all planets
        for i in range(nplanets):
            self.planets.append(SystemObject(infile, planet_ext + i, nplanets, i))

    def get_rv(self):
        """
        Return the RV data for the star and the planets
        """
        star_rv = np.array(self.star.vz)
        planet_rv = np.zeros((len(self.star.t), len(self.planets)))
        for i, planet in enumerate(self.planets):
            planet_rv[:, i] = planet.vz

        star_rv = star_rv * self.star.vz.unit
        planet_rv = planet_rv * self.planets[0].vz.unit
        return star_rv, planet_rv
