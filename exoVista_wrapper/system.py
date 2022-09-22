import numpy as np
from astropy.io.fits import getheader

from exoVista_wrapper.systemobject import SystemObject


class System:
    """
    Class for the whole stellar system
    """

    def __init__(self, infile):
        self.file = infile

        # fits file extensions, exoVista hard codes these
        star_ext = 3
        planet_ext = 4

        # Get the number of planets
        with open(infile, "rb") as f:
            # read header of first extension
            h = getheader(f, ext=0, memmap=False)
        n_ext = h["N_EXT"]  # get the largest extension
        nplanets = n_ext - 3

        # Create star object
        self.star = SystemObject(infile, star_ext)
        self.planets = []
        # loop over all planets
        for i in range(nplanets):
            self.planets.append(SystemObject(infile, planet_ext + i))

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
