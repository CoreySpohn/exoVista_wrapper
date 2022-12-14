import subprocess
from pathlib import Path

import astropy.io.fits as fits
import astropy.units as u
import numpy as np
import pandas as pd
from scipy import optimize
from tqdm import tqdm


def runcmd(cmd, verbose=False):
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)
    pass


def get_data(universes=np.arange(1, 13)):
    """
    This function gets all the exoVista data. It gets the csv file with the universe
    information and puts it in the "data/{universe_number}/target_database.csv".
    Then it goes through every url in the csv file and pulls that fits file into
    "data/{universe_number}/{file}.fits".
    """
    # Iterate over the different universes
    for n in tqdm(universes, position=0, desc="Universe", leave=False):
        universe_url = (
            "https://ckan-files.emac.gsfc.nasa.gov/"
            f"exovista/DEC21/{n}/target_database.csv"
        )
        Path(f"data/{n}/").mkdir(parents=True, exist_ok=True)
        if not Path(f"data/{n}/target_database.csv").exists():
            runcmd(f"wget --directory-prefix=data/{n} {universe_url}", verbose=False)
        df = pd.read_csv(f"data/{n}/target_database.csv", low_memory=False)
        for i in tqdm(
            np.arange(1, df.shape[0]), position=1, desc="System", leave=False
        ):
            # Get file url
            fit_url = df.at[i, "URL"]

            # Create file path
            file_path = Path(f"data/{n}/{fit_url.split('/')[-1]}")

            # If file doesn't exist then pull it
            if not file_path.exists():
                runcmd(
                    f"wget --directory-prefix=data/{n} {fit_url}",
                    verbose=False,
                )

            # Verify data and repull it if necessary
            pull_failure = True
            while pull_failure:
                pull_failure = check_data(file_path)
                if pull_failure:
                    runcmd(
                        f"wget --directory-prefix=data/{n} {fit_url}",
                        verbose=False,
                    )


def check_data(file):
    """
    This function verifies that all the fits files have the necessary data, sometimes
    they don't pull everything for some reason
    """
    system_info = fits.info(file, output=False)
    with open(file, "rb") as f:
        # read header of first extension
        h = fits.getheader(f, ext=0, memmap=False)
    n_ext = h["N_EXT"]  # get the largest extension
    failure = len(system_info) != n_ext + 1
    if failure:
        # Number of tables doesn't match the number of tables that the header
        # says exists, delete file
        file.unlink()
    return failure


def plot_planet_positions(ax, system, cmap, ind):
    """
    Put each planet in the scatter plot

    Args:
        ax (matplotlib axis):
            figure axis
        system (System object):
            Object holding all of the system data
        cmap (matplotlib colormap):
            colormap used for the plot
        ind (integer):
            The index value to plot

    Returns:
        ax (matplotlib axis):
            figure axis

    """
    cvals = np.linspace(0, 1, len(system.planets))
    for i, planet in enumerate(system.planets):
        color = cmap(cvals[i])
        ms = planet_marker_size(
            planet.vectors.at[ind, "z"],
            planet.vectors.z,
            base_size=5 + planet.radius.to(u.R_earth).value,
            factor=0.2,
        )
        ax.scatter(
            planet.vectors.x[ind] * (u.m.to(u.AU)),
            planet.vectors.y[ind] * (u.m.to(u.AU)),
            label=f"Planet {i}",
            color=color,
            s=ms,
        )
        ax.set_xlim([-10, 10])
        ax.set_ylim([-10, 10])
        ax.set_xlabel("AU")
        ax.set_ylabel("AU")
    return ax


def planet_marker_size(z, all_z, base_size=5, factor=0.5):
    """
    Make the planet marker smaller when the planet is behind the star in its orbit
    """
    z_range = np.abs(max(all_z) - min(all_z))

    # Want being at max z to correspond to a factor of 1 and min z
    # to be a factor of negative 1
    scaled_z = 2 * (z - min(all_z)) / z_range - 1

    marker_size = base_size * (1 + factor * scaled_z)

    return marker_size
