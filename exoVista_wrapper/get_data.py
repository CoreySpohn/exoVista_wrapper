import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
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
    # Iterate over the different universes
    for n in tqdm(universes, position=0, desc="Universe", leave=False):
        universe_url = f"https://ckan-files.emac.gsfc.nasa.gov/exovista/DEC21/\
                {n}/target_database.csv"
        Path(f"exoVista_data/{n}/").mkdir(parents=True)
        if not Path(f"exoVista_data/{n}/target_database.csv").exists():
            runcmd(f"wget --directory-prefix={n} {universe_url}", verbose=False)
        df = pd.read_csv(f"exoVista_data/{n}/target_database.csv", low_memory=False)
        for i in tqdm(
            np.arange(1, df.shape[0]), position=1, desc="Planet", leave=False
        ):
            fit_url = df.at[i, "URL"]
            # fit_url = row['URL']
            if not Path(f"exoVista_data/{n}/{fit_url.split('/')[-1]}").exists():
                runcmd(
                    f"wget --directory-prefix=exoVista_data/{n} {fit_url}",
                    verbose=False,
                )
