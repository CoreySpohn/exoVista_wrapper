import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


def runcmd(cmd, verbose=False, *args, **kwargs):
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)
    pass


def get_exoVista_data():
    # Iterate over the different universes
    for n in tqdm(np.arange(1, 13), position=0, desc="Universe", leave=False):
        universe_url = f"https://ckan-files.emac.gsfc.nasa.gov/exovista/DEC21/{n}/target_database.csv"
        Path(f"{n}/").mkdir(exist_ok=True)
        if not Path(f"{n}/target_database.csv").exists():
            runcmd(f"wget --directory-prefix={n} {universe_url}", verbose=False)
        df = pd.read_csv(f"{n}/target_database.csv", low_memory=False)
        for i in tqdm(
            np.arange(1, df.shape[0]), position=1, desc="Planet", leave=False
        ):
            fit_url = df.at[i, "URL"]
            # fit_url = row['URL']
            if not Path(f"{n}/{fit_url.split('/')[-1]}").exists():
                runcmd(f"wget --directory-prefix={n} {fit_url}", verbose=False)
