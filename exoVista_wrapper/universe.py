"""
This is the class file to hold an entire exoVista universe in a single object
"""
import pickle
from pathlib import Path

from tqdm import tqdm

from exoVista_wrapper.system import System


class Universe:
    """
    Class for the whole exoVista universe
    """

    def __init__(self, path, cache=False):
        """
        Args:
            path (str or Path):
                Location of all the system files. Should be something like "data/1/"
        """
        if cache:
            cache_base = Path(".cache", path.split("/")[1])
            if not cache_base.exists():
                cache_base.mkdir(parents=True)
        self.path = path

        # Load all systems
        p = Path(path).glob("*.fits")
        system_files = [x for x in p if x.is_file]
        self.systems = []
        for system_file in tqdm(
            system_files, desc="Loading systems", position=0, leave=False
        ):
            if cache:
                cache_file = Path(cache_base, system_file.stem + ".p")
                if cache_file.exists():
                    with open(cache_file, "rb") as f:
                        system = pickle.load(f)
                else:
                    system = System(system_file)
                    if system is not None:
                        continue
                    with open(cache_file, "wb") as f:
                        pickle.dump(system, f)
                self.systems.append(system)
            else:
                system = System(system_file)
                if system is not None:
                    self.systems.append(system)

        # Get star ids
        self.ids = [system.star.exoVista_id for system in self.systems]
        self.HIP_ids = [system.star.HIP_id for system in self.systems]
