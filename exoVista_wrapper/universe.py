"""
This is the class file to hold an entire exoVista universe in a single object
"""
from pathlib import Path

from tqdm import tqdm

from exoVista_wrapper.system import System


class Universe:
    """
    Class for the whole exoVista universe
    """

    def __init__(self, path):
        """
        Args:
            path (str or Path):
                Location of all the system files. Should be something like "data/1/"
        """
        # Load all systems
        p = Path(path).glob("*.fits")
        system_files = [x for x in p if x.is_file]
        self.systems = []
        for system_file in tqdm(
            system_files, desc="Loading systems", position=0, leave=False
        ):
            system = System(system_file)
            if system is not None:
                self.systems.append(system)
