"""Top-level package for autoadsorbate."""

__author__ = """Fakoe Edvin"""
__email__ = "edvin.fako@basf.com"
__version__ = "0.2.0"

from os.path import dirname, basename, isfile, join
import glob
from .autoadsorbate import *

modules = glob.glob(join(dirname(__file__), "*.py"))
__all__ = [
    basename(f)[:-3] for f in modules if isfile(f) and not f.endswith("__init__.py")
]
