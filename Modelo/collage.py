import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from matplotlib.colors import Normalize
from astropy.io import fits
from astropy.visualization import ZScaleInterval

BASE_DIR = Path(__file__).parent
sys.path.insert(0, str(BASE_DIR))

COLLAGE_DIR = BASE_DIR / "collage"

