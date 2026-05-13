# file name: utils.py
# date: 13-05-2026

import tomli as tomllib
from pathlib import Path

BASE_DIR = Path(__file__).parent
TOML_PATH = BASE_DIR / "tomls" / "Ha_SUB.toml"

with open(TOML_PATH, "rb") as f:
    toml_data = tomllib.load(f)

galaxias = toml_data["galaxias"]["galaxias"]
sufijo   = toml_data["galaxias"]["archivo_sufijo"]


def get_config(galaxia):
    fits_file = str(BASE_DIR / "imagenes_Ha" / f"{galaxia}{sufijo}")
    return {
        "cosmology": {
            "H0": 70,
            "Om0": 0.3,
            "Tcmb0": 2.725
        },
        "image": {
            "fits_file": fits_file,
            "crop_size": 30  # en kpc
        },
        "sfr": {
            "total_sfr": 100
        },
        "radio": {
            "alpha": 0.7,
            "redshift": 1,
            "obs_frequency": 10,
            "rest_frequency": 1.4,
            "beam_major": "0.05arcsec",
            "beam_minor": "0.05arcsec",
            "beam_pa": "0deg"
        },
        "output_folder": str(BASE_DIR / "resultados" / galaxia),
        "galaxia": galaxia
    }