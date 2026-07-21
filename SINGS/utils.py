# file name: utils.py
# date: 13-05-2026

import tomli as tomllib
from pathlib import Path

BASE_DIR = Path(__file__).parent
TOML_PATH = BASE_DIR / "tomls" / "Ha_SUB.toml"

with open(TOML_PATH, "rb") as f:
    toml_data = tomllib.load(f)

sufijo   = toml_data["galaxias"]["archivo_sufijo"]

# Solo galaxias con FITS disponible en imagenes_Ha/
galaxias = [
    g for g in toml_data["galaxias"]["galaxias"]
    if (BASE_DIR / "imagenes_Ha" / f"{g}{sufijo}").exists()
]

# Tabla de redshift y SFR típico de la Secuencia Principal
# basada en Leslie et al. (2020), para M* = 1e10 M_sun
REDSHIFT_SFR_TABLE = [
    {"redshift": 0.4, "sfr": 1.5, "label": "z0.4"},
    {"redshift": 1.0, "sfr":   7, "label": "z1.0"},
    {"redshift": 2.0, "sfr":  40, "label": "z2.0"},
    {"redshift": 3.0, "sfr": 100, "label": "z3.0"},
    {"redshift": 5.0, "sfr": 300, "label": "z5.0"},
]

# Estrellas de campo mal sustraídas en el FITS Hα crudo de SINGS DR4 (residuo
# dipolar: pico positivo muy intenso pegado a un pozo negativo, en imágenes
# de continuo restado). No son formación estelar real: si no se enmascaran,
# dominan la normalización de SFR y aparecen como una fuente compacta falsa
# tras la convolución del beam, más brillante que la galaxia real. Ver
# README.md "Máscaras manuales de estrellas mal sustraídas". Coordenadas en
# el sistema de la imagen ORIGINAL (sin recortar), (x=columna, y=fila, radio
# en píxeles nativos), encontradas por inspección del pico de flujo dentro
# del recorte y su perfil radial.
STAR_MASKS = {
    "ngc3627": [(1403, 1531, 15)],
    "ngc4594": [(1724, 852, 35)],
    "ngc5055": [(1601, 1218, 45)],
}


def get_config(galaxia, redshift=1.0, sfr=7, label="z1.0"):
    fits_file = str(BASE_DIR / "imagenes_Ha" / f"{galaxia}{sufijo}")
    return {
        "cosmology": {
            "H0": 70,
            "Om0": 0.3,
            "Tcmb0": 2.725
        },
        "image": {
            "fits_file": fits_file,
            "crop_size": 30,
            "crop_fraction": 0.7,
            "star_masks": STAR_MASKS.get(galaxia, [])
        },
        "sfr": {
            "total_sfr": sfr
        },
        "radio": {
            "alpha": 0.7,
            "redshift": redshift,
            "obs_frequency": 10,
            "rest_frequency": 1.4,
            "beam_major": "0.05arcsec",
            "beam_minor": "0.05arcsec",
            "beam_pa": "0deg"
        },
        "output_folder": str(BASE_DIR / "resultados" / label / galaxia),
        "galaxia": galaxia
    }