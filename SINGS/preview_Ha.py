# preview_Ha.py — genera PNGs individuales + collage de las imágenes Hα
# Ejecutar desde SINGS/: python preview_Ha.py

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ZScaleInterval, AsinhStretch, ImageNormalize

from utils import galaxias

OUTPUT_DIR = "imagenes_png"
os.makedirs(OUTPUT_DIR, exist_ok=True)

interval = ZScaleInterval()
stretch  = AsinhStretch()


def load_image(galaxia):
    path = os.path.join("imagenes_Ha", f"{galaxia}_HA_SUB_dr4.fits")
    with fits.open(path) as h:
        data = h[0].data.astype(float)
    # squeeeze por si tiene ejes extra (CASA)
    data = np.squeeze(data)
    return data


def save_individual(galaxia, data, norm):
    fig, ax = plt.subplots(figsize=(4, 4), dpi=150)
    ax.imshow(data, origin="lower", cmap="magma", norm=norm)
    ax.set_title(galaxia, fontsize=10, color="white", pad=4)
    ax.axis("off")
    fig.patch.set_facecolor("black")
    plt.tight_layout(pad=0.3)
    fig.savefig(
        os.path.join(OUTPUT_DIR, f"{galaxia}.png"),
        dpi=150,
        facecolor="black",
        bbox_inches="tight",
    )
    plt.close(fig)


# --- individuales ---
print(f"Generando {len(galaxias)} imágenes individuales...")
for i, galaxia in enumerate(galaxias, 1):
    data = load_image(galaxia)
    print(f"  [{i:2d}/{len(galaxias)}] {galaxia:14s}  min={data.min():.3e}  max={data.max():.3e}")
    vmin, vmax = interval.get_limits(data)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
    save_individual(galaxia, data, norm)

# --- collages por colormap ---
COLORMAPS = ["magma", "viridis"]

ncols = 9
nrows = int(np.ceil(len(galaxias) / ncols))

for cmap in COLORMAPS:
    print(f"\nGenerando collage ({cmap})...")
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 3), dpi=120)
    fig.patch.set_facecolor("black")

    for idx, ax in enumerate(axes.flat):
        if idx < len(galaxias):
            galaxia = galaxias[idx]
            data = load_image(galaxia)
            vmin, vmax = interval.get_limits(data)
            norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch)
            ax.imshow(data, origin="lower", cmap=cmap, norm=norm)
            ax.set_title(galaxia, fontsize=7, color="white", pad=2)
        ax.axis("off")

    plt.tight_layout(pad=0.5)
    fig.savefig(
        os.path.join(OUTPUT_DIR, f"collage_Ha_{cmap}.png"),
        dpi=120,
        facecolor="black",
        bbox_inches="tight",
    )
    plt.close(fig)

print(f"\nListo. PNGs en '{OUTPUT_DIR}/'  (+ collage_Ha_magma.png, collage_Ha_viridis.png)")
