
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
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval

BASE_DIR = Path(__file__).parent
sys.path.insert(0, str(BASE_DIR))
from utils import galaxias, REDSHIFT_SFR_TABLE

COLLAGE_DIR = BASE_DIR / "collage"

# Configuraciones
CONFIGS = [
    {
        "name":   "ngVLA-A",
        "folder": "ngVLA_config_A",
        "prefix": "ngVLA_A",
        "beam":   0.293,
        "noise":  "28.4 nJy/beam",
        "cmap":   "plasma",
        "color":  "#60a5fa",   # azul
    },
    {
        "name":   "ngVLA-B",
        "folder": "ngVLA_config_B",
        "prefix": "ngVLA_B",
        "beam":   0.961,
        "noise":  "32.7 nJy/beam",
        "cmap":   "plasma",
        "color":  "#34d399",   # verde
    },
    {
        "name":   "VLA-A",
        "folder": "VLA_config_A",
        "prefix": "VLA_A",
        "beam":   0.293,
        "noise":  "0.56 µJy/beam",
        "cmap":   "hot",
        "color":  "#fb923c",   # naranja
    },
    {
        "name":   "VLA-B",
        "folder": "VLA_config_B",
        "prefix": "VLA_B",
        "beam":   0.961,
        "noise":  "0.56 µJy/beam",
        "cmap":   "hot",
        "color":  "#f472b6",   # rosa
    },
]

# Colores de cada fila
Z_COLORS = ["#93c5fd", "#6ee7b7", "#fcd34d", "#f9a8d4", "#d8b4fe"]

BG_COLOR   = "#0d0d0d"
PANEL_BG   = "#111111"
GRAY_TEXT  = "#888888"
WHITE      = "#ffffff"


# FITS

def read_fits(path):

    p = Path(path)
    if not p.exists():
        return None, None
    try:
        with fits.open(p) as hdul:
            data   = np.squeeze(hdul[0].data).astype(float)
            header = hdul[0].header
        return data, header
    except Exception as exc:
        print(f"    [!] {p.name}: {exc}")
        return None, None


def auto_norm(data, mode="percentile"):
    
    if data is None:
        return Normalize(0, 1)
    v = data[np.isfinite(data)]
    if v.size == 0:
        return Normalize(0, 1)
    if mode == "zscale":
        lo, hi = ZScaleInterval(contrast=0.25).get_limits(v)
        return Normalize(max(0.0, float(lo)), float(hi))
    else:  # percentile
        hi = float(np.nanpercentile(data, 99.5))
        return Normalize(0.0, max(hi, 1e-30))


def beam_in_pixels(header):
    
    if not header:
        return None
    bmaj = header.get("BMAJ")               # grados
    cdelt = abs(header.get("CDELT2", 0.0))  # grados/pixel
    if bmaj and cdelt > 0:
        return float(bmaj) / float(cdelt)
    return None


def sfr_tag(sfr):
    
    return str(int(sfr)) if sfr == int(sfr) else str(sfr)


# Collage

def draw_panel(ax, data, norm, cmap, title="", unit="",
               title_color=WHITE, beam_px=None, note=None, aspect="auto",
               contour=False):

    ax.set_facecolor(PANEL_BG)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor("#2a2a2a")

    if data is not None:
        d = np.where(np.isfinite(data), data, 0.0)
        im = ax.imshow(d, origin="lower", cmap=cmap, norm=norm, aspect=aspect)

        # Colorbar
        cb = plt.colorbar(im, ax=ax, pad=0.01, fraction=0.046)
        cb.ax.tick_params(labelsize=5, colors=GRAY_TEXT)
        cb.ax.yaxis.label.set_color(GRAY_TEXT)
        if unit:
            cb.set_label(unit, fontsize=5, color=GRAY_TEXT)

        # Contorno de detección: RMS por sigma-clipping (descarta la propia
        # fuente para no sesgar el fondo hacia arriba) y contornos a
        # 3/5/10/20 sigma, mismo criterio y color que plot_with_beam() en
        # los scripts de telescopio. Sin reproyección WCS acá (imshow plano
        # en pixeles), así que ax.contour() sobre el mismo array basta -- se
        # alinea exacto con la imagen sin necesidad de aplpy. Ver README.md
        # "Contornos de detección (3 sigma)".
        if contour:
            _, _, rms = sigma_clipped_stats(d, sigma=3.0, maxiters=5)
            levels = [lvl for lvl in (3*rms, 5*rms, 10*rms, 20*rms) if d.max() > lvl]
            if levels:
                ax.contour(d, levels=levels, colors="#00ffcc", linewidths=0.7)

        # Elipse del beam (esquina inferior izquierda)
        if beam_px and beam_px > 0:
            ny, nx = d.shape
            margin = min(nx, ny) * 0.09
            e = Ellipse(
                (margin, margin),
                width=beam_px, height=beam_px,
                edgecolor=WHITE, facecolor=WHITE,
                alpha=0.75, zorder=5,
            )
            ax.add_patch(e)
    else:
        ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                transform=ax.transAxes, color="#444", fontsize=9)

    if title:
        ax.set_title(title, fontsize=8, fontweight="bold",
                     color=title_color, pad=3)
    if note:
        ax.text(0.5, -0.025, note, ha="center", va="top",
                transform=ax.transAxes, fontsize=5.5, color=GRAY_TEXT)


# titulo de columna
def draw_header_panel(ax, text, color=WHITE):
    
    ax.set_facecolor(BG_COLOR)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0.5, 0.5, text, ha="center", va="center",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold",
            color=color, multialignment="center")


# Collage por galaxia

def make_collage(galaxy):
    print(f"  {galaxy} ...", end="", flush=True)

    nz = len(REDSHIFT_SFR_TABLE)   # 5 redshifts
    nc = 1 + len(CONFIGS)           # 1 SFR + 4 telescopios = 5 columnas datos

    
    ha_data, _ = read_fits(
        BASE_DIR / "imagenes_Ha" / f"{galaxy}_HA_SUB_dr4.fits"
    )
    ha_norm = auto_norm(ha_data, mode="zscale")

    # Figura
    fw = 3.8 + nc * 3.1 + 0.3
    fh = 0.45 + nz * 3.2 + 0.7    # 0.7 para el suptitle
    fig = plt.figure(figsize=(fw, fh), facecolor=BG_COLOR)

    fig.suptitle(
        f"{galaxy.upper()}  —  Simulación Radio  ngVLA / VLA",
        color=WHITE, fontsize=13, fontweight="bold", y=0.998,
    )

    
    outer = gridspec.GridSpec(
        1, 2, figure=fig,
        width_ratios=[1.25, nc],
        left=0.01, right=0.99,
        top=0.965, bottom=0.02,
        wspace=0.04,
    )


    ax_ha = fig.add_subplot(outer[0, 0])

    
    inner = gridspec.GridSpecFromSubplotSpec(
        nz + 1, nc,
        subplot_spec=outer[0, 1],
        height_ratios=[0.14] + [1] * nz,
        hspace=0.30,
        wspace=0.18,
    )

    
    draw_panel(
        ax_ha, ha_data, ha_norm, "afmhot",
        title="Hα  original",
        unit="counts",
        note="SINGS DR4  ·  imagen de entrada",
        aspect="equal",
    )

    # Cabeceras
    headers = (
        ["Mapa radio\n(30 kpc crop)"]
        + [f"{c['name']}\n{c['beam']}\"  ·  {c['noise']}" for c in CONFIGS]
    )
    col_colors = [WHITE] + [c["color"] for c in CONFIGS]

    for j, (hdr_txt, hdr_col) in enumerate(zip(headers, col_colors)):
        draw_header_panel(fig.add_subplot(inner[0, j]), hdr_txt, color=hdr_col)

    # Filas de datos
    for i, rz in enumerate(REDSHIFT_SFR_TABLE):
        z, sfr, label = rz["redshift"], rz["sfr"], rz["label"]
        zc = Z_COLORS[i]          # color de fila
        row = i + 1               # fila real en inner (fila 0 = cabecera)

        # Mapa SFR/radio
        sfr_file = BASE_DIR / "resultados" / label / galaxy / \
                   f"mapa_{label}_30kpc_SFR{sfr_tag(sfr)}.fits"
        sfr_data, _ = read_fits(sfr_file)
        sfr_norm = auto_norm(sfr_data)

        draw_panel(
            fig.add_subplot(inner[row, 0]),
            sfr_data, sfr_norm, "viridis",
            title=f"z = {z}  ·  SFR = {sfr} M☉/yr",
            unit="Jy",
            title_color=zc,
        )

        # 4 configuraciones
        for j, cfg in enumerate(CONFIGS):
            noise_file = (
                BASE_DIR / cfg["folder"] / label
                / f"{cfg['prefix']}_{galaxy}_noise.fits"
            )
            t_data, t_hdr = read_fits(noise_file)
            t_norm = auto_norm(t_data)
            bpx    = beam_in_pixels(t_hdr)

            draw_panel(
                fig.add_subplot(inner[row, j + 1]),
                t_data, t_norm, cfg["cmap"],
                title=f"z = {z}",
                unit="Jy/beam",
                title_color=zc,
                beam_px=bpx,
                contour=True,
            )

    
    out_path = COLLAGE_DIR / f"{galaxy}.png"
    plt.savefig(out_path, dpi=140, bbox_inches="tight",
                facecolor=BG_COLOR, edgecolor="none")
    plt.close(fig)
    print(f" ok → {out_path.name}")


# Main

def main():
    COLLAGE_DIR.mkdir(exist_ok=True)
    print(f"Generando {len(galaxias)} collages en {COLLAGE_DIR}\n")
    for galaxy in galaxias:
        make_collage(galaxy)
    print(f"\nFinalizado. {len(galaxias)} collages guardados en:\n  {COLLAGE_DIR}/")


if __name__ == "__main__":
    main()
