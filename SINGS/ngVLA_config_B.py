import os
import sys
sys.path.insert(0, '.')

import numpy as np
from astropy.io import fits
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
import aplpy
from casatasks import importfits, imsmooth, exportfits

from utils import galaxias, REDSHIFT_SFR_TABLE


# Ruido correlacionado
def add_correlated_noise(fits_file, beam_arcsec, noise_nJy):

    print("\nAñadiendo ruido a:", fits_file)

    hdul = fits.open(fits_file)

    data = hdul[0].data
    header = hdul[0].header

    # conversión de ruido
    sigma_beam = noise_nJy * 1e-9

    # escala de pixel
    pixel_deg = abs(header["CDELT1"])
    pixel_arcsec = pixel_deg * 3600.0

    # áreas
    beam_area = 1.1331 * beam_arcsec**2
    pixel_area = pixel_arcsec**2

    # pixeles por beam
    N_pix_beam = beam_area / pixel_area

    # ruido por pixel
    sigma_pixel = sigma_beam / np.sqrt(N_pix_beam)

    print("Beam:", beam_arcsec, "arcsec")
    print("Noise:", noise_nJy, "nJy/beam")
    print("Pixel:", pixel_arcsec, "arcsec")
    print("Pixels/beam:", N_pix_beam)
    print("Sigma pixel:", sigma_pixel, "Jy")

    # Ruido blanco
    white_noise = np.random.normal(
        loc=0.0,
        scale=sigma_pixel,
        size=data.shape
    )

    # Relación FWHM
    fwhm_pix = beam_arcsec / pixel_arcsec

    sigma_kernel_pix = fwhm_pix / (
        2.0 * np.sqrt(2.0 * np.log(2.0))
    )

    print("FWHM [pix]:", fwhm_pix)
    print("Sigma kernel [pix]:", sigma_kernel_pix)

    kernel_size = int(6 * sigma_kernel_pix)

    t = np.linspace(
        -kernel_size,
        kernel_size,
        2*kernel_size + 1
    )

    x, y = np.meshgrid(t, t)

    kernel = np.exp(
        -(x**2 + y**2) /
        (2 * sigma_kernel_pix**2)
    )

    kernel /= kernel.sum()

    # Ruido correlacionado
    correlated_noise = fftconvolve(
        white_noise,
        kernel,
        mode="same"
    )

    # imagen final
    final_data = data + correlated_noise

    # beam
    beam_deg = beam_arcsec / 3600.0

    header["BMAJ"] = beam_deg
    header["BMIN"] = beam_deg
    header["BPA"] = 0.0

    header.comments["BMAJ"] = "Beam major axis [deg]"
    header.comments["BMIN"] = "Beam minor axis [deg]"
    header.comments["BPA"] = "Beam position angle [deg]"

    # guardar fits
    output = fits_file.replace(
        ".fits",
        "_noise.fits"
    )

    fits.writeto(
        output,
        final_data,
        header,
        overwrite=True
    )

    print("Archivo final:", output)

    return output


# dibujar beam
def plot_with_beam(fits_file):

    print("\nGenerando figura con beam...")

    fig = plt.figure(figsize=(8, 8))

    f = aplpy.FITSFigure(
        fits_file,
        figure=fig
    )

    # mostrar imagen
    f.show_colorscale(
        cmap="inferno"
    )

    # añadir beam
    f.add_beam()

    f.beam.set_corner("bottom left")

    f.beam.set_edgecolor("white")
    f.beam.set_facecolor("none")

    f.beam.set_linewidth(2)

    f.axis_labels.hide()
    f.tick_labels.hide()

    output_png = fits_file.replace(
        ".fits",
        ".png"
    )

    plt.savefig(
        output_png,
        dpi=300,
        bbox_inches="tight"
    )

    print("Figura guardada:", output_png)

    #plt.show()


# config ngVLA B
beam      = 0.961
noise_nJy = 32.68
beam_str  = "%sarcsec" % beam

for entry in REDSHIFT_SFR_TABLE:
    z     = entry["redshift"]
    sfr   = entry["sfr"]
    label = entry["label"]

    output_dir = os.path.join("ngVLA_config_B", label)
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*55}")
    print(f"  ngVLA-B  z={z}  SFR={sfr} M☉/yr  [{label}]")
    print(f"{'='*55}")

    for galaxia in galaxias:
        print(f"\n  → {galaxia}")

        input_fits   = os.path.join("resultados", label, galaxia,
                                    f"mapa_z{z}_30kpc_SFR{sfr}.fits")
        casa_image   = "temp_ngvla_b.image"
        smooth_image = "temp_ngvla_b_smooth.image"
        ngvla_fits   = os.path.join(output_dir, f"ngVLA_B_{galaxia}.fits")

        importfits(fitsimage=input_fits, imagename=casa_image, overwrite=True)
        imsmooth(imagename=casa_image, kernel="gauss",
                 major=beam_str, minor=beam_str, pa="0deg",
                 outfile=smooth_image, overwrite=True)
        exportfits(imagename=smooth_image, fitsimage=ngvla_fits, overwrite=True)

        final_fits = add_correlated_noise(ngvla_fits, beam, noise_nJy)
        plot_with_beam(final_fits)

print("\nSimulaciones ngVLA B terminadas.")
