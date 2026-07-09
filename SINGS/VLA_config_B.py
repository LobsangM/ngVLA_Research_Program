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
def add_correlated_noise(fits_file, beam_arcsec, noise_uJy):

    print("\nAñadiendo ruido a:", fits_file)

    hdul = fits.open(fits_file)

    data = hdul[0].data
    header = hdul[0].header

    # conversion de ruido
    sigma_beam = noise_uJy * 1e-6

    # escala de pixel
    pixel_deg = abs(header["CDELT1"])
    pixel_arcsec = pixel_deg * 3600.0

    print("Beam:", beam_arcsec, "arcsec")
    print("Noise:", noise_uJy, "uJy/beam")
    print("Pixel:", pixel_arcsec, "arcsec")
    print("Sigma beam:", sigma_beam, "Jy")

    # imsmooth conserva el flujo total (Jy/pixel), diluyendo el pico por
    # pixel en proporción al número de pixeles dentro del beam. Hay que
    # reescalar a Jy/beam (multiplicando por el área del beam en pixeles)
    # antes de sumar el ruido, que sí está definido en Jy/beam.
    beam_area_pix = (np.pi / (4.0 * np.log(2.0))) * (beam_arcsec / pixel_arcsec) ** 2
    data = data * beam_area_pix

    print("Pixeles/beam:", beam_area_pix)

    # Ruido blanco con std = sigma_beam. Al convolucionar con un kernel
    # normalizado en L2 (ver abajo), el ruido correlacionado resultante
    # conserva std ~ sigma_beam por pixel, consistente con Jy/beam.
    white_noise = np.random.normal(
        loc=0.0,
        scale=sigma_beam,
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

    # Normalización L2 (no L1): preserva el std del ruido tras la
    # convolución. kernel /= kernel.sum() conserva flujo (correcto para
    # suavizar señal) pero diluye la varianza del ruido.
    kernel /= np.sqrt((kernel**2).sum())

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
    header["BUNIT"] = "Jy/beam"

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

    fig = plt.figure(figsize=(8,8))

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

    # figura limpia
    f.axis_labels.hide()
    f.tick_labels.hide()

    # guardar figura
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


# config VLA B
beam      = 0.961
noise_uJy = 0.56
beam_str  = "%sarcsec" % beam

for entry in REDSHIFT_SFR_TABLE:
    z     = entry["redshift"]
    sfr   = entry["sfr"]
    label = entry["label"]

    output_dir = os.path.join("VLA_config_B", label)
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*55}")
    print(f"  VLA-B  z={z}  SFR={sfr} M☉/yr  [{label}]")
    print(f"{'='*55}")

    for galaxia in galaxias:
        print(f"\n  → {galaxia}")

        input_fits   = os.path.join("resultados", label, galaxia,
                                    f"mapa_z{z}_30kpc_SFR{sfr}.fits")
        casa_image   = "temp_vla_b.image"
        smooth_image = "temp_vla_b_smooth.image"
        vla_fits     = os.path.join(output_dir, f"VLA_B_{galaxia}.fits")

        importfits(fitsimage=input_fits, imagename=casa_image, overwrite=True)
        imsmooth(imagename=casa_image, kernel="gauss",
                 major=beam_str, minor=beam_str, pa="0deg",
                 outfile=smooth_image, overwrite=True)
        exportfits(imagename=smooth_image, fitsimage=vla_fits, overwrite=True)

        final_fits = add_correlated_noise(vla_fits, beam, noise_uJy)
        plot_with_beam(final_fits)

print("\nSimulaciones VLA B terminadas.")