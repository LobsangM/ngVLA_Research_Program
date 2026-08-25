import os
import sys
sys.path.insert(0, '.')

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
import aplpy
from casatasks import importfits, imsmooth, exportfits

from utils import galaxias, REDSHIFT_SFR_TABLE, NOISE_CONFIG


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

    print("Beam:", beam_arcsec, "arcsec")
    print("Noise:", noise_nJy, "nJy/beam")
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

    # Contorno de detección: RMS estimado con sigma-clipping (descarta la
    # fuente para no sesgar el fondo hacia arriba) y contornos a 3/5/10/20
    # sigma. Distingue emisión real de picos de ruido correlacionado que se
    # ven "calientes" por el mapa de color sin superar el umbral estadístico.
    # Ver README.md "Contornos de detección (3 sigma)".
    data = fits.getdata(fits_file)
    _, _, rms = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    # show_contour falla si un nivel no tiene ningún pixel por encima (el
    # WCS transform de matplotlib no admite un contorno vacío) - se filtran
    # los niveles que la imagen realmente alcanza antes de dibujar.
    levels = [lvl for lvl in (3*rms, 5*rms, 10*rms, 20*rms) if data.max() > lvl]
    if levels:
        f.show_contour(
            fits_file,
            levels=levels,
            colors="#00ffcc",
            linewidths=1.2
        )
    print("RMS estimado (sigma-clipped):", rms, "Jy/beam")

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


# config ngVLA A (ngECT Subarray: Mid+Spiral+OuterCore, synthesized beam @1.4GHz)
# Beam angosto/marginal, como config A en el resto del README. Core puro
# descartado: su beam (21.55") superaba el FOV completo de la imagen
# recortada (crop_size=30 en utils.py) a z>=0.4, colgando imsmooth por un
# padding de FFT gigantesco y sin aportar señal físicamente interpretable.
beam      = 0.99
noise_nJy = NOISE_CONFIG["ngVLA_A"]
beam_str  = "%sarcsec" % beam

for entry in REDSHIFT_SFR_TABLE:
    z     = entry["redshift"]
    sfr   = entry["sfr"]
    label = entry["label"]

    output_dir = os.path.join("Resultados", "ngVLA_config_A", label)
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*55}")
    print(f"  ngVLA-A  z={z}  SFR={sfr} M☉/yr  [{label}]")
    print(f"{'='*55}")

    for galaxia in galaxias:
        print(f"\n  → {galaxia}")

        input_fits   = os.path.join("resultados", label, galaxia,
                                    f"mapa_z{z}_30kpc_SFR{sfr}.fits")
        casa_image   = "temp_ngvla_a.image"
        smooth_image = "temp_ngvla_a_smooth.image"
        ngvla_fits   = os.path.join(output_dir, f"ngVLA_A_{galaxia}.fits")

        importfits(fitsimage=input_fits, imagename=casa_image, overwrite=True)
        imsmooth(imagename=casa_image, kernel="gauss",
                 major=beam_str, minor=beam_str, pa="0deg",
                 outfile=smooth_image, overwrite=True)
        exportfits(imagename=smooth_image, fitsimage=ngvla_fits, overwrite=True)

        final_fits = add_correlated_noise(ngvla_fits, beam, noise_nJy)
        plot_with_beam(final_fits)

print("\nSimulaciones ngVLA A terminadas.")
