import os
import numpy as np
from astropy.io import fits
from scipy.signal import fftconvolve
import copy
import matplotlib.pyplot as plt
import aplpy
import sys

sys.path.append("/Users/lobsangmendez/Documents/ResearchP_ngVLA")

from sfr_model_CASA import config, SFRRadioModel


#Ruido correlacionado
def add_correlated_noise(fits_file, beam_arcsec, noise_uJy):

    print("\nAñadiendo ruido a:", fits_file)

    hdul = fits.open(fits_file)

    data = hdul[0].data
    header = hdul[0].header

    #conversión de ruido
    sigma_beam = noise_uJy * 1e-6

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
    print("Noise:", noise_uJy, "uJy/beam")
    print("Pixel:", pixel_arcsec, "arcsec")
    print("Pixels/beam:", N_pix_beam)
    print("Sigma pixel:", sigma_pixel, "Jy")

   #Ruido blanco
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

    #beam
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


# plot beam
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

    plt.show()


#config VLA A
beam = 0.293
noise_uJy = 0.56

cfg = copy.deepcopy(config)

cfg["output_folder"] = "VLA_A_mapas"

cfg["radio"]["beam_major"] = f"{beam}arcsec"
cfg["radio"]["beam_minor"] = f"{beam}arcsec"


# Simulación
model = SFRRadioModel(cfg)

# luminosidad -> flujo
model.luminosity_to_flux()

# guardar mapa base
model.add_wcs_and_save()

# convolución con CASA
model.casa_convolution()


#añadir ruido
fits_file = os.path.join(
    "VLA_A_mapas",
    f"convolucion_{beam}arcsec.fits"
)

final_fits = add_correlated_noise(
    fits_file,
    beam,
    noise_uJy
)


# imagen con beam
plot_with_beam(final_fits)

print("\nSimulación VLA A terminada.")