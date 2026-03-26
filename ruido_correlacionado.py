import numpy as np
from astropy.io import fits
from scipy.signal import fftconvolve

archivo_input = "convolucion_0.05arcsec_rebin.fits"
archivo_output = "convolucion_0.05arcsec_rebin_ruido_correlacionado.fits"

# config de ruido
sigma_ruido = 8e-12      # en Jy
sigma_kernel = 3         # correlación espacial en pixeles


print("Sigma del ruido =", sigma_ruido, "Jy")
print("Sigma del kernel 2D =", sigma_kernel, "pix")



hdul = fits.open(archivo_input)
data = hdul[0].data
header = hdul[0].header

ny, nx = data.shape

# ruido blanco
ruido_blanco = np.random.normal(
    loc=0.0,
    scale=sigma_ruido,
    size=data.shape
)


# kernel
t = np.linspace(-3*sigma_kernel, 3*sigma_kernel, 6*sigma_kernel + 1)
x, y = np.meshgrid(t, t)

kernel = np.exp(-(x**2 + y**2) / (2 * sigma_kernel**2))
kernel /= kernel.sum()   # normalizar energía


# convolucionar
ruido_correlacionado = fftconvolve(ruido_blanco, kernel, mode="same")



data_ruido = data + ruido_correlacionado



fits.writeto(
    archivo_output,
    data_ruido,
    header,
    overwrite=True
)

print("Ruido correlacionado generado y guardado en:")
print(archivo_output)