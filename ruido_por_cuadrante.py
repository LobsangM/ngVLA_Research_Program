import numpy as np
from astropy.io import fits

archivo_input = "convolucion_0.05arcsec_rebin.fits"
archivo_output = "convolucion_0.05arcsec_rebin_4cuadrantes_ruido.fits"

# ruido de cada cuadrante

sigma_Q1 = 1e-13   # cuadrante superior izquierdo
sigma_Q2 = 3e-13   # cuadrante superior derecho
sigma_Q3 = 5e-13   # cuadrante inferior izquierdo
sigma_Q4 = 8e-13   # cuadrante inferior derecho


hdul = fits.open(archivo_input)
data = hdul[0].data
header = hdul[0].header

ny, nx = data.shape
mitad_x = nx // 2
mitad_y = ny // 2

# datos de la imagen
data_ruido = np.copy(data)

# distintos ruidos
ruido_Q1 = np.random.normal(
    loc=0.0,
    scale=sigma_Q1,
    size=(mitad_y, mitad_x)
)
data_ruido[:mitad_y, :mitad_x] += ruido_Q1


ruido_Q2 = np.random.normal(
    0.0,
    sigma_Q2,
    size=(mitad_y, nx - mitad_x)
)
data_ruido[:mitad_y, mitad_x:] += ruido_Q2


ruido_Q3 = np.random.normal(
    0.0,
    sigma_Q3,
    size=(ny - mitad_y, mitad_x)
)
data_ruido[mitad_y:, :mitad_x] += ruido_Q3


ruido_Q4 = np.random.normal(
    0.0,
    sigma_Q4,
    size=(ny - mitad_y, nx - mitad_x)
)
data_ruido[mitad_y:, mitad_x:] += ruido_Q4



fits.writeto(
    archivo_output,
    data_ruido,
    header,
    overwrite=True
)

print("Imagen con ruido gaussiano distinto por cuadrante.")
print(archivo_output)