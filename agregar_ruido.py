import numpy as np
from astropy.io import fits


archivo_input = "convolucion_0.05arcsec_rebin.fits"
archivo_output = "convolucion_0.05arcsec_rebin_ruido.fits"

# ruido
sigma_ruido = 1e-36   # 

print("Agregando ruido gaussiano")
print("sigma =", sigma_ruido, "Jy")


hdul = fits.open(archivo_input)

data = hdul[0].data
header = hdul[0].header

# generar ruido
ruido = np.random.normal(
    loc=0.0,
    scale=sigma_ruido,
    size=data.shape
)


data_ruido = data + ruido

fits.writeto(
    archivo_output,
    data_ruido,
    header,
    overwrite=True
)

hdul.close()

print("Imagen con ruido guardada:")
print(archivo_output)
