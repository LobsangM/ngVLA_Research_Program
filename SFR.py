#File: SFR.py
#Date: 2025-02-20

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Archivo fits
fits_file = "ngc6946Ha_I_Ha_ksb2004.fits"  
hdul = fits.open(fits_file)
data = hdul[0].data
hdul.close()

# suma total de H_alpha
total_Ha_flux = np.sum(data)

# maximo
max_Ha_flux = np.max(data)

# SFR esperado
SFR_total = 500  # Msun/yr

# Factor de conversión
scaling_factor = SFR_total / total_Ha_flux

# Aplicar el reescalamiento lineal
SFR_map = data * scaling_factor

# Aplicar logarítmica 
SFR_map_log = np.log10(1 + SFR_map)  

# Verificar que la suma total sigue siendo 500 Msun/yr
SFR_check = np.sum(SFR_map)

print(f"Factor de escalamiento: {scaling_factor}")
print(f"Verificación SFR total: {SFR_check} Msun/yr")

# Mostrar imagen 
plt.imshow(SFR_map_log, origin="lower", cmap="inferno")
plt.colorbar(label="log10(1 + SFR) [Msun/yr per pixel]")
plt.title("Mapa de Formación Estelar (Escala Logarítmica)")
plt.show()
