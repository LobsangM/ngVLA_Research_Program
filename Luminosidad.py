import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

#definimos cosmo
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

# Constantes a utilizar
alpha = 0.7
redshifts = [1, 2, 3, 4, 5]


# Cargar imagen
sfr_fits = fits.open("ngc6946Ha_I_Ha_ksb2004.fits")  
sfr_data = sfr_fits[0].data  


scaling_factor =  500 / np.sum(sfr_data)    #====================================
SFR_map = scaling_factor * sfr_data         #    Esto es la imagen original
SFR_map_log = np.log10(1 + SFR_map)         #====================================



# Conversión a luminosidad en 1.4 GHz y extraemos datos
lum_1_4GHz = sfr_data / 4.87e-29
lum_data = lum_1_4GHz[0].data 

# Valor de la luminosidad Total:
lum_total = np.sum(lum_1_4GHz)
print(f"Luminosidad total: {lum_total:.2e} erg/s/Hz")

# ====================================================================================================

# Mostrar imagen orginal en SFR

plt.imshow(SFR_map_log, origin="lower", cmap="inferno", norm=LogNorm(vmin=4e-5 , vmax=4e-3))
plt.colorbar(label="log10(1 + SFR) [Msun/yr per pixel]")
plt.title("Mapa de Formación Estelar (Escala Logarítmica)")
plt.show()

# Imagen convertida a Luminsidad

plt.imshow(lum_1_4GHz, origin="lower", cmap="inferno", norm=plt.Normalize(vmin=np.percentile(lum_1_4GHz, 5), vmax=np.percentile(lum_1_4GHz, 95)))
plt.colorbar(label="Luminosidad [erg/s/Hz]")
plt.title("Mapa de Luminosidad a 1.4 GHz")
plt.show()

# ====================================================================================================


# Convertir de luminosidad a densidad de flujo a 10 GHz
def luminosity_to_flux(lum_1_4GHz, z, alpha):

    D_L = cosmo.luminosity_distance(z).to('cm').value

    factor = (4 * np.pi * D_L * z**2) / ((1+z)**(1-alpha) * (1.4/10)**(-alpha))
    S_10GHz = lum_1_4GHz / factor

    return S_10GHz



flujo_total = np.sum(luminosity_to_flux(lum_1_4GHz, 1, alpha))
print(f"Densidad de flujo total a 10 GHz: {flujo_total} erg/s/cm^2/Hz")


# #Le aplicamos el reescalamiento logaritmico a los nuevos datos

factor_flux = 500 / np.sum(lum_data)
S_10GHz_map = factor_flux * luminosity_to_flux(lum_1_4GHz, 1, alpha)
S_10GHz_log = np.log10(1 + S_10GHz_map)

#Grafica con diferentes redshifts
plt.imshow(S_10GHz_log, origin='lower', cmap='inferno', norm=LogNorm(vmin=1e-5, vmax=1e-3))
plt.colorbar(label="Densidad de flujo [erg/s/cm^2/Hz]")
plt.title("Densidad de Flujo a 10 GHz para z=4")
plt.show()







# # Convertir para diferentes valores de z
# for i, z in enumerate(redshifts):
#     S_10GHz = luminosity_to_flux(luminosidad, z, alpha)
    
#     ax = axes[i]
#     im = ax.imshow(S_10GHz, origin='lower', cmap='inferno')
#     ax.set_title(f"z = {z}")
#     plt.colorbar(im, ax=ax, orientation='vertical')

#     print(f"Archivo flux_10GHz_z{z}.fits guardado correctamente.")

# plt.suptitle("Densidad de Flujo a 10 GHz para Diferentes z")
# plt.tight_layout()
# plt.show()