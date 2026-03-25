from astropy.io import fits

# ===============================
# Parámetros
# ===============================

archivo = "flux_map_z1_SFR100.fits"

tamano_fisico_kpc = 10      # 10 kpc totales
npix = 2900
kpc_per_arcsec = 8.12       # para z = 1 , sacado con la calculadora

# ===============================
# Cálculos
# ===============================

kpc_per_pix = tamano_fisico_kpc / npix
arcsec_per_pix = kpc_per_pix / kpc_per_arcsec
deg_per_pix = arcsec_per_pix / 3600

print("kpc/pix:", kpc_per_pix)
print("arcsec/pix:", arcsec_per_pix)
print("deg/pix:", deg_per_pix)

# ===============================
# Modificación del header
# ===============================

hdul = fits.open(archivo)
data = hdul[0].data
header = hdul[0].header.copy()

# WCS
header['CTYPE1'] = 'RA---SIN'
header['CTYPE2'] = 'DEC--SIN'
header['CUNIT1'] = 'deg'
header['CUNIT2'] = 'deg'

header['CRPIX1'] = npix / 2
header['CRPIX2'] = npix / 2

header['CRVAL1'] = 150.0
header['CRVAL2'] = 2.0

header['CDELT1'] = -deg_per_pix
header['CDELT2'] = deg_per_pix

header['COMMENT'] = "Simulacion z=1, 10 kpc x 10 kpc"

fits.writeto("mapa_z1_10kpc_SFR100.fits", data, header, overwrite=True)
hdul.close()
