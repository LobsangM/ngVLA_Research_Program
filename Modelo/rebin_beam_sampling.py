# Imrebin 
# Pixel objetivo ≈ 0.01 arcsec


# Imagen original (FITS)
fits_input = "resultados/mapa_base/convolucion_0.05arcsec.fits"

# Imagen en formato CASA
casa_image = "convolucion_0.05arcsec.image"

# Imagen rebin
rebin_image = "convolucion_0.05arcsec_rebin.image"

# Imagen final FITS
fits_output = "resultados/mapa_base/convolucion_0.05arcsec_rebin.fits"


print("Factor rebin ≈ 23")

#----

importfits(
    fitsimage=fits_input,
    imagename=casa_image,
    overwrite=True
)



#-----

imrebin(
    imagename=casa_image,
    outfile=rebin_image,
    factor=[23,23]
)



# ----

exportfits(
    imagename=rebin_image,
    fitsimage=fits_output,
    overwrite=True
)


print(fits_output)
