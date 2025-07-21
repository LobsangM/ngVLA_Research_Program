# EN ESTE ARCHIVO SE MUESTRAN ALGUNOS CÓDIGOS QUE NO PUDIERON SER MOSTRADOS
# EN EL ARCHIVO aplpy_test.ipynb COMO LAS GRÁFICAS 

# File: aplpy_test.py
# Date: 2025-07-18
import aplpy
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

data = fits.getdata('ngc6946Ha_I_Ha_ksb2004.fits')
vmin = np.nanmin(data)
vmax = np.nanmax(data)



fig = aplpy.FITSFigure('ngc6946Ha_I_Ha_ksb2004.fits')

fig.show_grayscale()
#fig.hide_grayscale()
fig.show_colorscale(cmap='nipy_spectral', stretch='log' , vmin=0.01, vmax=vmax)
#fig.hide_colorscale()




hdul = fits.open('ngc6946Ha_I_Ha_ksb2004.fits')
wcs = WCS(hdul[0].header)

data = hdul[0].data
ny, nx = data.shape 
x_center = nx / 2
y_center = ny / 2

ra_dec = wcs.pixel_to_world(x_center, y_center)
# print(ra_dec.ra.deg, ra_dec.dec.deg)


fig.recenter(ra_dec.ra.deg, ra_dec.dec.deg, radius=0.1)  # degrees

#fig.recenter( ra_dec.ra.deg, ra_dec.dec.deg , width=0.3, height=0.2)  # degrees

fig.show_contour(colors='white', linewidths=1, levels = [0.01, 0.02, 0.05, 0.1, 0.2], smooth=1)

plt.show()