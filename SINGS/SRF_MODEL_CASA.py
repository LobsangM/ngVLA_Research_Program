# name: SFR_MODEL_CASA.py
# date: 13-05-2026

import os
import sys
sys.path.insert(0, '.')

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from casatasks import importfits, imsmooth, exportfits
from utils import get_config, galaxias, REDSHIFT_SFR_TABLE


# -----------------

class SFRRadioModel:

    def __init__(self, config):

        self.config = config

        self._setup_cosmology()
        self._load_image()
        self._scale_sfr()
        self._compute_luminosity()

    # -------------------------------------------------
    def _setup_cosmology(self):

        cosmo_cfg = self.config["cosmology"]

        self.cosmo = FlatLambdaCDM(
            H0=cosmo_cfg["H0"] * u.km / u.s / u.Mpc,
            Om0=cosmo_cfg["Om0"],
            Tcmb0=cosmo_cfg["Tcmb0"] * u.K
        )

    # -------------------------------------------------
    def _load_image(self):

        file = self.config["image"]["fits_file"]

        hdul = fits.open(file)
        data = hdul[0].data
        self.original_header = hdul[0].header

        self.image = self._crop_to_galaxy(data, self.original_header)

        print("Imagen cargada con shape:", self.image.shape)

    # -------------------------------------------------
    def _crop_to_galaxy(self, data, header):
        # Recorta al núcleo real (OBJCTX/OBJCTY del header SINGS) en vez de
        # usar el frame nativo completo. Sin esto, artefactos de borde
        # (rayos cósmicos, estrellas mal sustraídas, gradientes de cielo)
        # se suman al total y dominan el SFR normalizado. Ver README.md
        # sección "Recorte centrado en el núcleo real".
        ny, nx = data.shape
        frac = self.config["image"].get("crop_fraction", 0.7)

        obj_x, obj_y = header.get("OBJCTX"), header.get("OBJCTY")
        try:
            obj_x, obj_y = float(obj_x), float(obj_y)
            valid = 0 <= obj_x <= nx and 0 <= obj_y <= ny
        except (TypeError, ValueError):
            valid = False

        cx, cy = (obj_x, obj_y) if valid else (nx / 2.0, ny / 2.0)
        origen = "OBJCTX/OBJCTY" if valid else "centro geométrico (sin OBJCTX/OBJCTY válido)"

        half = int(min(nx, ny) * frac / 2.0)
        x0, x1 = int(np.clip(cx - half, 0, nx)), int(np.clip(cx + half, 0, nx))
        y0, y1 = int(np.clip(cy - half, 0, ny)), int(np.clip(cy + half, 0, ny))

        print(f"Recorte centrado en ({cx:.1f}, {cy:.1f}) [{origen}]: filas {y0}:{y1}, columnas {x0}:{x1}")

        return data[y0:y1, x0:x1]

    # -------------------------------------------------
    def _scale_sfr(self):

        total_sfr_target = self.config["sfr"]["total_sfr"]

        # Solo el flujo positivo es Hα real; los píxeles negativos vienen de
        # una resta de continuo imperfecta y no representan "SFR negativo".
        # Si se normaliza con la suma cruda (incluye negativos), en imágenes
        # donde el fondo negativo domina (ver README.md "Normalización con
        # flujo positivo") el factor de escala se invierte y todo el mapa
        # queda con el signo al revés.
        positive_flux = np.clip(self.image, 0, None)
        scaling_factor = total_sfr_target / np.sum(positive_flux)

        self.SFR_map = scaling_factor * positive_flux
        self.SFR_map_log = np.log10(1 + self.SFR_map)

        print("SFR total final: %.2f Msun/yr" % np.sum(self.SFR_map))

    # -------------------------------------------------
    def _compute_luminosity(self):

        self.lum_map = self.SFR_map / (4.87e-29)
        self.lum_total = np.sum(self.lum_map)

        print("Luminosidad total: %.2e erg/s/Hz" % self.lum_total)

    # -------------------------------------------------
    def luminosity_to_flux(self):

        z = self.config["radio"]["redshift"]
        alpha = self.config["radio"]["alpha"]
        obs_freq = self.config["radio"]["obs_frequency"]
        rest_freq = self.config["radio"]["rest_frequency"]

        D_L = self.cosmo.luminosity_distance(z).to('cm').value

        factor = (4 * np.pi * D_L**2) / ((1+z)**(1-alpha)) \
                 * (rest_freq/obs_freq)**(-alpha)

        self.flux_map = self.lum_map / factor
        self.flux_total = np.sum(self.flux_map)

        print("Flujo total: %.2e Jy" % (self.flux_total*1e23))

    # -------------------------------------------------
    def plot_sfr(self):

        plt.imshow(self.SFR_map_log, origin="lower",
                   cmap="inferno", norm=LogNorm())
        plt.colorbar(label="log10(1 + SFR)")
        plt.title("Mapa SFR")
        plt.xlim(350, 2450)
        plt.ylim(500, 2650)
        plt.show()

    # -------------------------------------------------
    def plot_flux(self):

        plt.imshow(self.flux_map, origin="lower",
                   cmap="inferno", norm=LogNorm(vmin=1e-35, vmax=1e-32))
        #plt.colorbar(label="Flux Density [erg/s/cm^2/Hz]")
        #plt.title("Flux Map z=%s" % self.config["radio"]["redshift"])
        #plt.xlim(350, 2450)
        #plt.ylim(500, 2650)
        #plt.show()

    # -------------------------------------------------
    def print_statistics(self):

        print("----- Estadísticas -----")
        print("Max flux: %.2e" % np.max(self.flux_map))
        print("Min flux: %.2e" % np.min(self.flux_map))

    # -------------------------------------------------
    def save_flux_to_fits(self):

        z = self.config["radio"]["redshift"]
        total_sfr = self.config["sfr"]["total_sfr"]
        galaxia = self.config["galaxia"]

        folder = "mapas_SFR100"
        os.makedirs(folder, exist_ok=True)

        output_name = os.path.join(folder, "mapa_SFR%s_%s.fits" % (total_sfr, galaxia))

        # conversión de erg/s/cm^2/Hz a Jy
        flux_jy = self.flux_map * 1e23

        hdu = fits.PrimaryHDU(flux_jy)
        hdu.header["REDSHIFT"] = z
        hdu.header["SFR"] = total_sfr
        hdu.header["BUNIT"] = "Jy"

        hdu.writeto(output_name, overwrite=True)

        print("Mapa guardado como:", output_name)



# -------------------------------------------------
    def add_wcs_and_save(self):


        z = self.config["radio"]["redshift"]
        total_sfr = self.config["sfr"]["total_sfr"]

        tamano_fisico_kpc = self.config["image"]["crop_size"]   # 10 kpc
        npix = self.flux_map.shape[0]

      
        # calculadora
        kpc_per_arcsec = self.cosmo.kpc_proper_per_arcmin(z).value / 60.0

        kpc_per_pix = tamano_fisico_kpc / npix
        arcsec_per_pix = kpc_per_pix / kpc_per_arcsec
        deg_per_pix = arcsec_per_pix / 3600.0

        #conversión a jy antes de guardar:
        flux_jy = self.flux_map * 1e23

        hdu = fits.PrimaryHDU(flux_jy)
        header = hdu.header

        print("kpc/arcsec:", kpc_per_arcsec)
        print("kpc/pix:", kpc_per_pix)
        print("arcsec/pix:", arcsec_per_pix)
        print("deg/pix:", deg_per_pix)

       
        # Crear header

        header['CTYPE1'] = 'RA---SIN'
        header['CTYPE2'] = 'DEC--SIN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'

        header['CRPIX1'] = npix / 2.0
        header['CRPIX2'] = npix / 2.0

        header['CRVAL1'] = 150.0
        header['CRVAL2'] = 2.0

        header['CDELT1'] = -deg_per_pix
        header['CDELT2'] = deg_per_pix

        header['BUNIT'] = "Jy"
        header['REDSHIFT'] = z
        header['SFR'] = total_sfr

        header['COMMENT'] = "Simulacion z=%s, %s kpc x %s kpc" % (z, tamano_fisico_kpc, tamano_fisico_kpc)

        folder = self.config.get("output_folder", ".")
        os.makedirs(folder, exist_ok=True)

        output_name = os.path.join(
            folder,
            "mapa_z%s_%skpc_SFR%s.fits" % (z, tamano_fisico_kpc, total_sfr)
        )

        hdu.writeto(output_name, overwrite=True)

        print("Imagen con header guardado como:", output_name)



    # -------------------------------------------------
    def casa_convolution(self):

        z = self.config["radio"]["redshift"]
        size = self.config["image"]["crop_size"]
        total_sfr = self.config["sfr"]["total_sfr"]

        major = self.config["radio"]["beam_major"]
        minor = self.config["radio"]["beam_minor"]
        pa = self.config["radio"]["beam_pa"]

        output_folder = self.config.get("output_folder", ".")
        input_fits = os.path.join(output_folder, "mapa_z%s_%skpc_SFR%s.fits" % (z, size, total_sfr))

        casa_image = "temp_%s.image" % z
        smooth_image = "temp_%s_smooth.image" % z

        # quitar el nombre para el label
        beam_label = major.replace(" ", "")
        galaxia = self.config["galaxia"]

        folder = "Convoluciones"
        os.makedirs(folder, exist_ok=True)

        output_fits = os.path.join(
            folder,
            "convolucion_%s_%s.fits" % (beam_label, galaxia)
        )

        print("Convolución:")
        print("Beam mayor:", major)
        print("Beam menor:", minor)

        # FITS → CASA
        importfits(fitsimage=input_fits,
                imagename=casa_image,
                overwrite=True)

        # Convolución
        imsmooth(imagename=casa_image,
                kernel="gauss",
                major=major,
                minor=minor,
                pa=pa,
                outfile=smooth_image,
                overwrite=True)

        # CASA → FITS
        exportfits(imagename=smooth_image,
                fitsimage=output_fits,
                overwrite=True)

        print("Archivo final:", output_fits)

# MAIN

if __name__ == "__main__":
    for entry in REDSHIFT_SFR_TABLE:
        z     = entry["redshift"]
        sfr   = entry["sfr"]
        label = entry["label"]

        print(f"\n{'='*55}")
        print(f"  z={z}  SFR={sfr} M☉/yr  [{label}]  ({len(galaxias)} galaxias)")
        print(f"{'='*55}")

        for galaxia in galaxias:
            print(f"\n  → {galaxia}")

            config = get_config(galaxia, redshift=z, sfr=sfr, label=label)
            model  = SFRRadioModel(config)
            model.luminosity_to_flux()
            model.print_statistics()
            model.add_wcs_and_save()

    print("\nTodos los mapas base generados.")