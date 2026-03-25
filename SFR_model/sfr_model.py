import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import yaml


class SFRRadioModel:

    def __init__(self, config_file):

        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

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
        self.image = hdul[0].data
        self.original_header = hdul[0].header
        print(f"Imagen cargada con shape: {self.image.shape}")


    # -------------------------------------------------
    def _scale_sfr(self):

        total_sfr_target = self.config["sfr"]["total_sfr"]

        scaling_factor = total_sfr_target / np.sum(self.image)

        self.SFR_map = scaling_factor * self.image
        self.SFR_map_log = np.log10(1 + self.SFR_map)

        print(f"SFR total final: {np.sum(self.SFR_map):.2f} Msun/yr")

    # -------------------------------------------------
    def _compute_luminosity(self):

        self.lum_map = self.SFR_map / (4.87e-29)
        self.lum_total = np.sum(self.lum_map)

        print(f"Luminosidad total: {self.lum_total:.2e} erg/s/Hz")

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

        print(f"Flujo total: {self.flux_total*1e23:.2e} Jy")

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
        plt.colorbar(label="Flux Density [erg/s/cm^2/Hz]")
        plt.title(f"Flux Map z={self.config['radio']['redshift']}")
        plt.xlim(350, 2450)
        plt.ylim(500, 2650)
        plt.show()

    # -------------------------------------------------
    def print_statistics(self):

        print("----- Estadísticas -----")
        print(f"Max flux: {np.max(self.flux_map):.2e}")
        print(f"Min flux: {np.min(self.flux_map):.2e}")

    # -------------------------------------------------
    def save_flux_to_fits(self):

        z = self.config["radio"]["redshift"]
        total_sfr = self.config["sfr"]["total_sfr"]

        output_name = f"flux_map_z{z}_SFR{total_sfr}.fits"

        hdu = fits.PrimaryHDU(self.flux_map)
        hdu.header["REDSHIFT"] = z
        hdu.header["SFR"] = total_sfr
        hdu.header["BUNIT"] = "erg/s/cm^2/Hz"

        hdu.writeto(output_name, overwrite=True)

        print(f"Mapa guardado como: {output_name}")