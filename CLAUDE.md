# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Radio astronomy simulation pipeline that predicts how SINGS (Spitzer Infrared Nearby Galaxies Survey) galaxies would appear when observed with the ngVLA and VLA radio telescopes. The core workflow transforms observed Hα/MIPS flux maps into synthetic radio images at z=1 with realistic beam convolution and noise.

## Running the Pipeline

Scripts must be run **from within their own directory** because they use relative paths for FITS files.

**Single-galaxy prototype (`Modelo/`):**
```bash
cd Modelo
python sfr_model_CASA.py        # generates convolution for NGC 6946 (mapa_base)
python ngVLA_config_A.py        # simulate ngVLA-A observation + noise
python ngVLA_config_B.py        # simulate ngVLA-B observation + noise
python VLA_config_A.py          # simulate VLA-A observation + noise
python VLA_config_B.py          # simulate VLA-B observation + noise
python config_main.py           # orchestrator: runs all four telescope scripts sequentially
```

**SINGS batch pipeline (52 galaxies):**
```bash
cd SINGS
python SRF_MODEL_CASA.py        # step 1: generate flux maps for all galaxies
python ngVLA_config_A.py        # step 2a: convolve + noise, ngVLA-A config
python ngVLA_config_B.py        # step 2b: convolve + noise, ngVLA-B config
python VLA_config_A.py          # step 2c: convolve + noise, VLA-A config
python VLA_config_B.py          # step 2d: convolve + noise, VLA-B config
python preview_Ha.py            # optional: generate PNG previews + collage of raw Hα inputs
```

**CASA image rebinning (run inside CASA environment):**
```bash
# Modelo/rebin_beam_sampling.py contains CASA task calls (importfits, imrebin, exportfits)
# Must be executed via: casa --nogui -c rebin_beam_sampling.py
```

## Environment and Dependencies

- Python 3.10 virtual environment: `venv310/`
- Activate: `source venv310/bin/activate`
- Key packages: `numpy`, `astropy`, `matplotlib`, `scipy`, `aplpy`, `tomli`
- **CASA dependency**: Scripts that convolve images require `casatasks` (from CASA 6+). These scripts call `importfits`, `imsmooth`, `exportfits` — they will fail outside a CASA Python environment.

## Architecture

### Processing pipeline (per galaxy)

```
Input FITS (Hα/MIPS image)
    ↓  scale pixel values → SFR map [M☉/yr]       (_scale_sfr)
    ↓  SFR → radio luminosity [erg/s/Hz]           (_compute_luminosity)
    ↓  luminosity → flux density [Jy] at z=1       (luminosity_to_flux)
    ↓  add WCS header (kpc-based pixel scale)       (add_wcs_and_save)
    ↓  Gaussian beam convolution via CASA           (casa_convolution / imsmooth)
    ↓  add beam-correlated noise                    (add_correlated_noise)
Output: *_noise.fits  +  *_noise.png
```

### Key physics constants

- Radio-SFR relation: `L [erg/s/Hz] = SFR / 4.87e-29`
- Spectral index: α = 0.7
- Observing frequency: 10 GHz (ngVLA band); rest frequency: 1.4 GHz
- Flux: `F = L / (4π D_L² · (1+z)^(1-α) · (ν_rest/ν_obs)^(-α))`
- Correlated noise: white noise convolved with Gaussian kernel sized to match beam FWHM
- `crop_size` is 30 kpc in both `Modelo/` and `SINGS/`

### Telescope configurations

| Config    | Beam (arcsec) | Noise          |
|-----------|--------------|----------------|
| ngVLA-A   | 0.293        | 28.44 nJy/beam |
| ngVLA-B   | 0.961        | 32.68 nJy/beam |
| VLA-A     | 0.293        | 0.56 µJy/beam  |
| VLA-B     | 0.961        | 0.56 µJy/beam  |

### Directory layout

```
ResearchP_ngVLA/
├── Modelo/
│   ├── sfr_model_CASA.py          # SFRRadioModel class + single-galaxy run (NGC 6946)
│   ├── ngVLA_config_{A,B}.py      # telescope simulations (import sfr_model_CASA)
│   ├── VLA_config_{A,B}.py
│   ├── config_main.py             # orchestrator: runs all four telescope scripts
│   ├── rebin_beam_sampling.py     # CASA pixel rebinning (run via casa --nogui)
│   ├── ngc6946Ha_I_Ha_ksb2004.fits  # input FITS for single-galaxy prototype
│   └── resultados/                # output per telescope config
└── SINGS/
    ├── SRF_MODEL_CASA.py      # batch pipeline, iterates all galaxias
    ├── utils.py               # loads TOML → galaxias list + get_config()
    ├── tomls/
    │   ├── Ha_SUB.toml        # 52 galaxies, Hα band
    │   ├── MIPS70.toml        # same galaxies, 70µm band
    │   └── MIPS160.toml       # same galaxies, 160µm band
    ├── ngVLA_config_{A,B}.py  # batch telescope simulations
    ├── VLA_config_{A,B}.py
    ├── preview_Ha.py          # generates individual PNGs + collage of raw Hα inputs → imagenes_png/
    ├── imagenes_Ha/           # input: *_HA_SUB_dr4.fits per galaxy
    ├── imagenes_png/          # output of preview_Ha.py: {galaxy}.png + collage_Ha_{cmap}.png
    ├── mapas_SFR100/          # intermediate: mapa_SFR100_{galaxy}.fits (save_flux_to_fits)
    ├── resultados/{galaxy}/   # intermediate: mapa_z1_30kpc_SFR100.fits (add_wcs_and_save)
    ├── ngVLA_config_A/        # output: convolucion_{beam}arcsec_{galaxy}_noise.fits + .png
    ├── ngVLA_config_B/
    ├── VLA_config_A/
    ├── VLA_config_B/
    └── Convoluciones/         # convolved images: convolucion_0.05arcsec_{galaxy}.fits
```

### SFRRadioModel class (two versions)

Both `Modelo/sfr_model_CASA.py` and `SINGS/SRF_MODEL_CASA.py` define `SFRRadioModel`. The `Modelo/` version uses a single hardcoded config dict and NGC 6946 as input. The SINGS version adds per-galaxy config via `get_config()` (from `utils.py`) and writes outputs to `resultados/{galaxy}/`. The `Modelo/` telescope scripts import `SFRRadioModel` and `config` directly from `sfr_model_CASA` in the same directory.

### Switching photometric bands

The active band is controlled by `TOML_PATH` in `SINGS/utils.py`. Change it to `MIPS70.toml` or `MIPS160.toml` instead of `Ha_SUB.toml`, then update the `imagenes_Ha/` path and file suffix in `get_config()`.

### Code duplication in telescope config scripts

`add_correlated_noise` and `plot_with_beam` are copy-pasted verbatim into all four `Modelo/` telescope scripts and all four `SINGS/` telescope scripts. Changes to the noise model or plot style must be applied to all eight files.

### Known issues

- **CASA temp directories**: `SRF_MODEL_CASA.py` and the telescope scripts leave behind `temp_*.image` and `temp_*_smooth.image` CASA directories. These are safe to delete after a run.
- **Noise units differ**: ngVLA configs pass noise in **nJy** (`noise_nJy * 1e-9`) while VLA configs pass noise in **µJy** (`noise_uJy * 1e-6`). The `add_correlated_noise` signatures differ accordingly — they are not interchangeable despite looking identical.
