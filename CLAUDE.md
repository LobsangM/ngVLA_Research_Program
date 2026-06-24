# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Radio astronomy simulation pipeline that predicts how SINGS (Spitzer Infrared Nearby Galaxies Survey) galaxies would appear when observed with the ngVLA and VLA radio telescopes. The core workflow transforms observed Hα/MIPS flux maps into synthetic radio images at z=1 with realistic beam convolution and noise.

## Running the Pipeline

Scripts must be run **from within their own directory** because they use relative paths for FITS files.

**Single-galaxy prototype (root directory):**
```bash
cd /path/to/ResearchP_ngVLA
python sfr_model_CASA.py        # generates mapa_z1_30kpc_SFR100.fits and convolution
python ngVLA_config_A.py        # simulate ngVLA-A observation + noise
python ngVLA_config_B.py        # simulate ngVLA-B observation + noise
python VLA_config_A.py          # simulate VLA-A observation + noise
python VLA_config_B.py          # simulate VLA-B observation + noise
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

**Refactored standalone model (YAML config):**
```bash
cd SFR_model
python main.py                  # reads config.yaml — NOTE: main.py has a broken import (see Known Issues)
```

**CASA image rebinning (run inside CASA environment):**
```bash
# rebin_beam_sampling.py contains CASA task calls (importfits, imrebin, exportfits)
# Must be executed via: casa --nogui -c rebin_beam_sampling.py
```

## Environment and Dependencies

- Python 3.10 virtual environment: `venv310/`
- Activate: `source venv310/bin/activate`
- Key packages: `numpy`, `astropy`, `matplotlib`, `scipy`, `aplpy`, `tomli`, `yaml`
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
- `crop_size` is 30 kpc in the SINGS batch pipeline and 10 kpc in the `SFR_model/` prototype

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
├── sfr_model_CASA.py          # single-galaxy pipeline (root prototype)
├── ngVLA_config_{A,B}.py      # root-level telescope simulations
├── VLA_config_{A,B}.py
├── rebin_beam_sampling.py     # CASA pixel rebinning script
├── SFR_model/
│   ├── main.py                # entry point using config.yaml
│   ├── sfr_model.py           # SFRRadioModel class (YAML-based version)
│   └── config.yaml            # cosmology / image / sfr / radio params
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
    ├── imagenes_MIPS70/       # input: *_mips70_image_v5-0.fits
    ├── imagenes_MIPS160/      # input: *_mips160_image_v5-0.fits
    ├── imagenes_png/          # output of preview_Ha.py: {galaxy}.png + collage_Ha_{cmap}.png
    ├── mapas_SFR100/          # intermediate: mapa_SFR100_{galaxy}.fits (save_flux_to_fits)
    ├── resultados/{galaxy}/   # intermediate: mapa_z1_30kpc_SFR100.fits (add_wcs_and_save)
    ├── ngVLA_config_A/        # output: {TelescopeConfig}_{galaxy}.fits + _noise.fits + .png
    ├── ngVLA_config_B/
    ├── VLA_config_A/
    ├── VLA_config_B/
    └── Convoluciones/         # convolved images: convolucion_0.05arcsec_{galaxy}.fits
```

### Switching photometric bands

The active band is controlled by `TOML_PATH` in `SINGS/utils.py`. Change the import to use `MIPS70.toml` or `MIPS160.toml` instead of `Ha_SUB.toml`, then update `imagenes_Ha/` references in `get_config()` to the correct image folder and file suffix.

### SFRRadioModel class (two versions)

Both `sfr_model_CASA.py` (root) and `SINGS/SRF_MODEL_CASA.py` define a `SFRRadioModel` class. The SINGS version adds per-galaxy config via `get_config()` and writes outputs to `resultados/{galaxy}/`. The `SFR_model/sfr_model.py` is an older YAML-based variant that does not include CASA convolution.

### Code duplication in SINGS telescope config scripts

`add_correlated_noise` and `plot_with_beam` are copy-pasted verbatim into all four SINGS telescope scripts (`ngVLA_config_{A,B}.py`, `VLA_config_{A,B}.py`). They are not shared via a module. Changes to the noise model or plot style must be applied to all four files.

### Root-level experiment scripts

`agregar_ruido.py`, `ruido_correlacionado.py`, and `ruido_por_cuadrante.py` are standalone one-off experiments used during early development (operate on `convolucion_0.05arcsec_rebin.fits`). They are not part of the active pipeline. `SFR_model/suavizado.py` is a one-off WCS header utility that patches an existing FITS file with kpc-based pixel scale.

### Known issues

- **`SFR_model/main.py` has a broken import**: it does `from rebin_beam_sampling import SFRRadioModel` but `rebin_beam_sampling.py` does not define `SFRRadioModel`. The correct import should be `from sfr_model import SFRRadioModel`.
- **CASA temp directories**: `SRF_MODEL_CASA.py` and the telescope scripts leave behind `temp_*.image` and `temp_*_smooth.image` CASA directories in `SINGS/`. These are safe to delete after a run.
- The noise units differ between scripts: ngVLA configs pass noise in **nJy** (`noise_nJy * 1e-9`) while VLA configs pass noise in **µJy** (`noise_uJy * 1e-6`). The function signatures differ accordingly — `add_correlated_noise` is not interchangeable across the four scripts despite looking identical.
