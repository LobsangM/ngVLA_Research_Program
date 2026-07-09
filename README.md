# ngVLA_Research_Program

Pipeline de simulación en radioastronomía que predice cómo se verían galaxias del catálogo **SINGS** (*Spitzer Infrared Nearby Galaxies Survey*) si fueran observadas con los radiotelescopios **ngVLA** y **VLA**. El flujo de trabajo transforma mapas de flujo observados (Hα / MIPS) en imágenes de radio sintéticas a distintos corrimientos al rojo (z = 0.4, 1.0, 2.0, 3.0, 5.0), aplicando convolución de haz (*beam*) realista y ruido correlacionado.

## Contenido

- [Requisitos y entorno](#requisitos-y-entorno)
- [Estructura del repositorio](#estructura-del-repositorio)
- [Cómo ejecutar el pipeline](#cómo-ejecutar-el-pipeline)
- [Arquitectura del procesamiento](#arquitectura-del-procesamiento)
- [Física del modelo](#física-del-modelo)
- [Configuraciones de telescopio](#configuraciones-de-telescopio)
- [Tabla z × SFR](#tabla-z--sfr)
- [Notas y problemas conocidos](#notas-y-problemas-conocidos)

## Requisitos y entorno

- Python 3.10 (entorno virtual en `venv310/`)
- **CASA** (`casatasks`/`casatools` 6+): los scripts que convolucionan imágenes llaman a `importfits`, `imsmooth`, `exportfits` y fallarán fuera de un entorno Python de CASA.
- Paquetes clave: `numpy`, `astropy`, `matplotlib`, `scipy`, `aplpy`, `tomli` (ver `requirements.txt` para versiones exactas).

```bash
source venv310/bin/activate
pip install -r requirements.txt
```

## Estructura del repositorio

```
ResearchP_ngVLA/
├── Modelo/                         # Prototipo de una sola galaxia (NGC 6946)
│   ├── sfr_model_CASA.py           # Clase SFRRadioModel + ejecución para NGC 6946
│   ├── ngVLA_config_{A,B}.py       # Simulación ngVLA (A: 0.293", B: 0.961")
│   ├── VLA_config_{A,B}.py         # Simulación VLA (A: 0.293", B: 0.961")
│   ├── config_main.py              # Orquestador: corre todo el pipeline en secuencia
│   ├── rebin_beam_sampling.py      # Rebinning de píxeles vía CASA (importfits/imrebin/exportfits)
│   └── ngc6946Ha_I_Ha_ksb2004.fits # Imagen Hα de entrada
└── SINGS/                           # Pipeline en lote (14 galaxias con FITS disponible)
    ├── SRF_MODEL_CASA.py           # Genera mapas de flujo: todas las galaxias × todos los z
    ├── utils.py                    # Carga TOML → lista `galaxias`, `get_config()`, `REDSHIFT_SFR_TABLE`
    ├── main.py                     # Orquestador: corre las 5 etapas en secuencia
    ├── collage.py                  # Collage PNG por galaxia (RMS + contornos 3σ)
    ├── preview_Ha.py               # Previews PNG de las imágenes Hα crudas
    ├── tomls/                      # Configuración por galaxia (Ha_SUB.toml, MIPS70.toml, MIPS160.toml)
    ├── imagenes_Ha/                # Entrada: {galaxia}_HA_SUB_dr4.fits (14 archivos)
    ├── resultados/{z}/{galaxia}/   # Mapas de flujo con WCS: mapa_{z}_30kpc_SFR{sfr}.fits
    ├── ngVLA_config_A/{z}/         # Salida: ngVLA_A_{galaxia}_noise.fits + .png
    ├── ngVLA_config_B/{z}/         # Salida: ngVLA_B_{galaxia}_noise.fits + .png
    ├── VLA_config_A/{z}/           # Salida: VLA_A_{galaxia}_noise.fits + .png
    ├── VLA_config_B/{z}/           # Salida: VLA_B_{galaxia}_noise.fits + .png
    ├── collage/                    # Salida de collage.py: {galaxia}.png
    └── Convoluciones/              # Imágenes convolucionadas: convolucion_{beam}_{galaxia}.fits
```

## Cómo ejecutar el pipeline

Los scripts deben ejecutarse **desde su propio directorio**, ya que usan rutas relativas hacia los archivos FITS.

### Prototipo de una sola galaxia (`Modelo/`)

```bash
cd Modelo
python sfr_model_CASA.py        # genera la convolución para NGC 6946 (mapa_base)
python ngVLA_config_A.py        # simula observación ngVLA-A + ruido
python ngVLA_config_B.py        # simula observación ngVLA-B + ruido
python VLA_config_A.py          # simula observación VLA-A + ruido
python VLA_config_B.py          # simula observación VLA-B + ruido
python config_main.py           # orquestador: corre sfr_model_CASA.py + los cuatro scripts de telescopio
```

### Pipeline en lote SINGS (14 galaxias)

```bash
cd SINGS
python main.py                  # orquestador: corre los cinco scripts siguientes en secuencia
# o paso a paso:
python SRF_MODEL_CASA.py        # paso 1: mapas de flujo para todas las galaxias × todos los z
python ngVLA_config_A.py        # paso 2a: convolución + ruido, config ngVLA-A
python ngVLA_config_B.py        # paso 2b: convolución + ruido, config ngVLA-B
python VLA_config_A.py          # paso 2c: convolución + ruido, config VLA-A
python VLA_config_B.py          # paso 2d: convolución + ruido, config VLA-B
python collage.py               # opcional: collage PNG por galaxia (con RMS + contornos 3σ)
python preview_Ha.py            # opcional: previews PNG de las imágenes Hα crudas
```

### Rebinning de imágenes CASA (dentro del entorno CASA)

```bash
# Modelo/rebin_beam_sampling.py contiene tareas de CASA (importfits, imrebin, exportfits)
# Debe ejecutarse con: casa --nogui -c rebin_beam_sampling.py
```

## Arquitectura del procesamiento

Pipeline por galaxia × por redshift:

```
FITS de entrada (imagen Hα/MIPS)
    ↓  escalar valores de píxel → mapa de SFR [M☉/año]        (_scale_sfr)
    ↓  SFR → luminosidad de radio [erg/s/Hz]                   (_compute_luminosity)
    ↓  luminosidad → densidad de flujo [Jy] a z                 (luminosity_to_flux)
    ↓  agregar cabecera WCS (escala de píxel basada en kpc)      (add_wcs_and_save)
         → resultados/{z_label}/{galaxia}/mapa_{z_label}_30kpc_SFR{sfr}.fits
    ↓  convolución de haz gaussiano vía CASA                     (casa_convolution / imsmooth)
         → Convoluciones/convolucion_{beam}_{galaxia}.fits
    ↓  reescalar Jy/píxel → Jy/beam, luego agregar ruido           (add_correlated_noise)
Salida: {config}/{z_label}/{prefijo}_{galaxia}_noise.fits  +  *_noise.png
```

El bucle externo de `SRF_MODEL_CASA.py` (versión SINGS) itera: `for entry in REDSHIFT_SFR_TABLE: for galaxy in galaxias`.

La clase `SFRRadioModel` existe en dos versiones: `Modelo/sfr_model_CASA.py` (config fija, una sola galaxia: NGC 6946) y `SINGS/SRF_MODEL_CASA.py` (config por galaxia vía `get_config()`, escribe a `resultados/{z_label}/{galaxia}/`). Los scripts de telescopio en `Modelo/` importan `SFRRadioModel` y `config` directamente desde `sfr_model_CASA` en el mismo directorio.

## Física del modelo

- Relación radio-SFR: `L [erg/s/Hz] = SFR / 4.87e-29`
- Índice espectral: α = 0.7
- Frecuencia de observación: 10 GHz (banda ngVLA); frecuencia en reposo: 1.4 GHz
- Flujo: `F = L / (4π D_L² · (1+z)^(1-α) · (ν_rest/ν_obs)^(-α))`
- Ruido correlacionado: ruido blanco con `std = sigma_beam`, convolucionado con un kernel gaussiano normalizado en L2 (`kernel /= sqrt(sum(kernel**2))`), de forma que el ruido correlacionado resultante conserva `std ≈ sigma_beam` por píxel
- `crop_size` es 30 kpc tanto en `Modelo/` como en `SINGS/`

### Conversión de unidades crítica en `add_correlated_noise`: Jy/píxel → Jy/beam

`imsmooth` de CASA conserva la densidad de flujo total (no sabe que la entrada es un mapa por píxel sin haz restaurador, así que trata la convolución como conservadora de flujo). Esto diluye el brillo superficial de pico por un factor ~(área del beam / área del píxel) — a menudo 10⁴–10⁵ para estas imágenes, ya que los mapas están fuertemente sobremuestreados respecto al FWHM del beam. Los valores de ruido del telescopio (`noise_nJy`/`noise_uJy`), en cambio, están definidos **por beam** (convención estándar en interferometría radio: flujo integrado sobre un beam sintetizado).

Antes de agregar ruido, `add_correlated_noise` debe reescalar la señal convolucionada desde esta cantidad diluida tipo Jy/píxel a Jy/beam:

```python
beam_area_pix = (np.pi / (4.0 * np.log(2.0))) * (beam_arcsec / pixel_arcsec) ** 2
data = data * beam_area_pix
```

Sin este reescalado, la señal termina 3-5 órdenes de magnitud por debajo del piso de ruido y la fuente es completamente invisible, aunque el ruido en sí se genere correctamente. Esta lógica está duplicada en los ocho scripts de configuración de telescopio (cuatro en `Modelo/`, cuatro en `SINGS/`) — si se modifica la función de ruido, hay que aplicar el mismo fix en todos.

Incluso con el fix aplicado, las detecciones pueden ser genuinamente marginales: la config A (beam angosto, 0.293″) sacrifica sensibilidad por resolución, por lo que puede detectar marginalmente solo el nudo de formación estelar más brillante de una galaxia (SNR ~3-4σ) en vez del disco completo, mientras que la config B (beam ancho, 0.961″) integra más flujo por beam y muestra SNR más saludable (~7-8σ) para la misma fuente. Este es el comportamiento interferométrico esperado (un beam fino tiene poca sensibilidad de brillo superficial a emisión extendida), no un bug — ver `reportes/dudas.tex` para el razonamiento estadístico completo (RMS vía sigma-clipping, por qué los contornos 3σ son la forma correcta de separar emisión real de ruido, y por qué un estiramiento de color lineal min/max ingenuo en `plot_with_beam` es engañoso para detecciones marginales).

## Configuraciones de telescopio

| Config    | Beam (arcsec) | Ruido          |
|-----------|---------------|----------------|
| ngVLA-A   | 0.293         | 28.44 nJy/beam |
| ngVLA-B   | 0.961         | 32.68 nJy/beam |
| VLA-A     | 0.293         | 0.56 µJy/beam  |
| VLA-B     | 0.961         | 0.56 µJy/beam  |

## Tabla z × SFR

Definida en `SINGS/utils.py` como `REDSHIFT_SFR_TABLE`, basada en Leslie et al. 2020 (M★ = 10¹⁰ M☉):

| Label | z   | SFR (M☉/año) |
|-------|-----|--------------|
| z0.4  | 0.4 | 1.5          |
| z1.0  | 1.0 | 7            |
| z2.0  | 2.0 | 40           |
| z3.0  | 3.0 | 100          |
| z5.0  | 5.0 | 300          |

Esta tabla no es arbitraria: registra la SFR típica de secuencia principal a masa estelar fija en cada época, por lo que los bins a mayor z son intrínsecamente mucho más luminosos. Esa ganancia de luminosidad compensa aproximadamente (pero no exactamente) el atenuamiento cosmológico, razón por la cual el flujo integrado total se mantiene en el mismo orden de magnitud a través de los cinco redshifts para una galaxia dada.

## Notas y problemas conocidos

- **Directorios temporales de CASA**: `SRF_MODEL_CASA.py` y los scripts de telescopio dejan directorios `temp_{z}.image` y `temp_{z}_smooth.image`. Son seguros de borrar después de una corrida.
- **Unidades de ruido distintas**: las configs ngVLA pasan el ruido en **nJy** (`noise_nJy * 1e-9`) mientras que las configs VLA lo pasan en **µJy** (`noise_uJy * 1e-6`). Las firmas de `add_correlated_noise` difieren en consecuencia — no son intercambiables aunque se vean idénticas.
- **14 vs 52 galaxias**: `Ha_SUB.toml` lista 52 galaxias, pero `galaxias` en `utils.py` se filtra al cargar a solo las que tienen un FITS coincidente en `imagenes_Ha/` — actualmente 14.
- **`mapas_SFR100/` y `Convoluciones/`**: escritos por `save_flux_to_fits()` y `casa_convolution()` respectivamente, pero estos métodos no se llaman en el bucle principal actual (`__main__` solo llama a `add_wcs_and_save()`). Esos directorios son intermedios de versiones anteriores del pipeline.
- **Duplicación de código**: `add_correlated_noise` y `plot_with_beam` están copiados literalmente en los cuatro scripts de telescopio de `Modelo/` y los cuatro de `SINGS/`. Cambios al modelo de ruido o al estilo de gráfico deben aplicarse en los ocho archivos.
