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

### Recorte centrado en el núcleo real (`OBJCTX`/`OBJCTY`)

**Problema detectado**: en varios collages de `SINGS/` (ngc0628, ngc1097, ngc3521, ngc3627, ngc4254, ngc4321, ngc4536, ngc4594, ngc5194, ngc7331) la emisión aparecía desplazada del centro (a veces en una esquina) o la config ngVLA-A se veía dominada por ruido en vez de parecerse a ngVLA-B.

**Causa raíz**: `_scale_sfr()` normaliza la SFR total sumando el frame Hα nativo **completo** (`scaling_factor = total_sfr_target / np.sum(self.image)`), sin recortar nunca la imagen alrededor de la galaxia. `crop_size` (30 kpc) nunca recortaba nada — solo se usaba después, en `add_wcs_and_save()`, para etiquetar la escala de píxel del array completo. Esto significa que cualquier artefacto fuera del disco (rayos cósmicos, estrellas de campo mal sustraídas en la resta de continuo, gradientes residuales de cielo) se suma al mismo total que la emisión real de la galaxia y compite por la misma "porción" de SFR asignada.

Verificación sobre los FITS crudos de `imagenes_Ha/`: al recortar al 70% central del frame (centrado en el núcleo real), las galaxias reportadas como problemáticas retienen mucho menos flujo dentro del recorte que las que ya se veían bien:

| Galaxia   | % de flujo dentro del recorte 70% | Estado reportado |
|-----------|-----------------------------------|-------------------|
| ngc7331   | 25.2%                              | emisión indistinguible |
| ngc5194   | 40.5%                              | mucho ruido |
| ngc1097   | 47.6%                              | emisión en la esquina |
| ngc4536   | 49.8%                              | fuera de centro |
| ngc3627   | 55.6%                              | mucho ruido, fuera de centro |
| ngc0628   | 58.5%                              | mucho ruido en config A |
| ngc4594   | 65.8%                              | fuera de centro (además: `OBJCTX`/`OBJCTY` corruptos en el header, ver abajo) |
| ngc3521   | 63.8%                              | mucho ruido, fuera de centro |
| ngc4254   | 72.7%                              | fuera de centro |
| ngc4321   | 69.0%                              | fuera de centro |
| ngc5055   | 75.3%                              | **bien** |
| ngc3034   | 66.3%                              | **perfecta** (núcleo tan brillante que domina igual) |
| ngc3184   | 92.2%                              | **bien** |
| ngc3938   | 98.4%                              | **bien** |

Es decir: el 25-55% del "presupuesto de SFR" de las galaxias con problemas se estaba gastando fuera del disco real, diluyendo la señal verdadera (especialmente notorio en la config ngVLA-A, cuyo beam angosto ya de por sí tiene poca sensibilidad de brillo superficial — ver sección anterior).

**Fix aplicado**: `_load_image()` (en `Modelo/sfr_model_CASA.py` y `SINGS/SRF_MODEL_CASA.py`) ahora llama a un nuevo método `_crop_to_galaxy()` que recorta el array de entrada a una región cuadrada centrada en el núcleo real de la galaxia, de tamaño `crop_fraction * min(nx, ny)` (nuevo parámetro en config, default `0.7`), **antes** de que `_scale_sfr()` normalice el flujo. El centro se toma de las cabeceras `OBJCTX`/`OBJCTY` del FITS de entrada (posición real del núcleo en el header original SINGS/CTIO); si faltan o caen fuera del frame (como en `ngc4594`, cuyo header trae valores corruptos: `OBJCTX=-32698.4`, `OBJCTY=292992.4`, muy fuera de los límites de la imagen 2060×2062), se usa el centro geométrico del array como fallback. El resto del pipeline (`_scale_sfr`, `luminosity_to_flux`, `add_wcs_and_save`, convolución, ruido) no necesita cambios: opera sobre `self.image`/`self.flux_map` ya recortado, y el WCS sigue etiquetando el array resultante como `crop_size` kpc, ahora correctamente centrado en la galaxia.

`crop_fraction=0.7` es un valor conservador elegido para no recortar de más galaxias con disco genuinamente extendido (ej. ngc5055 retiene 75.3% de su flujo real con este valor; con `crop_fraction=0.5` cae a 29.4%, señal de que empieza a cortar disco real). Es un parámetro ajustable por galaxia si algún caso lo requiere (ej. ngc7331, muy elongada, o ngc5194/M51 con su compañera NGC 5195 dentro del mismo frame) — no hay una distancia real por galaxia en el pipeline para convertir `crop_size` en kpc físicos verdaderos, así que el recorte se define en fracción del frame nativo, no en kpc reales de la galaxia local.

### Normalización con flujo positivo (bug de inversión de signo)

Tras aplicar el recorte de arriba, ngc4254 y ngc4321 seguían mal: el "Mapa radio" se veía invertido — el fondo aparecía brillante y el disco de la galaxia (visible y bien centrado en Hα original) aparecía como un hueco oscuro.

**Causa raíz**: `_scale_sfr()` normalizaba con `scaling_factor = total_sfr_target / np.sum(self.image)`, usando la suma cruda del array (positivos **y** negativos). Los FITS Hα-restados de continuo tienen fondo de cielo ligeramente negativo por una resta de continuo imperfecta; en la mayoría de las galaxias esto es una fracción pequeña que no cambia el signo del total, pero en ngc4254 y ngc4321 el fondo negativo domina lo suficiente como para que la suma completa del frame sea **negativa**, incluso sin recortar:

```
ngc4254: sum = -1.324e+05   (97.4% de los píxeles negativos)
ngc4321: sum = -7.746e+03   (68.9% de los píxeles negativos)
```

Con `scaling_factor` negativo, todo `SFR_map` queda con el signo invertido: los nudos de formación estelar reales (positivos en el crudo) terminan con SFR negativo (se muestran como el "hueco"), mientras el ruido de fondo (negativo en el crudo) se vuelve positivo y domina la imagen. Este bug ya existía antes del fix del recorte — simplemente quedaba oculto porque, sin recortar, el resto del frame (fuera de la galaxia) a veces compensaba la suma lo suficiente para que no cambiara de signo.

**Fix aplicado**: en `_scale_sfr()` (ambas copias), se recorta el array a flujo no-negativo antes de normalizar y de construir `SFR_map`:

```python
positive_flux = np.clip(self.image, 0, None)
scaling_factor = total_sfr_target / np.sum(positive_flux)
self.SFR_map = scaling_factor * positive_flux
```

Esto es correcto físicamente (el Hα real no puede ser negativo; los píxeles negativos son ruido de resta de continuo, no "SFR negativo") y no afecta a las demás galaxias, cuyo total ya era positivo y seguirá siéndolo.

### Máscaras manuales de estrellas mal sustraídas (ngc3627, ngc4594, ngc5055)

Tras aplicar el recorte y la normalización con flujo positivo, tres galaxias seguían mostrando una fuente compacta espuria separada del núcleo real: **ngc5055**, y — al revisar el resto del lote con el mismo criterio — también **ngc3627** y **ngc4594**. En los tres casos el pico de flujo dentro del recorte es un múltiplo extremo (196x a 1431x) del percentil 99.9 del resto de la imagen, y cae exactamente donde el propio equipo de SINGS marcó un círculo sobre la imagen Hα original (una estrella de campo mal sustraída, ya señalada como mala en la reducción). Los valores crudos alrededor muestran el patrón dipolar típico de una resta de continuo mal calibrada para esa estrella (ejemplo, ngc5055):

```
  9960   59272  55340  -3800
  -933   27861   8150  -4314
```

**Por qué no se automatizó con una regla estadística**: se probó sigma-clipping genérico (rechazar píxeles > N-sigma sobre la distribución positiva) pero resultó demasiado agresivo — en imágenes con núcleos reales muy compactos y brillantes (ej. ngc3034/M82, que funciona perfecto) el mismo filtro estadístico recorta la fuente real junto con el artefacto, porque ambos son "outliers" de la misma manera. No hay forma de distinguir "estrella mal sustraída" de "núcleo real muy brillante" con una regla global sin arriesgar romper las galaxias que ya están bien.

**Fix aplicado**: máscara manual por galaxia, con coordenadas fijas en el sistema de la imagen **original** (sin recortar), aplicada en `_load_image()` **antes** del recorte (`_mask_known_artifacts()`, `SINGS/SRF_MODEL_CASA.py`) — así la máscara no depende de `crop_fraction`. Las coordenadas y radios están en `SINGS/utils.py` (`STAR_MASKS`), encontrados localizando el pico de flujo dentro del recorte y verificando con un perfil radial dónde el valor vuelve al nivel de fondo (para no recortar de más ni de menos):

| Galaxia  | (x, y) original | radio (px) | reducción de flujo positivo |
|----------|------------------|------------|------------------------------|
| ngc3627  | (1403, 1531)     | 15         | 18.3% |
| ngc4594  | (1724, 852)      | 35         | 36.5% |
| ngc5055  | (1601, 1218)     | 45         | 44.8% |

Verificación post-máscara: el nuevo pico de cada galaxia baja a una razón mucho más razonable frente al percentil 99.9 (ngc5055: de 1431x a 8.1x, comparable a galaxias ya sanas como ngc3938 a 11x). En ngc4594 el nuevo máximo cae casi exactamente en el centro geométrico del recorte — consistente con el núcleo real de M104 (Sombrero), que es intrínsecamente muy compacto y brillante, no un artefacto. En ngc3627 queda un resto (92, era 231 antes de enmascarar) pero es un **único píxel aislado rodeado de ceros en las cuatro direcciones** — un rayo cósmico o píxel caliente puntual, de naturaleza distinta a la estrella mal sustraída (que ocupaba ~47 píxeles contiguos): tras la convolución con el beam se diluye igual que cualquier pico real de un solo píxel, así que no se le agregó máscara.

Esta máscara es parte del flujo normal de ejecución: como se aplica dentro de `SRF_MODEL_CASA.py`, `main.py` la ejecuta automáticamente para las 14 galaxias (para las 11 restantes, `STAR_MASKS.get(galaxia, [])` devuelve una lista vacía y no hace nada). No es una regla general — es específica a estas 3 galaxias por coordenada fija; si aparecieran artefactos similares en otro FITS de SINGS, habría que repetir el proceso (localizar el pico, verificar el perfil radial, agregar una entrada a `STAR_MASKS`).

### Galaxias sin causa identificada: ngc4536, ngc5194

Para estas dos no se encontró bug de signo ni artefacto puntual (suma positiva, sin outliers extremos, columna "Mapa radio" ya centrada y con estructura correcta tras el recorte). La ausencia de fuente visible en los paneles ngVLA/VLA en estos casos es consistente con el mismo fenómeno de sensibilidad interferométrica a brillo superficial extendido ya documentado más arriba (config A vs B), no con un defecto de código — a diferencia de M82 (núcleo compacto), estas galaxias reparten el mismo flujo total fijo (idéntico para las 14 galaxias a un z dado, por construcción de `_scale_sfr`) sobre un disco más extendido, resultando en menor brillo superficial por beam.

### Contornos de detección (3 sigma)

`plot_with_beam()` (duplicada en los ocho scripts de telescopio, `Modelo/` y `SINGS/`) generaba la figura por galaxia con `aplpy.FITSFigure` mostrando solo el mapa de color (`show_colorscale`) y el beam (`add_beam`). El problema: con mapas de color como `inferno`/`plasma`/`hot`, cualquier pico de ruido correlacionado cerca del percentil alto de la imagen se ve igual de "caliente" que una detección real — el color solo no permite distinguir señal de ruido.

**Fix**: se agregó una estimación de RMS por sigma-clipping (`astropy.stats.sigma_clipped_stats`, descarta iterativamente los píxeles que se alejan >3σ de la media para no dejar que la propia fuente sesgue el fondo hacia arriba) y se dibujan contornos en cian (`#00ffcc`) a 3σ, 5σ, 10σ y 20σ con `f.show_contour()`:

```python
data = fits.getdata(fits_file)
_, _, rms = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
f.show_contour(fits_file, levels=[3*rms, 5*rms, 10*rms, 20*rms],
               colors="#00ffcc", linewidths=1.2)
```

Verificado visualmente contra un `_noise.fits` ya generado: el contorno traza exactamente el borde de la región que supera el umbral estadístico, dejando sin marcar el resto del mapa de color aunque se vea "caliente" — que es el punto de agregarlo. Si no hay contorno visible en una combinación galaxia/z/config, es una no-detección o detección marginal real (ver sección de config A vs B más arriba), no un fallo del script.

Esto solo se aplicó a los paneles de telescopio (imágenes ya convolucionadas con el beam + ruido instrumental). El mapa de "SFR" pre-instrumento (`resultados/{z}/{galaxia}/mapa_*.fits`) no tiene un contorno equivalente: es ruido pixel a pixel sin suavizar y puede heredar gradiente de fondo del Hα de entrada, así que un contorno de silueta limpio ahí requeriría trabajo aparte de caracterización de fondo.

`reportes/dudas.tex` ya documentaba esta misma tarea (pedida antes por el asesor) con el razonamiento estadístico completo detrás de por qué 3σ, cómo cambia el número de contornos con el redshift, y un bug de ruido que se había encontrado y corregido en `add_correlated_noise` — pero el código de los contornos nunca había quedado aplicado en el repo; quedó implementado ahora siguiendo esa misma especificación.

**Bug encontrado al correr `main.py` con este cambio**: `f.show_contour()` de aplpy falla (`ValueError: Expected 2 world coordinates, got 0`) si alguno de los niveles pedidos no tiene ningún píxel por encima — el transform WCS de `astropy.visualization.wcsaxes` no admite un contorno vacío. Como el SNR varía mucho entre galaxia/z/config, no todas las imágenes alcanzan 10σ o 20σ. Fix: filtrar los niveles a los que la imagen realmente supera antes de llamar a `show_contour` (`[lvl for lvl in (3,5,10,20)*rms if data.max() > lvl]`), en los ocho scripts de telescopio.

**`collage.py`**: se agregó el mismo contorno (3/5/10/20σ, sigma-clipped, cian) a los cuatro paneles de telescopio del grid combinado, vía `ax.contour()` de matplotlib en vez de `aplpy.FITSFigure` — como esos paneles ya son `imshow` plano sin proyección WCS (no hay reproyección de por medio), el contorno se alinea exacto con la imagen sin necesitar aplpy, y se evita reescribir el motor de dibujo compuesto de `collage.py` (grid irregular con `GridSpec` anidado, colorbars y estilos manuales). El panel de "Mapa radio" sigue sin contorno, por la misma razón de arriba.

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
