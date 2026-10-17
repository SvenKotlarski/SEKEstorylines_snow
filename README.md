# SEKE Storylines Snowfall Analysis

This directory contains scripts to calculate and rank a snowfall proxy for
Swiss climate-scenario data. The proxy is daily precipitation on days when
daily mean temperature is at or below 2 °C. For each winter year, the scripts
calculate the mean over November-April and then average spatially over valid
grid cells above 1000 m.

## Configuration and inputs

The scripts currently use absolute paths configured in the scripts:

- Daily tas and pr NetCDF files:
  `/net/stratus/c2sm-data/CH2025/ogd-climate-scenarios-ch2025-grid/{tas,pr}`
- Model list and topography files: `/highres/svenk/SEKEstorylines_snow/META`
- Intermediate CDO files: `/highres/svenk/SEKEstorylines_snow/WORK`
- Results: `/highres/svenk/SEKEstorylines_snow/RESULTS`
- Plots: `/highres/svenk/SEKEstorylines_snow/PLOTS`

The model list is `META/model_list.txt`. Both analysis scripts use period
`gwl3.0`. The topography mask
`META/topo.swiss1_ch01r.swiss.lv95_gec1000.nc` must already exist; the CDO
script contains a commented command showing how to create it from the
topography file.

## Analysis

Run one of the two analysis implementations:

```sh
ksh 01_analysis_cdo.sh
Rscript 01_analysis_R.R
```

Both implementations process daily data, select November-April, exclude the
incomplete leading and trailing winter periods, assign the remaining days to
winter years, and calculate daily-value-based annual means. They apply the
1000 m topography mask and calculate the spatial mean. The CDO implementation
uses `ifthen` for the topography mask so zero snowfall outside the mask remains
missing rather than entering the spatial mean.

The CDO outputs are written directly to `RESULTS`:

- `result_<model>_gwl3.0.txt` and `.nc`: spatially averaged time series
- `result_field_<model>_gwl3.0.nc`: annual gridded fields above 1000 m
- `result_field_<model>_gwl3.0_nonmasked.nc`: annual gridded fields without
  the elevation mask

The R implementation produces equivalent outputs with `_R` in the filename,
for example `result_<model>_gwl3.0_R.txt` and
`result_field_<model>_gwl3.0_nonmasked_R.nc`. If an output with the same name
already exists, it chooses the next available numbered suffix to avoid
overwriting it.

The R analysis requires the `ncdf4` package. The CDO analysis requires CDO
with support for the `delete,timestep` selection used by the script.

## Sorting

After running an analysis, create a tab-separated ranking of its results:

```sh
ksh 02_sort_snowfall.sh
ksh 02_sort_snowfall_R.sh
```

The regular sorter reads the CDO `result_*_gwl3.0.txt` files and writes
`RESULTS/snowfall_sorted.txt`. The R sorter reads the `_R.txt` results and
writes `RESULTS/snowfall_sorted_R.txt`. Both sort using the unrounded
snowfall values, then write the snowfall column rounded to two decimal places.
The output columns are `snowfall`, `year`, and `model`.

## Plotting

Run the ranking plot from this directory after running the regular sorter:

```sh
Rscript 03_plot_snowfall.R
```

It reads `RESULTS/snowfall_sorted.txt` and writes
`PLOTS/snowfall_sorted.png`, highlighting the five lowest-snowfall cases.

`04_plot_field_shaded_multimodel.R` plots the five lowest-snowfall winters and
the single highest-snowfall winter in `RESULTS/snowfall_sorted.txt`. Each
panel uses the matching
`result_field_<model>_gwl3.0_nonmasked.nc` timestep, so the map includes the
full spatial distribution, including elevations below 1000 m. The plot uses
topographic shading and a common snowfall colour scale:

```sh
Rscript 04_plot_field_shaded_multimodel.R [snowfall_sorted.txt] [out.png]
```

By default it reads `RESULTS/snowfall_sorted.txt` and writes
`PLOTS/snowfall_lowest5_unmasked.png`. It requires `ncdf4` and the topography
file `META/topo.swiss1_ch01r.swiss.lv95.nc`.

`04_plot_field_shaded_multimodel_3D.R` draws the same six winters as 3D
perspective maps (base R `persp`, no extra packages): the unmasked snowfall
field is draped over the terrain, with hillshade and 4x vertical exaggeration.
It writes `PLOTS/snowfall_lowest5_unmasked_3D.png` by default. Camera angle,
exaggeration and grid thinning are set by `VIEW_THETA`, `VIEW_PHI`,
`VERTICAL_EXAG` and `GRID_STEP` at the top of the script.
