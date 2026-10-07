# ANALYSIS STEPS

Data sit in /net/stratus/c2sm-data/CH2025/ogd-climate-scenarios-ch2025-grid. Orography: /net/stratus/c2sm-data/CH2025/ogd-climate-scenarios-ch2025-grid/meta

Copied two files to /highres/svenk/SEKEstorylines_snow/data_temp: ogd-climate-scenarios-ch2025-grid_ch_pr_smhi-rca-noresm_gwl3.0.nc --> pr.nc

## Topo: Maskout cells <1000 m
cdo -setctomiss,0 -gec,1000 topo.swiss1_ch01r.swiss.lv95.nc topo.swiss1_ch01r.swiss.lv95_gec1000.nc

## Mask temperature: All days<=2°C to 1, all others to missval
cdo lec,2 tas.nc tas_masked.nc

## 3 Multiply daily pr with daily temperature mask (0 --> days with zero precip or t>2°C, missval --> cells outside CH)
cdo mul pr.nc tas_masked.nc snow.nc

###

## 1 Compute monthly means
cdo monmean snow.nc snow_monmean.nc

## 2 Select Nov-Apr
cdo selmon,11,12,1,2,3,4 snow_monmean.nc snow_monmean_selmon.nc

## 3 Get rid of first four months (incomplete winter)
cdo delete,timestep=1/4 snow_masked_monmean_selmon.nc snow_masked_monmean_selmon_seltimestep.nc

## 4 Get rid of last two months (incomplete winter)
cdo seltimestep,1/-3 snow_masked_monmean_selmon_seltimestep.nc snow_masked_monmean_selmon_seltimestep_seltimestep.nc

## 5 Shift time by two months
cdo shifttime,+2months snow_masked_monmean_selmon_seltimestep_seltimestep.nc snow_masked_monmean_selmon_seltimestep_seltimestep_shifftime.nc

## 6 Annual mean
cdo yearmean snow_masked_monmean_selmon_seltimestep_seltimestep_shifttime.nc snow_masked_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean.nc



