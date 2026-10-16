#!/bin/ksh

# *** SEKE STORYLINES SNOW: CDO Version ***
# ***
# *** Compute seasonal (NOV-APRIL) mean snowfall sum for CH above 1000m. Temperature threshold: 2°C.

set -ex

# *** Directories ***
DATADIR=/net/stratus/c2sm-data/CH2025/ogd-climate-scenarios-ch2025-grid
METADIR=/highres/svenk/SEKEstorylines_snow/META
WORKDIR=/highres/svenk/SEKEstorylines_snow/WORK
RESDIR=/highres/svenk/SEKEstorylines_snow/RESULTS

# *** Filenames ***
FILE_BASE=ogd-climate-scenarios-ch2025-grid_ch
#FILE_TAS=tas.nc
#FILE_PR=pr.nc
FILE_ORO=topo.swiss1_ch01r.swiss.lv95.nc

# *** FURTHER SETTINGS
## GWL / PERIOD
PERIOD='gwl3.0'

# *** MODEL LIST ***
MODEL_LIST=$(cat "${METADIR}/model_list.txt")

cd $WORKDIR

## Topo: Maskout cells <1000 m
#cdo -setctomiss,0 -gec,1000 ${METADIR}/${FILE_ORO} ${METADIR}/${FILE_ORO%.nc}_gec1000.nc

for MODEL in $MODEL_LIST; do
    echo "Processing model: $MODEL"

## 1 Copy model files to working directory
cp ${DATADIR}/tas/${FILE_BASE}_tas_${MODEL}_${PERIOD}.nc ${WORKDIR}/tas.nc
cp ${DATADIR}/pr/${FILE_BASE}_pr_${MODEL}_${PERIOD}.nc ${WORKDIR}/pr.nc

## 2 Mask temperature: All days<=2°C to 1, all others to missval
cdo lec,2 tas.nc tas_masked.nc
rm -f tas.nc

## 3 Multiply daily pr with daily temperature mask (0 --> days with zero precip or t>2°C, missval --> cells outside CH)
cdo mul pr.nc tas_masked.nc snow.nc
rm -f pr.nc tas_masked.nc

## 4 Select Nov-Apr
cdo selmon,11,12,1,2,3,4 snow.nc snow_selmon.nc
rm -f snow.nc

## 5 get rid of first four months (incomplete winter --> timesteps 1-120) and last two months (incomplete winter --> timesteps 5370-5430)
cdo delete,timestep=1/120,5370/5430 snow_selmon.nc snow_selmon_seltimestep.nc
rm -f snow_selmon.nc

## 6 Shift time by two months and compute annual mean
cdo yearmean -shifttime,+2months snow_selmon_seltimestep.nc snow_selmon_seltimestep_yearmean.nc
rm -f snow_selmon_seltimestep.nc

## 7 Mask out cells <1000 m
## Note: use "ifthen" rather than "mul" here. With "mul", a cell with a
## true (non-missing) value of exactly 0 multiplied by a missing mask
## cell yields 0 instead of missing, so cells below 1000 m / outside CH
## with zero winter snowfall would incorrectly leak into the result.
## "ifthen" correctly propagates missingness from the mask regardless of
## the data value.
cdo ifthen ${METADIR}/${FILE_ORO%.nc}_gec1000.nc snow_selmon_seltimestep_yearmean.nc snow_selmon_seltimestep_yearmean_masked.nc

## 8 Compute spatial mean
cdo fldmean snow_selmon_seltimestep_yearmean_masked.nc snow_selmon_seltimestep_yearmean_masked_fldmean.nc

## 12 Output to textfile
cdo outputtab,year,value snow_selmon_seltimestep_yearmean_masked_fldmean.nc > snow_selmon_seltimestep_yearmean_masked_fldmean.txt

## 13 move relevant files to RESULT directory
mv snow_selmon_seltimestep_yearmean_masked_fldmean.txt  ${RESDIR}/result_${MODEL}_${PERIOD}.txt
mv snow_selmon_seltimestep_yearmean_masked_fldmean.nc  ${RESDIR}/result_${MODEL}_${PERIOD}.nc
mv snow_selmon_seltimestep_yearmean_masked.nc  ${RESDIR}/result_field_${MODEL}_${PERIOD}.nc
mv snow_selmon_seltimestep_yearmean.nc  ${RESDIR}/result_field_${MODEL}_${PERIOD}_nonmasked.nc

done

exit

