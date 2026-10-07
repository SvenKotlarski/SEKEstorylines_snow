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

## 4 Compute monthly means
cdo monmean snow.nc snow_monmean.nc
rm -f snow.nc

## 5 Select Nov-Apr
cdo selmon,11,12,1,2,3,4 snow_monmean.nc snow_monmean_selmon.nc
rm -f snow_monmean.nc

## 6 Get rid of first four months (incomplete winter)
cdo delete,timestep=1/4 snow_monmean_selmon.nc snow_monmean_selmon_seltimestep.nc
rm -f snow_monmean_selmon.nc

## 7 Get rid of last two months (incomplete winter)
cdo seltimestep,1/-3 snow_monmean_selmon_seltimestep.nc snow_monmean_selmon_seltimestep_seltimestep.nc
rm -f snow_monmean_selmon_seltimestep.nc

## 8 Shift time by two months
cdo shifttime,+2months snow_monmean_selmon_seltimestep_seltimestep.nc snow_monmean_selmon_seltimestep_seltimestep_shifttime.nc
rm -f snow_monmean_selmon_seltimestep_seltimestep.nc

## 9 Annual mean
cdo yearmean snow_monmean_selmon_seltimestep_seltimestep_shifttime.nc snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean.nc
rm -f snow_monmean_selmon_seltimestep_seltimestep_shifttime.nc

## 10 Mask out cells <1000 m
cdo mul snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean.nc ${METADIR}/${FILE_ORO%.nc}_gec1000.nc snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean_masked.nc
rm -f snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean.nc

## 11 Compute spatial mean
cdo fldmean snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean_masked.nc snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean_masked_fldmean.nc

## 12 Output to textfile
cdo outputtab,year,value snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean_masked_fldmean.nc > snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean_masked_fldmean.txt

## 13 move relevant files to RESULT directory
mv snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean_masked_fldmean.txt  ${RESDIR}/result_${MODEL}_${PERIOD}.txt
mv snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean_masked_fldmean.nc  ${RESDIR}/result_${MODEL}_${PERIOD}.nc
mv snow_monmean_selmon_seltimestep_seltimestep_shifttime_yearmean_masked.nc  ${RESDIR}/result_field_${MODEL}_${PERIOD}.nc

done

exit

