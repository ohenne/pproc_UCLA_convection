#!/bin/bash
set -ex

RUNNAME=plus4K_detail

PFAD=`pwd`

OUTPATH=${modelo}/${RUNNAME}/level1/

echo $PFAD

# path to the output data
cd ${modelo}/${RUNNAME}/level0/


nn=$( ls -l ${RUNNAME}.ts.????????.nc | wc -l )
nx=$( ls -l ${RUNNAME}.ts.0000????.nc | wc -l )
ny=$( ls -l ${RUNNAME}.ts.????0000.nc | wc -l )

${PFAD}/reducets ${RUNNAME} $nx $ny ${OUTPATH}
${PFAD}/reduceps_geostr ${RUNNAME} $nx $ny ${OUTPATH}

exit
