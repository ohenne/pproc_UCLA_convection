#!/bin/bash
set -ex


#PFAD=`pwd`

FSTEM=plus4K_detail 
#ystart=0
interval=79
yend=0
#yend=2
mkdir -p $modelo/$FSTEM/level2

for ystart in 0 80 160 240; do
 yend=$(($ystart + $interval))
 for N in $(seq $ystart $yend) ; do
   NSTRING=$(printf %01d $N)
   dimstr=$(printf %01d $dim)  #neccesary?
   python calc_cape.py ${FSTEM} ${NSTRING} 1 &
#   read -p "Continue (y/n)?" choice
#   if [[ "$choice" == y ]]; then
#   echo "hm"
#   else
#   exit 1
#   fi
  #python calculate_cape_rivers_jan.py ${NSTRING}  &
 done 
 wait
done 
cd $modelo/$FSTEM/level2
cdo gather $modelo/$FSTEM/level2/${FSTEM}.cape.river_y_???.nc $modelo/$FSTEM/level2/${FSTEM}.cape.river.nc #_${ystart}_${yend}.nc 
ncks -v $modelo/$FSTEM/level2/${FSTEM}.cape.river_y_000.nc $modelo/$FSTEM/level2/time.nc
ncks -A $modelo/$FSTEM/level2/time.nc $modelo/$FSTEM/level2/${FSTEM}.cape.river.nc #_${ystart}_${yend}.nc
rm $modelo/$FSTEM/level2/time.nc
#rm ${FSTEM}.cape.river_y_???.nc
