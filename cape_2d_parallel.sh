#!/bin/bash
set -ex


#PFAD=`pwd`

FSTEM=test1plus4K #plus4K_detail
ystart=0
yend=100
mkdir -p $modelo/$FSTEM/level2


for N in $(seq $ystart $yend) ; do
  NSTRING=$(printf %01d $N)
  dimstr=$(printf %01d $dim)  #neccesary?
  python calc_cape.py ${FSTEM} ${NSTRING} 1 &
#  read -p "Continue (y/n)?" choice
#  if [[ "$choice" == y ]]; then
#  echo "hm"
#  else
#  exit 1
#  fi
  #python calculate_cape_rivers_jan.py ${NSTRING}  &
done 
wait 
cdo gather test1plus4K.cape.river_y_???.nc test1plus4K.cape.river_${ystart}_${yend}.nc 

