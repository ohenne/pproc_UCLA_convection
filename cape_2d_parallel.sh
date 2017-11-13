#!/bin/bash
set -ex


#PFAD=`pwd`

FSTEM=test1plus4K #plus4K_detail
mkdir -p $modelo/$FSTEM/level2
dim=320
#INPFAD=/mnt/lustre01/scratch/m/m214106/${FSTEM}/
#OUTPFAD=/mnt/lustre01/scratch/m/m214106/${FSTEM}/paste_d8-11/


#cd $INPFAD

#for N in $(seq 0 ${dim}) ; do #$(seq 1 ${dim}); do
for N in $(seq 0 150) ; do
#for N in $(seq 4) ; do
  NSTRING=$(printf %01d $N)
  dimstr=$(printf %01d $dim)  #neccesary?
#  python calculate_cape_rivers_corrected_171104.py ${NSTRING} ${FSTEM} ${dimstr} &
  python calc_cape_jan_olga.py ${FSTEM} ${NSTRING} 1 &
#  read -p "Continue (y/n)?" choice
#  if [[ "$choice" == y ]]; then
#  echo "hm"
#  else
#  exit 1
#  fi
  #python calculate_cape_rivers_jan.py ${NSTRING}  &
done 


