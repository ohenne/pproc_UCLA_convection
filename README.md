# pproc_UCLA_convection
post processing scripts for UCLA model output

prepare:
define $modelo in your basrc

-------------------------------
merge patchy modeloutput:
for 1d and 2d files (ts, ps)
./reduce_stat.sh

  uses:
  reducets.f90
  reduceps_geostr.f90

  compile on ocean with:
  gfortran reducets.f90 -o reducets `/groups/ocean/software/netcdf4/gcc/4.4.1.1/bin/nc-config --fflags --flibs`
  gfortran reduceps_geostr.f90 -o reduceps_geostr `/groups/ocean/software/netcdf4/gcc/4.4.1.1/bin/nc-config --fflags --flibs`
  creates:simname.ps.nc simname.ts.nc
-------------------------------

for 3d modeloutput:
./gather_3d_1_timeperiod_v2.sh
  check if cdo and nco is loaded
  creates files for every variable simname.out.vol.varname.nc in subdir paste from the output files
  can be done for selected range of timesteps
  selected timeranged can be pasted together with
connectfiles.sh

------------------------------

calculate variables from output:
3d_to_2d_olga_oct_v1.py
  calculates 2d variables: lwp, convergence(integrated), advection ...
  creates file:  test_budget.nc

calculate convective indices:
cape_2d_parallel.sh  (is $fstem used in program? need to be set again in py?)
  starts programm for single slices parallel
  set dimension (works only for xdim=ydim)

  uses:
  calculate_cape_rivers_jan.py >> simname.cape.river.ny.nc
  CAPE calculateion for single slices
  dimension need to be adapted; changed to be given in shell

merge files afterwards with:
  paste_rivers_v1.py >> simname.cape.river.nc

