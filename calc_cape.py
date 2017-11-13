import numpy
import matplotlib.pyplot as plt
import sys
import os #to use path from bashrc
from matplotlib import colors, cm
from matplotlib.ticker import MaxNLocator
from netCDF4 import Dataset
import os, errno
font = {'family' : 'serif',
	'color'  : 'black',
	'weight' : 'normal',
	'size'   : 12,
	}
#stem='test1'
stem=(sys.argv[1])
subdir='paste/'

step = 1
if len(sys.argv)>3:
    step=int(sys.argv[3])

print '---------------------------------------------------------'
print 'run name: \t',stem
print 'cross section at: y = ',int(sys.argv[2])
print 'step size in x direction: step = ',step, ' (default:1)'
print '---------------------------------------------------------'

# input path for horizontally averaged ps file
data_in_path_ps   = os.environ.get('modelo')+stem+'/level1/'

# input path for 4d volume files 
data_in_path      = os.environ.get('modelo')+stem+'/level1/'

# output path for final "river" file
data_out_path     = os.environ.get('modelo')+stem+'/level2/'

# defining the relevant 4d input filenames
data_in_t_name    = stem+'.out.vol.t.nc'
data_in_q_name    = stem+'.out.vol.q.nc'
data_in_l_name    = stem+'.out.vol.l.nc'

# defining the ps filename
data_in_ps_name   = stem+'.ps.nc'  # profile file for the pressure

# defining the cross section to be used
y_cross           = int(sys.argv[2])

# the final output filename
data_out_name     = stem+'.cape.river_y_'+str(y_cross).zfill(3)+'.nc'
out_filename      = data_out_path+data_out_name

# show plot or save plot?
l_show = False
l_save = True

# constants
Rd  = 287.04		# specific gas constant of dry air in J/(kg.K)
Rv  = 461.50		# specific gas constant of water vapor in J/(kg.K)
cpd = 1005.7		# specific heat capacity at const pressure of dry air in J/(kg.K)
cpv = 1870.0		# specific heat capacity at const pressure of water vapor in J/(kg.K)
g   = 9.81		# falling acceleration in m/(s.s)
p0  = 1.e5		# basis state pressure in Pa
Lv  = 2.5e6		# latent heat of vaporization at 0 deg C in J/kg
Gd  = 0.0098            # dry adiabatic lapse rate in K/m
eps = Rd/Rv             # ratio of gas constants (~0.62)


#------------------------------------------------------------------
# open netcdf dataset from UCLA-LES model: 
#------------------------------------------------------------------

data_in_ps   = Dataset(data_in_path_ps+data_in_ps_name, 'r')
data_in_z    = numpy.array(data_in_ps.variables["zt"][:])
data_in_zm   = numpy.array(data_in_ps.variables["zm"][:])

data_in_t    = Dataset(data_in_path+data_in_t_name, 'r')
data_in_q    = Dataset(data_in_path+data_in_q_name, 'r')
data_in_l    = Dataset(data_in_path+data_in_l_name, 'r')

data_in_time = data_in_t.variables['time']

# extracting the dimensions of the datafile
dim_x        = 1 #data_in_t.dimensions['xt'].size
dim_y        = data_in_t.dimensions['yt'].size
dim_z        = data_in_t.dimensions['zt'].size
dim_time     = data_in_t.dimensions['time'].size


# generate output file
try: 
    os.remove(out_filename)
except OSError:
    pass

data_out     = Dataset(out_filename, 'w', format='NETCDF4')
xt           = data_out.createDimension('xt',dim_x)
xt           = y_cross
yt           = data_out.createDimension('yt',dim_y)
time         = data_out.createDimension('time', None)


time_out     = data_out.createVariable('time','f8',('time',))
time_out.units   = 'seconds since 2010-01-01 00:00:00'

# create output variables
cape_out     = data_out.createVariable('cape', 'f4',('time', 'yt', 'xt'))
cin_out      = data_out.createVariable('cin',  'f4',('time', 'yt', 'xt'))
tot_cape_out = data_out.createVariable('total cape', 'f4',('time', 'yt', 'xt'))
tot_cin_out  = data_out.createVariable('total cin',  'f4',('time', 'yt', 'xt'))
scape_out    = data_out.createVariable('scape',  'f4',('time', 'yt', 'xt'))
lscape_out   = data_out.createVariable('lscape',  'f4',('time', 'yt', 'xt'))
lcl_out      = data_out.createVariable('lcl',  'f4',('time', 'yt', 'xt'))
lfc_out      = data_out.createVariable('lfc',  'f4',('time', 'yt', 'xt'))
loc_out      = data_out.createVariable('loc',  'f4',('time', 'yt', 'xt'))
ncross_out   = data_out.createVariable('ncross',  'f4',('time', 'yt', 'xt'))

# attributes long_name and unit
time_out.long_name = 'time'
time_out.units = 's'

cape_out.long_name = 'Convective available energy'
cape_out.units = 'J kg-1'

cin_out.long_name = 'Convective inhibition'
cin_out.units = ''

tot_cape_out.long_name = 'Total positive buoyant energy'
tot_cape_out.units = ''

tot_cin_out.long_name = 'Total negative buoyant energy'
tot_cin_out.units = ''

scape_out.long_name = 'Surface CAPE'
scape_out.units = ''

lscape_out.long_name = 'Height of positive buoyant surface air'
lscape_out.units = 'm'

lcl_out.long_name = 'Lifting condensation level'
lcl_out.units = ''

lfc_out.long_name = 'Level of free convection'
lfc_out.units = 'm'

loc_out.long_name = 'Limit of convection'
loc_out.units = 'm'

ncross_out.long_name = 'Number of crossings between background profile and test parcel'
ncross_out.units = '[-]'
# the air pressure is assumed to be relative weakly varying, therefore it is extracted from the horizontally averaged profile of pressure for only the first timestep
data_p       = numpy.array(data_in_ps.variables["p"][0,:])  # air pressure

# the exner function can be computed once for all horizontal coordinates
exner        = (data_p[:]/p0)**(Rd/cpd)
data_in_ps.close()

# the range of time steps and in index output_ti, initially set to zero
time_list = range(0,dim_time,step)
output_ti = 0

#for ti in range(dim_time):
for ti in time_list:

    time_out[ti]= data_in_time[ti]
    output_ti          = output_ti+1

# defining 2D vectors (the slices) for liquid water potential temperature thl, total water mixing ration qtot, and liquid water mixing ration ql.
    data_thl           = numpy.zeros(shape=(dim_y,dim_z))
    data_qtot          = numpy.zeros(shape=(dim_y,dim_z))
    data_ql            = numpy.zeros(shape=(dim_y,dim_z))

# filling the vectors with slices from the 3D volume files
    data_thl[:,:]      = numpy.array(data_in_t.variables["t"][ti,:,y_cross,:])
    data_qtot[:,:]     = numpy.array(data_in_q.variables["q"][ti,:,y_cross,:])
    data_ql[:,:]       = numpy.array(data_in_l.variables["l"][ti,:,y_cross,:])

# defining fields for the output variables: LCL, LFC, LOC, CAPE, and CIN
    z_lscape           = numpy.zeros(shape=(dim_y))
    z_lcl              = numpy.zeros(shape=(dim_y))
    z_free_conv        = numpy.zeros(shape=(dim_y))
    z_limit_of_conv    = numpy.zeros(shape=(dim_y))
    cape               = numpy.zeros(shape=(dim_y))
    cin                = numpy.zeros(shape=(dim_y))
    tot_cape           = numpy.zeros(shape=(dim_y))
    tot_cin            = numpy.zeros(shape=(dim_y))
    scape              = numpy.zeros(shape=(dim_y))
    nc                 = numpy.zeros(shape=(dim_y))

    nc[:]              = 0.
    scape[:]           = 0.
    z_lscape[:]        = 0.
    z_free_conv[:]     = numpy.nan
    z_limit_of_conv[:] = numpy.nan
    
    print(ti)

    for yi in range(0,dim_y):

	# potential temperature, temperature and virtual potential temperature 
	# from liquid water potential temperature
        tpot  = data_thl[yi,:] + Lv*data_ql[yi,:]/cpd/exner
        temp  = tpot*exner        
        # calculation of the virtual potential temperature, tvpot, note that this is approximate, as solid water components are neglected:
        # tvpot = tpot(1 + 0.61r_v-r_l), where r_l is liquid+solid. Here, only liquid is used.
        tvpot = tpot*(1.+(Rv/Rd-1.)*(data_qtot[yi,:]-data_ql[yi,:])-data_ql[yi,:])

	# water vapor mixing ratio and specific humidity 
	# in lowest model level from total water mixing ratio
        rv    = data_qtot[yi,1] - data_ql[yi,1]
        shum  = rv/(1.-rv)
        
	# find pressure at lifted condensation level (lcl)
        pmin = 0.5*p0
        pmax = p0
        for iter in range(10):
            # defining a 'test pressure'
            plcl = (pmin+pmax)*0.5

            # computing a 'test temperature', i.e. the temperature at the 'test pressure':            
            t_lcl = temp[1]*(plcl/p0)**(Rd/cpd)

	    # saturation vapor pressure for the 'test temperature':
            esat = 610.78*numpy.exp(17.2694*(t_lcl-273.16)/(t_lcl-35.85))

	    # saturation mixing ratio at lcl:
            rlcl = eps * esat/(plcl-esat)

            if rlcl>rv:   # air undersaturated
                pmax = plcl
            else:	  # air oversaturated
                pmin = plcl

        #### end of iteration: tlcl and plcl are now the temperature and pressure at the lifting condensation level. ####

        # now calculate virtual potential temperature profile of air parcel
        # make empty arrays for tvpot_parcel and temp_parcel
        tvpot_parcel   = numpy.zeros_like(tvpot)
        temp_parcel    = numpy.zeros_like(tvpot)
        temp_parcel[1] = temp[1]
        tpot_parcel    = temp[1]/exner[1]

        # it is assumed that the first level is below the lifting condensation level
        below_lcl = True
        # cycling through the vertical levels
        for zlev in range(1,dim_z):
            if data_p[zlev] > plcl: # still below lcl? --> dry adiabatic lifting
                tvpot_parcel[zlev] = tpot_parcel*(1.+(Rv/Rd-1.)*rv)
                temp_parcel[zlev]  = tpot_parcel*exner[zlev]	    	
            else:
                if below_lcl:		    # first level above lcl?
                    below_lcl = False

        	    # linear interpolation of lifting condensation z-level:
                    # note that this is somewhat coarse, since p does not decrease linearly with height.
                    epsilon   = (data_p[zlev-1]-plcl)/(data_p[zlev-1]-data_p[zlev])
                    z_lcl[yi] = data_in_z[zlev-1] + epsilon*(data_in_z[zlev]-data_in_z[zlev-1])
                    p_low     = plcl
                    temp_low  = t_lcl
                    z_low     = z_lcl[yi]
                else: # i.e., no longer below the lcl. We now need the pseudoadiabatic temperature profile
                    p_low     = data_p[zlev-1]
                    temp_low  = temp_parcel[zlev-1]
                    z_low     = data_in_z[zlev-1]

		# defining saturation vapor pressure (esat), saturation mixing ratio (rsat) and the derivative d esat(T)/dT.
                esat = 610.78*numpy.exp(17.2694*(temp_low-273.16)/(temp_low-35.85))
                rsat = eps  *esat/(p_low-esat)
                dedT = 4098.*esat/(temp_low-35.85)**2

		# T_lapse = (Gd + ((Lv/cpd)*(Rd/Rv)*esat*g)/(p_low*Rd*temp_low)) / (1.+((Lv/cpd)*(Rd/Rv)*dedT)/p_low)
                # here, Rm is approximated by Rd
                # the derivation can be found David Tarboton lecture notes or in K. Emanual, Atm. Conv, p. 131           
                T_lapse = (Gd + ((Lv/cpd)*eps*esat*p_low*g)/((p_low-0.378*esat)**2.*Rd*temp_low)
	        	   ) / (1.+(Lv/cpd)*(eps/(p_low-0.378*esat)+0.378*eps*esat/(p_low-0.378*esat)**2)*dedT)

                temp_parcel[zlev]  = temp_low - T_lapse*(data_in_z[zlev]-z_low)
                tvpot_parcel[zlev] = temp_parcel[zlev]/exner[zlev]*(1.+(Rv/Rd-1.)*rsat)

        tvpot_parcel[0] = tvpot_parcel[1]

        ### end of loop on zlev: temp_parcel and tvpot_parcel are now determined. #####

	### integration of CAPE and CIN: this requires a final loop over zlev     #####
        below_free_convection     = True
        below_limit_of_convection = True
        Bhelp = (tvpot_parcel[2:dim_z]-tvpot[2:dim_z])*(tvpot_parcel[1:dim_z-1]-tvpot[1:dim_z-1])
        indexes = [index for index in range(len(Bhelp)) if Bhelp[index] < 0] 
        nc[yi]=0 if all(indexes)==0 else len(indexes) 
        surface_cape = True if nc[yi]>=3 else False
        B_old = 0.;        
        for zlev in range(1,dim_z-1):

            # compute the buoyancy of the test parcel relative to the environmental air at the same zlevel
            B = g/tvpot[zlev]*(tvpot_parcel[zlev]-tvpot[zlev])*(data_in_zm[zlev]-data_in_zm[zlev-1])            
            
            if B>0:
                tot_cape[yi] = tot_cape[yi] + B
            else:
                tot_cin[yi]  = tot_cin[yi]  - B 
            

            #if B*B_old<0.: # testing for a change of sign in buoyancy
             #   nc[yi] = nc[yi] + 1.
            
            B_old = B
            
            # computing surface cape. surface cape is defined as all positive buoyancy starting at the surface, until buoyancy becomes negative 
            # for the first time
            if surface_cape:
                if B>=0:
                    scape[yi] = scape[yi] + B
                else: # B < 0
                    surface_cape = False
                    # linear interpolation of z_scape:
                    epsilon         = (tvpot_parcel[zlev]-tvpot[zlev])/((tvpot_parcel[zlev+1]-tvpot[zlev+1])-(tvpot_parcel[zlev]-tvpot[zlev]))
                    z_lscape[yi]    = data_in_z[zlev] + epsilon*(data_in_z[zlev+1]-data_in_z[zlev])
    
            elif below_free_convection:
	    	# CIN integration:
                cin[yi] = cin[yi] - B
                
	    	# check if next level is above lcl AND above level of free convection
                if data_in_z[zlev+1]>=z_lcl[yi] and tvpot[zlev+1]<tvpot_parcel[zlev+1]:
                    below_free_convection = False

        	    # linear interpolation of z_free_conv:
                    epsilon         = (tvpot[zlev]-tvpot_parcel[zlev])/((tvpot_parcel[zlev+1]-tvpot[zlev+1])-(tvpot_parcel[zlev]-tvpot[zlev]))
                    z_free_conv[yi] = data_in_z[zlev] + epsilon*(data_in_z[zlev+1]-data_in_z[zlev])

            elif below_limit_of_convection:
	    	# CAPE integration:
                cape[yi] = cape[yi] + B

            	# check if next level is above limit of free convection
                if tvpot[zlev+1]>tvpot_parcel[zlev+1]:
                    below_limit_of_convection = False

        	    # linear interpolation of z_limit_of_conv:
                    epsilon             = (tvpot_parcel[zlev]-tvpot[zlev])/((tvpot[zlev+1]-tvpot_parcel[zlev+1])-(tvpot[zlev]-tvpot_parcel[zlev]))
                    z_limit_of_conv[yi] = data_in_z[zlev] + epsilon*(data_in_z[zlev+1]-data_in_z[zlev])

        #if below_free_convection: cin[yi] = numpy.nan  # no buoyant layer
    numpy.where(nc > 5, tot_cape-scape-cape, cape)
    scape_out[ti,:]    = scape
    cape_out[ti,:]     = cape
    cin_out[ti,:]      = cin
    tot_cape_out[ti,:] = tot_cape
    tot_cin_out[ti,:]  = tot_cin
    lscape_out[ti,:]   = z_lscape
    lcl_out[ti,:]      = z_lcl
    lfc_out[ti,:]      = z_free_conv
    loc_out[ti,:]      = z_limit_of_conv
    ncross_out[ti,:]   = nc

data_out.close()

#fig = plt.figure(figsize=(12, 12))

#a = fig.add_subplot(2,2,1, adjustable='box', aspect=1.0)
#plt.contourf(data_in_yt,data_in_xt,cape,extend='both',cmap=cm.RdBu)

#a = fig.add_subplot(2,2,2, adjustable='box', aspect=1.0)
#plt.contourf(data_in_yt,data_in_xt,cin,extend='both',cmap=cm.RdBu)

#a = fig.add_subplot(2,2,3, adjustable='box', aspect=1.0)
#plt.contourf(data_in_yt,data_in_xt,z_free_conv,extend='both',cmap=cm.RdBu)

#a = fig.add_subplot(2,2,4, adjustable='box', aspect=1.0)
#plt.contourf(data_in_yt,data_in_xt,z_limit_of_conv,extend='both',cmap=cm.RdBu)

#if l_show: plt.show()
#if l_save: plt.savefig('CAPE.river.png')
#plt.close()

data_in_t.close()
data_in_q.close()
data_in_l.close()
