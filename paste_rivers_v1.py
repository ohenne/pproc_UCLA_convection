import numpy
import matplotlib.pyplot as plt
import sys
from matplotlib import colors, cm
from matplotlib.ticker import MaxNLocator
from netCDF4 import Dataset
import os #to use path from bashrc
	
stem='test1plus4K' #_detail'
subdir='level2/'


n=5
data_in_path_ps   = os.environ.get('modelo')+'/'+stem+'/'
data_in_path      = data_in_path_ps+subdir+'/'

#data_out_name = stem+'.cape.river_y_'+sys.argv[1]+'.nc'
data_out_name = stem+'.cape.river.nc'

dim_x = 320
dim_y = 320


#dim_x = data_in_t.variables['xt'].size
#dim_y = data_in_t.variables['yt'].size
#dim_z = len(data_in_z)

# generate output file
print(data_in_path+data_out_name)
data_out = Dataset(data_in_path+data_out_name, 'w', format='NETCDF4')
xt = data_out.createDimension('xt',dim_x)
yt = data_out.createDimension('yt',dim_y)
time = data_out.createDimension('time', None)

time_out = data_out.createVariable('time','f8',('time',))
#xt_out = data_out.createVariable('xt','f8',('xt',))
#yt_out = data_out.createVariable('yt','f8',('yt',))

#xt_out[:] = data_in_xt
#yt_out[:] = data_in_yt

# create output variables
cape_out  = data_out.createVariable('cape', 'f4',('time', 'xt',  'yt'))
cin_out   = data_out.createVariable('cin',  'f4',('time', 'xt',  'yt'))
lcl_out   = data_out.createVariable('lcl',  'f4',('time', 'xt',  'yt'))
lfc_out   = data_out.createVariable('lfc',  'f4',('time', 'xt',  'yt'))
loc_out   = data_out.createVariable('loc',  'f4',('time', 'xt',  'yt'))

#data_p    = numpy.array(data_in_ps.variables["p"][0,:])  # air pressure
#data_in_ps.close()


#for ti in range(dim_time):
#for yi in range(0,dim_y):
for yi in range(0,(dim_y-1)):
    print(yi)
    data_in   = stem+'.cape.river_y_'+str(yi)+'.nc'
    print(data_in)
    print(data_in_path)
    
    data_in_t = Dataset(data_in_path+data_in, 'r')


    if yi==0: 
        data_in_time = numpy.array(data_in_t.variables["time"][:])
    
    in_cape            = numpy.array(data_in_t.variables["cape"][:])
    in_cin             = numpy.array(data_in_t.variables["cin"][:])
    in_lcl             = numpy.array(data_in_t.variables["lcl"][:])
    in_lfc             = numpy.array(data_in_t.variables["lfc"][:])
    in_loc             = numpy.array(data_in_t.variables["loc"][:])

    dim_time = len(data_in_time)
    time_list = range(dim_time)
    output_ti = 0

#    for ti in time_list:

    #time_out[output_ti] = data_in_time[ti]
    #output_ti = output_ti+1

 
        #lcl           = numpy.zeros(shape=(dim_y))
        #z_free_conv     = numpy.zeros(shape=(dim_y))
        #z_limit_of_conv = numpy.zeros(shape=(dim_y))
        #cape            = numpy.zeros(shape=(dim_y))
        #cin             = numpy.zeros(shape=(dim_y))

        #cape=data_in_cape[ti,:]
        #cin=data_in_cin[ti,:]
        #cape=data_in_cape[ti,:]
        #cape=data_in_cape[ti,:]
        #cape=data_in_cape[ti,:]

        #z_free_conv[:]     = numpy.nan
        #z_limit_of_conv[:] = numpy.nan
    
        #print(ti)

    cape_out[:,:,yi] = in_cape
    cin_out [:,:,yi] = in_cin
    lcl_out [:,:,yi] = in_lcl
    lfc_out [:,:,yi] = in_lfc
    loc_out [:,:,yi] = in_loc
        
data_out.close()

data_in_t.close()
#data_in_q.close()
#data_in_l.close()
