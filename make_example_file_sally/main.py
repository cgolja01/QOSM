
import numpy as np
from netCDF4 import Dataset

nc = Dataset( 'ozone_model_comparison_example_raw.nc' )
data = {
        k:{ 'vals':np.squeeze(nc.variables[k])
           ,'units':nc.variables[k].units
          } for k in nc.variables.keys()
       }
nc.close()

sza = np.array([143.80, 135.76, 124.93, 112.64, 99.580, 86.097,\
                72.389, 58.596, 44.869, 31.491, 19.395, 13.058,\
                19.391, 31.487, 44.866, 58.595, 72.391, 86.103,\
                99.591, 112.66, 124.95, 135.81, 143.86, 146.98])

nc = Dataset( 'ozone_model_comparison_example.nc', 'w' )

nc.createDimension( 'timestep', len(sza) )
nc.createVariable('sza', sza.dtype, ('timestep',))
nc.variables['sza'][:] = sza
nc.variables['sza'].units = 'Degrees, 0 is overhead, 90 horizon.'

nc.createDimension( 'layer', len(data['lyZ']['vals'] ) )
nc.createDimension( 'level', len(data['lvZ']['vals'] ) )

col_name_orig = 0
col_name_outp = 1
col_long_name = 2
keys = [('lyZ','z','altitude')
       ,('lyP','p','pressure')
       ,('lyT','T','temperature')
       ,('lyO','o3','ozone')
       ,('nm','nm','Total number density')
       ,('lyO_source','o3_s','Ozone source')
       ,('sink_nox','o3_l_nox','Ozone loss du to NOx')
       ,('sink_ox','o3_l_ox','Ozone loss du to Ox')
       ,('j2','j2','j2 (source) photolysis rate')
       ,('j2_srb','j2_srb','j2 (source) photolysis rate from Schumann-Runge')
       ,('j2_hz','j2_hz','j2 (source) photolysis rate from Herzberg')
       ,('j3','j3','j3 (sink) photolysis rate)')
       ,('lyNOx','nox','Odd-nitrogen family mixing ratio')]
keys.sort()

for k in keys:
    nc.createVariable( k[col_name_outp]
                     , data[k[col_name_orig]]['vals'].dtype
                     , ('layer') )
    nc.variables[k[col_name_outp]][:] = data[k[col_name_orig]]['vals']
    nc.variables[k[col_name_outp]].units = data[k[col_name_orig]]['units']
    nc.variables[k[col_name_outp]].long_name = k[col_long_name]

keys = [('lvP','p_lev','pressure on levels')]

for k in keys:
    nc.createVariable( k[col_name_outp]
                     , data[k[col_name_orig]]['vals'].dtype
                     , ('level') )
    nc.variables[k[col_name_outp]][:] = data[k[col_name_orig]]['vals']
    nc.variables[k[col_name_outp]].units = data[k[col_name_orig]]['units']
    nc.variables[k[col_name_outp]].long_name = k[col_long_name]



nc.close()

