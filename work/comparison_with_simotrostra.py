
# Basic running script

import typhon



import rce_pce_column

import numpy as np
from scipy.interpolate import interp1d
from netCDF4 import Dataset

infile = "initial_file_for_testing.nc"
col = rce_pce_column.column( infile )
col.change_value('totaltime' ,2)
col.change_value('recordtime',1)
col.change_value('timestepsize',1/24)
col.bool_ozone_integration = True
col.bool_run_rad = False

nc = Dataset( 'acefts_nox_trop_djf.nc', 'r' )
ace_nox = nc.variables['NOx'][:]
ace_lyz = nc.variables['z'][:]
nc.close()

#phlev = typhon.math.nlogspace(1000e2, 20, 195)  # [Pa]                  
#p = np.exp(0.5 * (np.log(phlev[1:]) + np.log(phlev[:-1])))
p = col.D['lyP']

t = typhon.physics.standard_atmosphere(p.data, coordinates='pressure')  # [K]
t = col.D['lyT']

z = typhon.physics.pressure2height(p, t)/1000  # [m]
zhlev = np.zeros(len(z)+1)
for i in range(0,len(z)): zhlev[i+1] = -zhlev[i] + 2*z[i]

nox = interp1d( ace_lyz.data, ace_nox.data, fill_value = 'extrapolate' )(z.data)
nox[np.isnan(nox)] = 0
nox[nox<0] = 0

initial_o3 = 3.6478 * (p*0.01)**0.83209 * np.exp(-(p*0.01) / 11.3515) * 1e-6

#col.change_value('lvP',phlev)
#col.change_value('lyP',p)
col.change_value('lyO',initial_o3)
#col.change_value('lyT',t)
#col.change_value('lyZ',z)
#col.change_value('lvZ',zhlev)

col.change_value('lyNOx',nox)

col.extra['lyO_source'][1] = 'ppv/day'
col.extra['sink_ox'][1] = 'ppv/day'
col.extra['sink_nox'][1] = 'ppv/day'

print( "Run the column with 'col.wrapper()'" )
print( "When finished, make an output file with "+
       "'write_to_file( <some filename> )'" )
print( "If you want to make an init file, which should be used to initialize "+
       "column objects, then set init_file=True in the write_to_file command." )
