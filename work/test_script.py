
# Basic running script

import rce_pce_column

import numpy as np
from netCDF4 import Dataset

infile = "initial_file_for_testing.nc"
col = rce_pce_column.column( infile )
col.change_value('totaltime' ,2)
col.change_value('recordtime',1)
col.change_value('timestepsize',1/24)
col.bool_ozone_integration = True
col.bool_run_rad = False

col.extra['lyO_source'][1] = 'ppv/day'
col.extra['sink_ox'][1] = 'ppv/day'
col.extra['sink_nox'][1] = 'ppv/day'

print( "Run the column with 'col.wrapper()'" )
print( "When finished, make an output file with "+
       "'write_to_file( <some filename> )'" )
print( "If you want to make an init file, which should be used to initialize "+
       "column objects, then set init_file=True in the write_to_file command." )
