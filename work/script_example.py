
#
# Here's an example script for running siracha.
#

# First, import libraries.
import rce_pce_column
import numpy as np

# Initialize the column object.
col = rce_pce_column.column( "intial_file_example.nc" )

# The model now could run at this point, but you 
# probably want to make some small changes now.
col.change_value('totaltime' , 300) # Simulation time to 300 days
col.change_value('recordtime',10)   # Output is mean of last 10 days
col.change_value('timestepsize',1/24)   # Timestep is 1 hour
col.temps_convection = False        # Let's turn off convection for
                                    # temperature, just 'cause.

col.extra['lyO_source'][1] = 'ppv/day'
col.extra['sink_ox'][1] = 'ppv/day'
col.extra['sink_nox'][1] = 'ppv/day'

# Run the model
col.wrapper()

# When done, write the output 
col.write_to_file( "output_file_example.nc" )
