"""
Takes dropsonde netcdf file path and list of variables to plot
and saves to png.  
If run without list of variables returns list of available variables 

python .../Dropsonde_FAAM_quickplot.py <netcdf dropsonde filepath> [var1, var2]

"""

import sys
import iris
import matplotlib.pyplot

import pdb

# Get the total number of args passed 
pdb.set_trace()

cmdargs = sys.argv[1:]
tot_args = len(cmdargs)

filepath = cmdargs[0]
cubes = iris.load(filepath)

if tot_args<2:
    print cubes

assert tot_args>=2, "Number of input arguments should be >2, if you want to plot anythingg"

# Get the arguments list

v_list = [cmdargs[1:]]

for v in v_list:

    plot_cube = cubes.extract(v)
    
    if plot_cube!=None:
        print plot_cube
        
    else:
        print 'cube %s does not exist in this file' % v

    pdb.set_trace()

# plot all
# plot section


    
