"""
Takes dropsonde netcdf file path and list of variables to plot and saves to png.  
If run without list of variables returns list of available variables 
If var 1 is 'all' will plot all variables in file


python .../Dropsonde_FAAM_quickplot.py <netcdf dropsonde filepath> var1, var2
"""
import sys
import iris
import matplotlib.pyplot as plt
import numpy as np

import pdb

# File with plot parameters

pp = ''
plot_params = np.readfromtxt(pp, delimiter = ',') # var name, min, max

# Get the total number of args passed 
pdb.set_trace()

tot_args = len(sys.argv)
cmdargs = str(sys.argv)
filepath = cmdargs[0]

cubes = iris.load(filepath)

if tot_args<2:
    print cubes

assert tot_args==2, "Number of input arguments wrong, should be 2"

# Get the arguments list

if cmdargs[1] != 'all':
    v_list = [cmdargs[1:]]

else:
    v_list = [str(cube.name()) for cube in cubes]

for v in v_list:

    plot_cube = cubes.extract(v)

    

    

    
