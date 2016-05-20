"""
Takes dropsonde netcdf file path and list of variables to plot and saves to png.  
If run without list of variables returns list of available variables 
If var 1 is 'all' will plot all variables in file

python .../Dropsonde_FAAM_quickplot.py <netcdf dropsonde filepath> var1, var2
"""xs
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

<<<<<<< HEAD
cmdargs = sys.argv[1:]
tot_args = len(cmdargs)

=======
tot_args = len(sys.argv)
cmdargs = str(sys.argv)
>>>>>>> refs/remotes/origin/master
filepath = cmdargs[0]

cubes = iris.load(filepath)

if tot_args<2:
    print cubes

assert tot_args>=2, "Number of input arguments should be >2, if you want to plot anythingg"

# Get the arguments list

<<<<<<< HEAD
v_list = [cmdargs[1:]]
=======
if cmdargs[1] != 'all':
    v_list = [cmdargs[1:]]

else:
    v_list = [str(cube.name()) for cube in cubes]
>>>>>>> refs/remotes/origin/master

for v in v_list:

    plot_cube = cubes.extract(v)
    
    if plot_cube!=None:
        print plot_cube
        
    else:
        print 'cube %s does not exist in this file' % v

    

    

# plot all
# plot section


    
