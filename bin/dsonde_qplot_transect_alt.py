"""
Takes dropsonde netcdf file path and list of variables to plot and saves to png.  
If run without list of variables returns list of available variables 
If var 1 is 'all' will plot all variables in file


python .../Dropsonde_FAAM_quickplot.py <netcdf dropsonde filepath> var1, var2
"""
import sys
import iris
import iris.unit as unit

import numpy as np
from scipy.interpolate import griddata

import matplotlib.pyplot as plt

#import pdb

def calculate_distance_from_latlon(first_station_lon, first_station_lat, lon, lat):

    fslat_rad = np.radians(first_station_lat)
    fslon_rad = np.radians(first_station_lon)
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)

    #Haversine Formula
    
    a = np.sin((lat_rad-fslat_rad)/2)**2 + np.cos(lat_rad) * np.cos(fslat_rad) * np.sin((lon_rad-fslon_rad)/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    d = 6371 * c
    
    return d

def ConvertSecondsSince1970ToDatetime(hours_since):
    """ Takes numpy array of 'hours since 1970' and string format
     and returns numpy array of date strings"""

    u = unit.Unit('seconds since 1970-01-01 00:00:00',calendar='gregorian')

    return np.array([u.num2date(da) for da in np.array(hours_since)])

#Things to change, possibly

# File with plot parameters
pp = '/nfs/see-fs-01_users/eepdw/Github/MetCalcs/bin/plot_parameters.csv'
# Vertical param
altvar = 'altitude above MSL'
save_dir = ''

figprops = dict(figsize=(12,8), dpi=72)

#altvar = 'gps reported altitude above MSL'

plot_params = np.genfromtxt(pp, delimiter = ',', dtype='str') # var name, min, max

# Get the total number of args passed 

cmdargs = sys.argv[1:]
tot_args = len(cmdargs)

assert tot_args>0, "Need some filenames"

filepath = cmdargs[0]
cubes = iris.load(filepath)

lat_cubes = cubes.extract('north latitude')
lon_cubes = cubes.extract('east longitude')
alt_cubes = cubes.extract(altvar)

#pdb.set_trace()

if tot_args==1:
    print cubes
    #[print '%s, %s, %s' % (plot_cube.name(), \
      #                     plot_cube.data.min(), plot_cube.data.max()) for plot_cube in cubes]

assert tot_args>=2, "Number of input arguments wrong, should be 2"

# Get the arguments list

if cmdargs[1] == 'all':
    v_list = np.unique([str(cube.name()) for cube in cubes])

else:
    v_list = cmdargs[1:]

v_list = np.array(v_list)[np.in1d(v_list, plot_params)]

print v_list

for v in v_list:

    var_cubes = cubes.extract(v)

    ma_init=0

    for cube in var_cubes:

        #pdb.set_trace()

        t_idx = np.where([np.all(np.in1d(cube.coord('time').points, c.coord('time').points)) for c in lat_cubes])[0]

        if ma_init==0:
            ma_init=1
            lats=lat_cubes[t_idx].data
            lons=lon_cubes[t_idx].data
            alts=alt_cubes[t_idx].data
            var=cube.data
            time=cube.coord('time').points
            sounding_start = np.nanmin(cube.coord('time').points)

        else:
            ma_init=1
            lats = np.ma.append(lats, lat_cubes[t_idx].data)
            lons = np.ma.append(lons, lon_cubes[t_idx].data)
            alts = np.ma.append(alts, alt_cubes[t_idx].data)
            var = np.ma.append(var, cube.data)
            time = np.ma.append(time, cube.coord('time').points)
            sounding_start = np.ma.append(sounding_start, np.nanmin(cube.coord('time').points))

    start = np.argmin(time)

    distances = calculate_distance_from_latlon(lons[start], lats[start], lons, lats)

    if np.any(np.diff(distances<0.)):
        print 'WARNING: Distances, not always increasing, you probably want to use time script instead, as data gridding may get messy'
      
    xi, yi= np.meshgrid(np.linspace(distances.min(), distances.max(), 500),
                                        np.linspace(alts.min(), alts.max(), 1000))

    mask = distances.mask * alts.mask * var.mask

    grid_data = griddata((distances.data[~mask], alts.data[~mask]), var.data[~mask], (xi, yi), method='linear', fill_value=np.nan)

    # Plot

    p =np.where(plot_params[:,0]==v)[0]

    plt.figure(**figprops)

    plt.contourf(xi, yi, grid_data, np.linspace(float(plot_params[p,1][0]), float(plot_params[p,2][0]), 24))

    cbar = plt.colorbar()
    cbar.set_label(cube.units)

    plt.xlabel('Distance from first measurement / km')
    plt.ylabel(alt_cubes[t_idx].name()+' / '+alt_cubes[t_idx].units.symbol)
    plt.title(v)

    # Plot and label where soundings started

    sd_idx = np.where(np.in1d(time, sounding_start))[0]
    d_ss = distances[sd_idx]
    a_ss = alts[sd_idx]

    ss = ConvertSecondsSince1970ToDatetime(sounding_start)
    
    plt.scatter(d_ss, a_ss, color='k')
    for label, x, y in zip(ss, d_ss, a_ss):
        plt.annotate(
            label.strftime('%H:%M'), 
            xy = (x, y), xytext = (0, 1),
            textcoords = 'offset points', ha = 'right', va = 'bottom')

    if save_dir=='':
        plt.show()
    else:
        plt.savefig('%s' % (save_dir))

    

    

    
