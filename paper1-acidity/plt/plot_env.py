#!/usr/bin/env python
import sys, os
sys.path.append("../../../aq_partmc/tool/")
import mpl_helper
import scipy.io, numpy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pylab
import csv
(figure, axes) = mpl_helper.make_fig(right_margin=0.8)
#axes2 = axes.twinx()

# input the file
#f1  = scipy.io.netcdf_file("out_right_low/cloud_parcel_process.nc",'r',mmap=False)
f1  = scipy.io.netcdf_file("../urban-plume/out_5h/urban_plume_aq_chem_b_process.nc",'r',mmap=False)
f2  = scipy.io.netcdf_file("../urban-plume/out_5h_avcomp/urban_plume_aq_chem_b_process.nc",'r',mmap=False)


# get the variables
time       = f2.variables["time"].data / 60
h2o1       = f2.variables["tot_h2o_conc"].data*1000 # liquid water content, initial unit: kg m-3
temp1      = f2.variables["temp"].data
dens1      = f2.variables["density"].data # unit: kg/m3
lwc1       = h2o1/dens1            # g/kg
rh1        = f2.variables["rh"].data * 100

print(time)
print(lwc1)
print(rh1)

# Input lwc value from the pycel model
#get the final value of lwc at 43.44 min, and we get the
#value from the lwc rate at around 0.05 g/kg per min
# plot
plt1 = axes.plot(time, rh1, "b-",   label   = "RH", linewidth=1.5)
axes.set_ylabel(r"RH", color='b')
#major_ticks = np.arange(0, 25, 5)
#minor_ticks = np.arange(0, 25, 1)
#axes.set_xticks(major_ticks)
#axes.set_xticks(minor_ticks, minor=True)
axes.tick_params('y', colors='b')
plt.grid(True)
ax2 = axes.twinx()
plt2 = ax2.plot(time, temp1, "g-",   label   = "Temperature", linewidth=1.5)
ax2.set_ylabel(r"T(K)", color='g')
ax2.tick_params('y', colors='g')

#set legend
plt4  = plt1+plt2 
labs = [l.get_label() for l in plt4]
axes.legend(plt4, labs, loc=7, frameon=False)
#plt.grid(True, which='both')

#major ticks every 20, minor ticks every 5
#major_ticks = np.arange(0, 3600, 200)
#minor_ticks = np.arange(0, 3600, 100)
#axes.set_yticks(major_ticks)
#axes.set_yticks(minor_ticks, minor=True)

# set labels
axes.set_xlabel(r"Time (min)")

#figure.set_size_inches(4,5)
figure.savefig("../figs/init_env_avcom.pdf")
