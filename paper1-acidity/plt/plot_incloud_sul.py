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
f2  = scipy.io.netcdf_file("../urban-plume/out_5h_avcomp_onebin/urban_plume_aq_chem_b_process.nc",'r',mmap=False)


# get the variables
time       = f1.variables["time"].data / 60 
dens1      = f1.variables["density"].data # unit: kg/m3
dens2      = f2.variables["density"].data # unit: kg/m3
sul1       = f1.variables["tot_so4_conc"].data*29/96*1e9 / dens1 + f1.variables["tot_HSO4m_conc"].data*29/97*1e9 / dens1
sul2       = f2.variables["tot_so4_conc"].data*29/96*1e9 / dens2 + f2.variables["tot_HSO4m_conc"].data*29/97*1e9 / dens2

# Input lwc value from the pycel model
#get the final value of lwc at 43.44 min, and we get the
#value from the lwc rate at around 0.05 g/kg per min
# plot
plt1 = axes.plot(time, sul1, "b-",   label   = "Particle-resolved", linewidth=1.5)
plt2 = axes.plot(time, sul2, "r--",   label   = "Composition-averaged", linewidth=1.5)
axes.set_ylabel(r"Sulfate (ppb)")
axes.set_xlabel(r"Time (min)")

ax2 = axes.twinx()
plt3 = ax2.plot(time, (sul2 - sul1)*100/sul1, "g-",   label   = "Difference", linewidth=1.5)
ax2.set_ylabel(r"Diff ratio ($\%$)", color='g')
ax2.tick_params('y', colors='g')

#major_ticks = np.arange(0, 10, 5)
#minor_ticks = np.arange(0, 25, 1)
#axes.set_xticks(major_ticks)
#axes.set_xticks(minor_ticks, minor=True)
#axes.tick_params('y', colors='b')
plt.grid(True)

#set legend
plt4  = plt1+plt2 + plt3
labs = [l.get_label() for l in plt4]
#plt.legend(loc=9, frameon=False)
axes.legend(plt4, labs, loc=9, frameon=False)
#plt.grid(True, which='both')

#major ticks every 20, minor ticks every 5
#major_ticks = np.arange(0, 3600, 200)
#minor_ticks = np.arange(0, 3600, 100)
#axes.set_yticks(major_ticks)
#axes.set_yticks(minor_ticks, minor=True)

# set labels

#figure.set_size_inches(4,5)
figure.savefig("../figs/bulk_sul_compare_5h_onebin.pdf")
