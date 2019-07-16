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
f1  = scipy.io.netcdf_file("../urban-plume/out/urban_plume_process.nc",'r',mmap=False)


# get the variables
time       = f1.variables["time"].data / 3600 
dens1      = f1.variables["density"].data # unit: kg/m3
sul1       = f1.variables["tot_so4_conc"].data*29/96*1e9 / dens1 + f1.variables["tot_HSO4m_conc"].data*29/97*1e9 / dens1
nh4        = f1.variables["tot_nh4_conc"].data * 29/18 * 1e9/dens1
no3        = f1.variables["tot_no3_conc"].data * 29/62 * 1e9/dens1
oc         = f1.variables["tot_oc_conc"].data *  29/12 * 1e9/dens1
print(sul1)
# Input lwc value from the pycel model
#get the final value of lwc at 43.44 min, and we get the
#value from the lwc rate at around 0.05 g/kg per min
# plot
plt.ylim(0,20)
plt1 = axes.plot(time, sul1, "b-",   label   = "S(VI)", linewidth=1.5)
plt2 = axes.plot(time, nh4,  "g-",   label   = r"$\rm NH_4^+$", linewidth=1.5)
plt3 = axes.plot(time, no3,  "r-",   label   = r"$\rm NO_3^-$", linewidth=1.5)
plt4 = axes.plot(time, oc,   "k-",   label   = "OC", linewidth=1.5)
#plt2 = axes.plot(time, sul2, "r-",   label   = "Composition-averaged", linewidth=1.5)
axes.set_ylabel(r"Species (ppb)")
axes.set_xlabel(r"Time (h)")

#ax2 = axes.twinx()
#plt3 = ax2.plot(time, (sul2 - sul1)*100/sul1, "g-",   label   = "Difference", linewidth=1.5)
#ax2.set_ylabel(r"Diff ratio ($\%$)", color='g')
#ax2.tick_params('y', colors='g')

#major_ticks = np.arange(0, 10, 5)
#minor_ticks = np.arange(0, 25, 1)
#axes.set_xticks(major_ticks)
#axes.set_xticks(minor_ticks, minor=True)
#axes.tick_params('y', colors='b')
plt.grid(True)

#set legend
#plt4  = plt1+plt2 + plt3
#labs = [l.get_label() for l in plt4]
plt.legend(loc=2, frameon=False)
#axes.legend(plt4, labs, loc=9, frameon=False)

#figure.set_size_inches(4,5)
figure.savefig("../figs/urban_species.pdf")
