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
#f1  = scipy.io.netcdf_file("../urban-plume/out_5h/urban_plume_aq_chem_b_process.nc",'r',mmap=False)
#f2  = scipy.io.netcdf_file("../urban-plume/out_5h/gas_state.nc",'r',mmap=False)
f1  = scipy.io.netcdf_file("../urban-plume/out_5h_avcomp/urban_plume_aq_chem_b_process.nc",'r',mmap=False)
f2  = scipy.io.netcdf_file("../urban-plume/out_5h_avcomp/gas_state.nc",'r',mmap=False)

# get the variables
time       = f1.variables["time"].data / 60
dens1      = f1.variables["density"].data # unit: kg/m3
ao3        = f1.variables["tot_aO3_conc"].data * 29/48 * 1e9/dens1
ah2o2      = f1.variables["tot_aH2O2_conc"].data * 29/34 * 1e9/dens1
s_iv       = f1.variables["tot_aSO2_conc"].data * 29/64 * 1e9/dens1 + f1.variables["tot_HSO3m_conc"].data * 29/81 * 1e9/dens1

gas = f2.variables["gas_mixing_ratio"].data
O3g  = gas[:,10]
H2O2g= gas[:,15]
SO2g = gas[:,17]
print(s_iv)
print(O3g)
print(H2O2g)
print(SO2g)
# Input lwc value from the pycel model
#get the final value of lwc at 43.44 min, and we get the
#value from the lwc rate at around 0.05 g/kg per min
# plot
#plt.ylim(0,20)
#plt1 = axes.plot(time, sul1, "b-",   label   = "S(VI)", linewidth=1.5)
plt1 = axes.plot(time, ao3,    "g-",   label   = r"$\rm aO_3$", linewidth=1.5)
plt2 = axes.plot(time, ah2o2,  "r-",   label   = r"$\rm aH_2O_2$", linewidth=1.5)
plt3 = axes.plot(time, s_iv,  "b-",   label   = r"$\rm aSO_2+HSO_3^-$", linewidth=1.5)
#plt2 = axes.plot(time, sul2, "r-",   label   = "Composition-averaged", linewidth=1.5)
axes.set_ylabel(r"Aqueous species (ppb)")
axes.set_xlabel(r"Time (min)")

ax2  = axes.twinx()
plt4 = ax2.plot(time, O3g,   "g--",   label   = r"Gas $\rm O_3$",    linewidth=1.5)
plt5 = ax2.plot(time, H2O2g, "r--",   label   = r"Gas $\rm H_2O_2$", linewidth=1.5)
plt6 = ax2.plot(time, SO2g,  "b--",   label   = r"Gas $\rm SO_2$", linewidth=1.5)
ax2.set_ylabel(r"Gas species (ppb)")
#ax2.set_ylabel(r"Diff ratio ($\%$)", color='g')
#ax2.tick_params('y', colors='g')

#major_ticks = np.arange(0, 10, 5)
#minor_ticks = np.arange(0, 25, 1)
#axes.set_xticks(major_ticks)
#axes.set_xticks(minor_ticks, minor=True)
#axes.tick_params('y', colors='b')
plt.grid(True)

#set legend
plt   = plt1+plt2 + plt3 + plt4 + plt5 + plt6
labs = [l.get_label() for l in plt]
#plt.legend(loc=2, frameon=False)
axes.legend(plt, labs, loc=4, frameon=False)
#figure.set_size_inches(7,4)
figure.savefig("../figs/aq_oxi_species_avcomp.pdf")
