#!/usr/bin/env python

import sys, os
sys.path.append("../../../aq_partmc/tool/")
import mpl_helper
import scipy.io, numpy
import matplotlib as mpl
#mpl.rcParams['font.size'] = 15

(figure, axes) = mpl_helper.make_fig(right_margin=0.8)
axes2 = axes.twinx()

ncf = scipy.io.netcdf_file("../urban-plume/out/urban_plume_process.nc")
time = ncf.variables["time"].data / 3600
num_conc = ncf.variables["tot_num_conc"].data / 1e6
#num_conc_err = ncf.variables["tot_num_conc_ci_offset"].data / 1e6
mass_conc = ncf.variables["tot_mass_conc"].data * 1e9
#mass_conc_err = ncf.variables["tot_mass_conc_ci_offset"].data * 1e9

plt1=axes.plot(time[0:24], num_conc[0:24], color="blue", label="num. conc.", linewidth = 1.5)
plt2=axes2.plot(time[0:24], mass_conc[0:24], color="green", label = "mass conc.", linewidth =1.5)
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"num. conc. / $\rm cm^{-3}$")
axes2.set_ylabel(r"mass conc. / $\rm \mu g\ m^{-3}$")
axes.grid(True)
plt4  = plt1+plt2
labs = [l.get_label() for l in plt4]
axes.legend(plt4, labs, loc=2, frameon=False)

figure.savefig("../figs/urban_plume_total.pdf")
