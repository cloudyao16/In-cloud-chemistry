#!/usr/bin/env python

import sys, os
sys.path.append("../../../aq_partmc/tool/")
import mpl_helper
import matplotlib
import partmc
import scipy.io, numpy

for (filename, index) in partmc.get_filename_list('../mono/out/', r'urban_plume_aq_chem_0000000([1-9]+)_process\.nc'):
    (figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)
   
    print(filename)
    ncf = scipy.io.netcdf_file(filename)
    ind = int(index) - 1
    diam_edges = ncf.variables["diam_edges"].data * 1e6
    ph_edges = ncf.variables["ph_edges"].data 
    diam_ph_dist = ncf.variables["diam_ph_dist_neut"].data * 1e-6

    p = axes.pcolor(diam_edges, ph_edges, diam_ph_dist,
                    norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e5), linewidths = 0.1)

    axes.annotate(r"$\rm %s\; min$"%ind, xy=(0.015, 8), xycoords="data",
                  va="center", ha="center",size=12)
    axes.set_xscale("log")
    axes.set_xlabel(r"Dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-2, 1e1)

    axes.set_yscale("linear")
    axes.set_ylabel(r"Sulfate mass fraction $w_{\rm sulfate}$ / \%")
    axes.set_ylim(0,10)

    axes.grid(True)
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                           orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number conc. $n(D_{\rm dry},w_{\rm SO_4^{2-}})$ / $\rm cm^{-3}$")

    #out_filename = "figure_16Nov/aq_parcel_diam_so4_dist_initial.pdf"
    out_filename = "../figs/ph_2d_%s.pdf"%ind

    figure.savefig(out_filename)
    print(out_filename)
