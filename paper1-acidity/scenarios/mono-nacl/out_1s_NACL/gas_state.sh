#!/bin/bash

module load gnu/cdo-1.6.3-gnu-4.4.6

for f in urban_plume_aq_chem_0001*.nc
do
	echo "Processing $f"
	cdo selvar,gas_mixing_ratio $f obsolete_$f 
done
cdo merge obsolete* gas_state.nc
rm obsolete*
