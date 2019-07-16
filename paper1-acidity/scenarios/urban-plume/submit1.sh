#!/bin/tcsh

#SBATCH --job-name=yuyao
#SBATCH --time=12:00:00
#SBATCH -c 12
#SBATCH -N 1
#SBATCH --mem=15000
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=yuyao3@illinois.edu
#SBATCH -o urban-ph.o%j # Output and error file name (%j expands to jobID)
echo -n "Job starting at "; date
# Change to the directory to run the executable
# Run the executable
/data/keeling/a/yuyao3/c/aq_partmc/build/partmc urban_plume_aq_chem_5h_avcomp_onebin.spec 

echo ' '
echo -n "Job finished at "; date
