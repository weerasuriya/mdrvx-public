#!/bin/bash
# Script to contatenate all eq files into single files
# Brace expansion: https://unix.stackexchange.com/questions/385357/cat-files-in-specific-order-based-on-number-in-filename

./meow.awk ../2_main.R

cp ../2.7_eq_transit.R 08_transit1.R
cp ../2.8_calc_local.R 10_calcloc.R
cp ../2.9_calc_fittargets.R 12_calcfit.R

rm 2_main_rebuild.R

# Various helper/init
rm 02_init.R
cat ../1.1_library.R \
    ../2.1_initialise.R \
    ../2.2_name.R \
    ../2.3_vectors_nh.R \
    ../2.4_vectors_mort.R \
    ../2.5_vectors_rx.R \
    ../2.5.2_vectors_vx.R > 02_init.R

# First step (fstep) epi equations
rm 04_fstep.R
cat ../2.6.0_eq_susceptible.R \
    ../2.6.1_eq_ntds.R \
    ../2.6.2_eq_ntdr.R \
    ../2.6.3_eq_ptds.R \
    ../2.6.4_eq_ptdr.R > 04_fstep.R

# First step vaccine equations
rm 06_fstepvax.R
cat ../2.6.0_eq_susceptible_vx.R \
    ../2.6.1_eq_ntds_vx.R \
    ../2.6.2_eq_ntdr_vx.R \
    ../2.6.3_eq_ptds_vx.R \
    ../2.6.4_eq_ptdr_vx.R > 06_fstepvax.R

# 08_ is transit eqs - self contained

# Final assembly
cat 0{1..9}_*.R {10..100}_*.R > 2_main_rebuild.R 2>/dev/null

exit 0