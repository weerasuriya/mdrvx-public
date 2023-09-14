#!/bin/bash
# Shell script for fitting baseline model - mdrvx - Chathika Weerasuriya <c.weerasuriya@lshtm.ac.uk>
#
#$ -N CN_rs_prepub
#$ -cwd -V
#$ -q long.q
#$ -m beas
#$ -M c.weerasuriya@lshtm.ac.uk
#$ -l mem_free=1G,h_vmem=1.5G
#$ -t 1-10000
#$ -j y
#$ -o /home/lsh1604836/clusterOE/

timestamp=$(date +%Y%m%d-%H%M -r "$HOME"/MDR-Aeras/mdrvx/code/3.4_MCMC_cluster_job.sh)
cd "$HOME"/MDR-Aeras/mdrvx/output/random_sampler/ || exit

export OPPATH=$(pwd)
mkdir -p co
R CMD BATCH $HOME/MDR-Aeras/mdrvx/code/3_RS_cluster_controller.R co/rs_${SGE_TASK_ID}.out
