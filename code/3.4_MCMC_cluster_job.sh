#!/bin/bash
# Shell script for fitting baseline model - mdrvx - Chathika Weerasuriya <c.weerasuriya@lshtm.ac.uk>
#
#$ -N CNparex2
#$ -cwd -V
#$ -q long.q
#$ -m beas
#$ -M c.weerasuriya@lshtm.ac.uk
#$ -l mem_free=1.2G,h_vmem=1.5G
#$ -t 1-1000
#$ -j y
#$ -o /home/lsh1604836/clusterOE/

cd "$HOME"/MDR-Aeras/mdrvx/output/mcmc_sampler/ || exit
mkdir -p co

module load "R/3.5.3"
export MCMC_LOG_LEVEL=INFO
R CMD BATCH "$HOME"/MDR-Aeras/mdrvx/code/3.2_MCMC_cluster_controller.R co/mcmc_"${SGE_TASK_ID}".out
