#!/bin/bash
# Calling scripts for vaccine vs baseline runs
#$ -cwd -V
#$ -q long.q
#$ -m beas
#$ -M c.weerasuriya@lshtm.ac.uk
#$ -N CN-OG-vcovR
#$ -l mem_free=2G,h_vmem=4G
#$ -t 1-1000
#$ -R y
#$ -j y
#$ -o /home/lsh1604836/clusterOE

mkdir -p "$HOME"/MDR-Aeras/mdrvx/output/vbp/co
cd "$HOME"/MDR-Aeras/mdrvx/output/vbp || exit

module load 'R/3.5.3'
export MCMC_LOG_LEVEL=INFO

R CMD BATCH "$HOME"/MDR-Aeras/mdrvx/code/"3.12_vbp_v2.R" co/vbp_"${SGE_TASK_ID}".out
