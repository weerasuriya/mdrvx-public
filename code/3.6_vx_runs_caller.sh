#!/bin/bash
# Calling scripts for vaccine vs baseline runs
#$ -cwd -V
#$ -q short.q
#$ -m beas
#$ -M c.weerasuriya@lshtm.ac.uk
#$ -N CN7wanefixvxe
#$ -l mem_free=2G,h_vmem=4G
#$ -t 1-1000
#$ -R y
#$ -j y
#$ -o /home/lsh1604836/clusterOE

mkdir -p "$HOME"/MDR-Aeras/mdrvx/output/vxe/co
cd "$HOME"/MDR-Aeras/mdrvx/output/vxe || exit

module load 'R/3.5.3'
export MCMC_LOG_LEVEL=INFO

R CMD BATCH "$HOME"/MDR-Aeras/mdrvx/code/"3.10_vxe_generator.R" co/vxe_"${SGE_TASK_ID}".out
