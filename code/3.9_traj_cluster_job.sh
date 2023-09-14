#!/bin/bash
#$ -cwd -V
#$ -q tiny.q
#$ -m beas
#$ -M c.weerasuriya@lshtm.ac.uk
#$ -N CN-SW-traj
#$ -l mem_free=2G,h_vmem=2G
#$ -R y
#$ -j y
#$ -o /home/lsh1604836/clusterOE
#$ -t 1-1000
# TRAJ-POL - in loop - with GC and multiple scenarios.

cd "$HOME"/MDR-Aeras/mdrvx/output/trajectories/CN || exit

module load 'R/3.5.3'
export MCMC_LOG_LEVEL=INFO

R CMD BATCH "$HOME"/MDR-Aeras/mdrvx/code/3.8_trajectories.R co/"${SGE_TASK_ID}"_co_traj.out
