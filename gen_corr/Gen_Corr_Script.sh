#!/bin/sh
#$ -cwd
#$ -l h_rt=24:00:00
#$ -V
#$ -l h_vmem=16G

. /etc/profile.d/modules.sh
module load R

R CMD BATCH Gen_Corr_$SGE_TASK_ID.R