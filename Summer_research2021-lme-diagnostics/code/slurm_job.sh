#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --partition=StudentSlurm
#SBATCH --nodelist=crug2
#SBATCH --cpus-per-task=12
#SBATCH --job-name=lme_resid
#SBATCH --error=error_resid_extra.out
#SBATCH --output=slurm_resid_extra.out
/opt/R/4.0.4/lib/R/bin/Rscript --vanilla main_resid_extra.R