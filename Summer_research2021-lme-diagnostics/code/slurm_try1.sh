#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --partition=StudentSlurm
#SBATCH --nodelist=bio338
#SBATCH --cpus-per-task=12
#SBATCH --job-name=extra_tryout
#SBATCH --error=error_resid_extra.out
#SBATCH --output=resid_extra_tryout_op.out
/opt/R/4.0.4/lib/R/bin/Rscript --vanilla main_resid_extra.R