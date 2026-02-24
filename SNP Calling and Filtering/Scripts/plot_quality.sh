#!/bin/bash
#SBATCH --job-name=plot_quality      
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.sow@ufl.edu    
#SBATCH --ntasks=1
#SBATCH --mem=1gb   
#SBATCH --time=00:20:00   # Walltime
#SBATCH --output=plot_quality.%j.out   
date; hostname; pwd

module load R

Rscript plot_quality.R

see here for notes on preparing R https://docs.rc.ufl.edu/software/apps/r/#job-script-examples
