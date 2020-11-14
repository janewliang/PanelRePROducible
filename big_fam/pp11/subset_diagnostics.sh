#!/bin/bash

#SBATCH -J subset_diagnostics #A single job name for the array
#SBATCH -o results/cluster/subset_diagnostics__%A_%a.out #Standard output
#SBATCH -e results/cluster/subset_diagnostics__%A_%a.err #Standard error
#SBATCH -p serial_requeue #Partition
#SBATCH -t 1440         #Runtime in minutes
#SBATCH --mem-per-cpu=10000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH --mail-type=END # Email
#SBATCH --mail-user=jwliang@g.harvard.edu

R CMD BATCH --no-restore --no-save ../../_scripts/pp11/subset_diagnostics.R results/cluster/subset_diagnostics.Rout
