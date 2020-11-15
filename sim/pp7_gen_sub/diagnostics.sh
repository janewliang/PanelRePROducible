#!/bin/bash

#SBATCH -J diagnostics #A single job name for the array
#SBATCH -o results/cluster/diagnostics__%A_%a.out #Standard output
#SBATCH -e results/cluster/diagnostics__%A_%a.err #Standard error
#SBATCH -p serial_requeue #Partition
#SBATCH -t 1440         #Runtime in minutes
#SBATCH --mem-per-cpu=40000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH --mail-type=END # Email
#SBATCH --mail-user=jwliang@g.harvard.edu

R CMD BATCH --no-restore --no-save ../../_deps/pp7/diagnostics.R results/cluster/diagnostics.Rout