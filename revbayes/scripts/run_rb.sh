#!/bin/bash

#SBATCH --nodes=1 # 1 node per job
#SBATCH --time=15-00:00:00 # no idea how long this will take, let's be safe 
#SBATCH --mem=16G

#SBATCH --output=output/jobs/job_%j.out
#SBATCH --error=output/jobs/job_%j.err

#SBATCH --job-name="fbdhisse_canidae"

#SBATCH --mail-user=petrucci@iastate.edu   # my e-mail
#SBATCH --mail-type=BEGIN # get notifications for all job cases
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# load modules
module purge
module load revbayes/dev-tp-sa-fix.tensorphylo

# go to the correct directory
cd /work/LAS/phylo-lab/petrucci/canids_chap2/revbayes

# source it, the parameter combination, and the actual script
rb scripts/fbdhisse_lik.Rev 
