#!/bin/bash

#SBATCH --nodes=1 # 1 node per job
#SBATCH --time=01:00:00 # just an hour--calling the array itself doesn't take long
#SBATCH --array=1-18 # 18 parameter combinations

#SBATCH --output=output/bisse/jobs/array/job_%A_%a.out
#SBATCH --error=output/bisse/jobs/array/job_%A_%a.err

#SBATCH --job-name="array_fbdsse_canids"

#SBATCH --mail-user=petrucci@iastate.edu   # my e-mail
#SBATCH --mail-type=BEGIN # get notifications for all job cases
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# go to the correct directory
cd /work/LAS/phylo-lab/petrucci/canids_chap2/revbayes

# source it, the parameter combination, and the actual script
rb ref_${SLURM_ARRAY_TASK_ID}.Rev master.Rev
