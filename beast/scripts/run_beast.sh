#!/bin/bash

#SBATCH --nodes=1 # 1 node per job
#SBATCH --time=31-00:00:00 # no idea how long this will take, let's be safe 

#SBATCH --output=output/jobs/job_%j.out
#SBATCH --error=output/jobs/job_%j.err

#SBATCH --job-name="4_canidae_run"

#SBATCH --mail-user=petrucci@iastate.edu   # my e-mail
#SBATCH --mail-type=BEGIN # get notifications for all job cases
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

cd /work/LAS/phylo-lab/petrucci/canids_chap2/beast

$HOME/beast/bin/beast scripts/canidae_run.xml
