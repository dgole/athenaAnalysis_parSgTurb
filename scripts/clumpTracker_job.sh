#!/bin/bash

#SBATCH -J ct_100
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -p development

# The following commands will be executed when this script is run.
module load gcc
module load python3

DIR=$SCRATCH/parSgTurb/data/prodRuns/run100/bin/

python3 clumpTracker.py $DIR 220 289







#
