#!/bin/bash

#SBATCH -J 3d
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 64
#SBATCH -p development

# The following commands will be executed when this script is run.
module load gcc
module load python3

DIR=$SCRATCH/parSgTurb/data/prodRuns/run100/bin/

python3 slices.py $DIR 2







#
