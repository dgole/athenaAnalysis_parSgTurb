#!/bin/bash

#SBATCH -J 3d_100
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH --tasks-per-node 1
#SBATCH -p development

# The following commands will be executed when this script is run.
module load gcc
module load python3


DIR=$SCRATCH/parSgTurb/data/prodRuns/run100/bin/
python3 pspec_justMP_2D.py $DIR 19 21 255 257






#
