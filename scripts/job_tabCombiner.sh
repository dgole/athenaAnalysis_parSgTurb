#!/bin/bash

#SBATCH -J tc
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -p development

# The following commands will be executed when this script is run.

export DIR=$SCRATCH/parSgTurb/data/newBCs_turbScaling

module load gcc
module load python3

python3 tabCombiner.py $DIR/run10/bin/ $DIR/run12/bin/ $DIR/run14/bin/ $DIR/run22/bin/ $DIR/run30/bin/




#
