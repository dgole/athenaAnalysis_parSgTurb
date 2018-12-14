#!/bin/bash

#SBATCH -J tcAll
#SBATCH --time=1:00:00
#SBATCH -N 2
#SBATCH --ntasks 120 
#SBATCH -p development

# The following commands will be executed when this script is run.

export DIR=$SCRATCH/parSgTurb/data/newBCs_turbScaling

module load gcc
module load python3

python3 tabCombiner.py 100 $DIR/run10/bin/ $DIR/run10_decay/bin/ $DIR/run12/bin/ $DIR/run14/bin/ $DIR/run22/bin/ $DIR/run30/bin/



#
