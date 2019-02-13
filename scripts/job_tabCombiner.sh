#!/bin/bash

#SBATCH -J tcAll
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH --ntasks 64 
#SBATCH -p development

# The following commands will be executed when this script is run.

export DIR=$SCRATCH/parSgTurb/data/prodRuns

module load gcc
module load python3

python3 tabCombiner.py 50 $DIR/run10/bin/ $DIR/run11/bin/ $DIR/run12/bin/ $DIR/run13/bin/ 



#
