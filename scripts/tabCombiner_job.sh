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

while true
do
	python3 tabCombiner_big.py 64 $DIR/run200/bin/
	python3 tabCombiner_big.py 64 $DIR/run201/bin/
done

#
