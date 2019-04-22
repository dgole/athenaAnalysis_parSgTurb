#!/bin/bash

#SBATCH -J tcAll
#SBATCH --time=2:00:00
#SBATCH -N 16
#SBATCH --ntasks 1024
#SBATCH -p normal

# The following commands will be executed when this script is run.

DIR=$SCRATCH/parSgTurb/data/prodRuns

module load gcc
module load python3

python3 tabCombiner_big.py 1024 $DIR/run100/bin/
python3 tabCombiner_big.py 1024 $DIR/run101/bin/


#
