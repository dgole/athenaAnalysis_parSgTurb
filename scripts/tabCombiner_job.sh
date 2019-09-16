#!/bin/bash

#SBATCH -J tcAll
#SBATCH --time=2:00:00
#SBATCH -N 8
#SBATCH --ntasks 512
#SBATCH -p development 

# The following commands will be executed when this script is run.

DIR=$SCRATCH/parSgTurb/data/prodRuns

module load gcc
module load python3

n=512
python3 tabCombiner_big.py $n $DIR/run602/bin/
#python3 tabCombiner_big.py $n $DIR/run102_noPar/bin/

n=128
python3 tabCombiner_big.py $n $DIR/run603/bin/

#python3 tabCombiner_big.py $n $DIR/run103/bin/
#python3 tabCombiner_big.py $n $DIR/run102/bin/
#python3 tabCombiner_big.py $n $DIR/run101/bin/
#python3 tabCombiner_big.py $n $DIR/run100/bin/






#
