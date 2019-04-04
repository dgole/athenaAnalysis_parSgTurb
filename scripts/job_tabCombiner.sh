#!/bin/bash

#SBATCH -J tcAll
#SBATCH --time=2:00:00
#SBATCH -N 1
#SBATCH --ntasks 64
#SBATCH -p development

# The following commands will be executed when this script is run.

export DIR=$SCRATCH/parSgTurb/data/kspaceTest

module load gcc
module load python3

while true
do
	python3 tabCombiner_big.py 64 $DIR/run120/bin/
        python3 tabCombiner_big.py 64 $DIR/run121/bin/
        python3 tabCombiner_big.py 64 $DIR/run122/bin/
        python3 tabCombiner_big.py 64 $DIR/run123/bin/

        python3 tabCombiner_big.py 64 $DIR/run150/bin/
        python3 tabCombiner_big.py 64 $DIR/run151/bin/
        python3 tabCombiner_big.py 64 $DIR/run152/bin/
        python3 tabCombiner_big.py 64 $DIR/run153/bin/
        python3 tabCombiner_big.py 64 $DIR/run154/bin/

        python3 tabCombiner_big.py 64 $DIR/run160/bin/
        python3 tabCombiner_big.py 64 $DIR/run161/bin/
        python3 tabCombiner_big.py 64 $DIR/run162/bin/
        python3 tabCombiner_big.py 64 $DIR/run163/bin/
        python3 tabCombiner_big.py 64 $DIR/run164/bin/

        python3 tabCombiner_big.py 64 $DIR/run170/bin/
        python3 tabCombiner_big.py 64 $DIR/run171/bin/
        python3 tabCombiner_big.py 64 $DIR/run172/bin/
        python3 tabCombiner_big.py 64 $DIR/run173/bin/
        python3 tabCombiner_big.py 64 $DIR/run174/bin/

        python3 tabCombiner_big.py 64 $DIR/run180/bin/
        python3 tabCombiner_big.py 64 $DIR/run181/bin/
        python3 tabCombiner_big.py 64 $DIR/run182/bin/
        python3 tabCombiner_big.py 64 $DIR/run183/bin/
        python3 tabCombiner_big.py 64 $DIR/run184/bin/



done

#
