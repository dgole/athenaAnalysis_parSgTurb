#!/bin/bash

#SBATCH -J tcAll
#SBATCH --time=6:00:00
#SBATCH -N 1
#SBATCH --ntasks 64
#SBATCH -p normal 

# The following commands will be executed when this script is run.

export DIR=$SCRATCH/parSgTurb/data/prodRuns

module load gcc
module load python3

while true
do
	python3 tabCombiner_big.py 64 $DIR/run100/bin/

	python3 tabCombiner_big.py 64 $DIR/run60/bin/ 
	python3 tabCombiner_big.py 64 $DIR/run61/bin/ 
	python3 tabCombiner_big.py 64 $DIR/run62/bin/
	python3 tabCombiner_big.py 64 $DIR/run63/bin/  	

	python3 tabCombiner_big.py 64 $DIR/run70/bin/
	python3 tabCombiner_big.py 64 $DIR/run71/bin/
	python3 tabCombiner_big.py 64 $DIR/run72/bin/

	python3 tabCombiner_big.py 64 $DIR/run80/bin/
	python3 tabCombiner_big.py 64 $DIR/run81/bin/
	python3 tabCombiner_big.py 64 $DIR/run82/bin/

        python3 tabCombiner_big.py 64 $DIR/run90/bin/
        python3 tabCombiner_big.py 64 $DIR/run91/bin/
        python3 tabCombiner_big.py 64 $DIR/run92/bin/

        python3 tabCombiner_big.py 64 $DIR/run110/bin/
        python3 tabCombiner_big.py 64 $DIR/run111/bin/
        python3 tabCombiner_big.py 64 $DIR/run112/bin/

done




#
