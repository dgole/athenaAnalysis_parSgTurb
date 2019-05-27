#!/bin/bash

#SBATCH -J ct_100
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -p normal

# The following commands will be executed when this script is run.
module load gcc
module load python3

DIR=./planNew2/run_plan_100/

python3 clumpTracker.py $DIR 220 289







#
