#!/bin/bash

#SBATCH -J ct_103_ht
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH -p normal

# The following commands will be executed when this script is run.
module load gcc
module load python3

#DIR=./planNew2/run_plan_100/
#python3 clumpTracker.py $DIR 220 275
#DIR=./plan_highDensThresh/run_plan_100/
#python3 clumpTracker.py $DIR 220 275
#DIR=./plan_lowDensThresh/run_plan_100/
#python3 clumpTracker.py $DIR 220 240

#DIR=./planNew2/run_plan_101/
#python3 clumpTracker.py $DIR 300 390
#DIR=./plan_highDensThresh/run_plan_101/
#python3 clumpTracker.py $DIR 300 390

#DIR=./planNew2/run_plan_103/
#python3 clumpTracker.py $DIR 220 290
DIR=./plan_highDensThresh/run_plan_103/
python3 clumpTracker.py $DIR 200 290







#
