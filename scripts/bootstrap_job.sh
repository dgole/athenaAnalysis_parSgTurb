#!/bin/bash

#SBATCH -J bs_bcpl_103
#SBATCH --time=12:00:00
#SBATCH -N 16
#SBATCH --ntasks 1024
#SBATCH -p normal

# The following commands will be executed when this script is run.
module load gcc
module load python3

#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_100/ 220 275 1000 spl
#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_100/ 220 275 1000 stpl
#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_100/ 220 275 1000 bpl
#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_100/ 220 275 1000 bcpl

#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_101/ 300 390 1000 spl
#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_101/ 300 390 1000 stpl
#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_101/ 300 390 1000 bpl
#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_101/ 300 390 1000 bcpl

#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_103/ 220 290 1000 spl
#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_103/ 220 290 1000 stpl
#python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_103/ 220 290 1000 bpl
python3 ct_fits_mcmc_bootstrap_parallel.py ./planNew2/run_plan_103/ 220 290 1000 bcpl









#
