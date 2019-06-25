#!/bin/bash

python3 planPlotter_advancedStats.py ../../data/prodRuns/run100/ 253 1000 > fits_100.txt
python3 planPlotter_advancedStats.py ../../data/prodRuns/run101/ 367 1000 > fits_101.txt
python3 planPlotter_advancedStats.py ../../data/prodRuns/run103/ 291 1000 > fits_103.txt

python3 clumpTracker_plots.py ../../data/prodRuns/run100/ 220 275 1000 > fits_ct_100.txt
python3 clumpTracker_plots.py ../../data/prodRuns/run101/ 300 390 1000 > fits_ct_101.txt
python3 clumpTracker_plots.py ../../data/prodRuns/run103/ 220 290 1000 > fits_ct_103.txt








#
