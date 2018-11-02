#!/bin/bash

for var in "$@"
do
    cp baseJob.sh job"$var".sh
		echo "python3 readSplitTabsAndExport.py $var" >> job"$var".sh
		sbatch -A ucb17_summit1 job"$var".sh
		rm job"$var".sh
done

#
