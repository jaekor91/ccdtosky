#!/bin/bash
start_idx=0
last_idx=3400000
step=1700000
for (( i=$start_idx, j=$start_idx+$step; i <= $last_idx; i+=$step, j+=$step ))
do
	echo "Calculating chunk [$i, $j]"
	python runscript_run_per_chunk.py $i $j >> Nside11.out
done