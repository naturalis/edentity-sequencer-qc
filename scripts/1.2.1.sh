#!/bin/bash

files=$(ls ./)
for file in $files
do
	nohup python ~/edentity-sequencer-qc/scripts/1.2.1.py $file >> ~/edentity-sequencer-qc/results/elements-1.2.1-homopolymer_decay.csv  &
done