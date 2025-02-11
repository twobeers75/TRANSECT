#!/usr/bin/env bash

##########################################################################################################
#### DE analysis by random selection in RECOUNT3 data
#### JToubia - January 2023
#########

#### MODIFY THE FOLLOWING PATH TO SUIT YOUR SETUP
TRANSECT_HOME="/path/to/TRANSECT"

#### Should not require any modifications below this line
n_genes="$1"

echo "Starting to run random tests"

echo ""
outdir="${TRANSECT_HOME}/output/Test_Random_n${n_genes}"

for i in {1..100}
do
	echo "Running Test $i"
	mkdir -p ${outdir}/Test${i}
	cd ${outdir}/Test${i}
	Rscript --vanilla ${TRANSECT_HOME}/bin/Aux/r3_random_edgeR_plus.R ${n_genes}
	find . -type d -empty -delete
done
echo ""

echo "Finished running random tests"
