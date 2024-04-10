#!/usr/bin/env bash

##########################################################################################################
#### DE analysis by random selection in RECOUNT3 data
#### JToubia - January 2023
#########

n_genes="$1" 

echo "Starting to run random tests"

echo ""
outdir="/data_b/tools/SCA/output/Test_Random_n${n_genes}"

for i in {1..100}
do
	echo "Running Test $i"
	mkdir -p ${outdir}/Test${i}
	cd ${outdir}/Test${i}
	Rscript --vanilla /data_b/tools/SCA/bin/Aux/r3_random_edgeR_plus.R ${n_genes}
	find . -type d -empty -delete
done
echo ""

echo "Finished running random tests"
