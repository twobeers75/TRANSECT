#!/usr/bin/env bash

##########################################################################################################
#### Post DE analysis folder organisation
#### JToubia - January 2023
#########

mkdir 01-Stratification
mv *.html GOI_exp_raw_OG.tsv GOI_exp_with_strat.tsv 01-Stratification

mkdir 03-Enrichment
if [ -d "GSEA" ]; then
	mv GSEA 03-Enrichment
fi
mv WebGestalt gene_normalised_expression_data_raw_gsea.cls gene_normalised_expression_data_raw_gsea.txt top_tags_ranked.rnk 03-Enrichment

mkdir 02-DE
mv glimma-plots *.png *.csv *.tsv 02-DE



