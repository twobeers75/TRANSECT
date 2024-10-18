#!/usr/bin/env bash

##########################################################################################################
#### Post DE analysis folder organisation
#### JToubia - January 2023
#### Crude organisation of results! Checking if file exists and moving it appropriately if so
#########

### Start with the files produced from the stratifiaction process
mkdir 01-Stratification
if [ -f "GOI_exp_raw_OG.tsv" ]; then
	mv GOI_exp_raw_OG.tsv 01-Stratification
fi
if [ -f "GOI_exp_with_strat.tsv" ]; then
	mv GOI_exp_with_strat.tsv 01-Stratification
fi
count_html=$(ls -1 *.html 2>/dev/null | wc -l)
if [ $count_html != 0 ]; then
	mv *.html 01-Stratification
fi

### Next, the files produced from the enrichment analyses
mkdir 03-Enrichment
if [ -d "GSEA" ]; then
	mv GSEA 03-Enrichment
fi
if [ -d "WebGestalt" ]; then
	mv WebGestalt 03-Enrichment
fi
if [ -f "gene_normalised_expression_data_raw_gsea.cls" ]; then
	mv gene_normalised_expression_data_raw_gsea.cls 03-Enrichment
fi
if [ -f "gene_normalised_expression_data_raw_gsea.txt" ]; then
	mv gene_normalised_expression_data_raw_gsea.txt 03-Enrichment
fi
if [ -f "top_tags_ranked.rnk" ]; then
	mv top_tags_ranked.rnk 03-Enrichment
fi

### Finally, everything left should be from the DE analysis
mkdir 02-DE
if [ -d "glimma-plots" ]; then
	mv glimma-plots 02-DE
fi
count_png=$(ls -1 *.png 2>/dev/null | wc -l)
if [ $count_png != 0 ]; then
	mv *.png 02-DE
fi
count_csv=$(ls -1 *.csv 2>/dev/null | wc -l)
if [ $count_csv != 0 ]; then
	mv *.csv 02-DE
fi
count_tsv=$(ls -1 *.tsv 2>/dev/null | wc -l)
if [ $count_tsv != 0 ]; then
	mv *.tsv 02-DE
fi



