#!/usr/bin/env Rscript

###*****************************************************************************
### Retrieve RECOUNT3 data
### JToubia - January 2023
###*****************************************************************************
message("Downloading data from RECOUNT3 repository")

###*****************************************************************************
### Import libraries ####
###*****************************************************************************
message("Loading required packages")
suppressMessages(if (!require("pacman")) install.packages("pacman"))
p_load(recount, recount3, DEFormats, data.table)

###*****************************************************************************
### Read in Args ####
###*****************************************************************************
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (ProjectID)", call.=FALSE)
}

pID <- args[1]
date_stamp <- args[2]

###*****************************************************************************
### Setup for retrieval of data ####
###*****************************************************************************
main_dir <- getwd()
file_out_suffix <- paste("R3-", pID, "-", date_stamp, sep="")

###*****************************************************************************
### Start ####
###*****************************************************************************
message(paste("RECOUNT3 data retrieval for", pID))

### get the available projects in recount3
human_projects <- available_projects()

# ### infor about TCGA and GTEx samples
# subset(human_projects, file_source == "tcga" & project_type == "data_sources")
# subset(human_projects, file_source == "gtex" & project_type == "data_sources")

### get information required for downloading a project using SRA ID or TCGA/GTEx code
proj_info <- subset(
  human_projects,
  project == pID & project_type == "data_sources"
)

### pull the data based on the info gathered above
message("Starting data retrieval")
rse_gene_data <- create_rse(proj_info)
message(paste("Raw RECOUNT3 data for", pID, "retrieved"))

### Scale the raw counts and create TPMs as well
message("Scale raw counts")
assay(rse_gene_data, "counts") <- transform_counts(rse_gene_data)
message("Generate TPMs")
assay(rse_gene_data, "TPM") <- recount::getTPM(rse_gene_data)
message("Following assays now retrieved")
message(paste(assayNames(rse_gene_data), collapse = ' & '))

message("Creating formatted copies of data")
### pull out count matrix from rse object
gene_counts_df <- as.data.frame(assay(rse_gene_data, "counts"))
gene_tpm_df <- as.data.frame(assay(rse_gene_data, "TPM"))
#gene_counts_df <- assays(rse_gene_data)$counts

### create a gene name lookup table
#rowRanges(rse_gene_data)[,c("gene_id","gene_name")]
gene_name_lookup <- mcols(rse_gene_data, use.names=FALSE)[,c("gene_id","gene_name")]
rownames(gene_name_lookup) <- gene_name_lookup$gene_id
## lets see how many duplicates there are in gene_name
# n_occur <- data.frame(table(gene_name_lookup$gene_name))
# n_occur[n_occur$Freq > 1,]
# dim(n_occur[n_occur$Freq > 1,])

### create a sample ID lookup table
if(proj_info$file_source == "tcga"){
  sample_name_lookup <- colData(rse_gene_data)[, c("rail_id","tcga.tcga_barcode")]
  sample_name_lookup$tcga.tcga_barcode <- substr(sample_name_lookup$tcga.tcga_barcode,1,nchar(sample_name_lookup$tcga.tcga_barcode)-12)
  colnames(sample_name_lookup) <- c("rail_id","sampleID")
} else if(proj_info$file_source == "gtex"){
  sample_name_lookup <- colData(rse_gene_data)[, c("rail_id","gtex.sampid")]
  colnames(sample_name_lookup) <- c("rail_id","sampleID")
} else {
  stop("Execution stopped. Project source is not TCGA nor GTEx")
}
## lets see how many duplicates there are in gene_name
# n_occur <- data.frame(table(sample_name_lookup$tcga.tcga_barcode))
# n_occur[n_occur$Freq > 1,]
# dim(n_occur[n_occur$Freq > 1,])
# DAMN! There is!! Have to deal with this too
# all(gene_counts_df[, c("bc7b5d4c-4668-4fba-99ab-88ce5ef041b7")] == gene_counts_df[, c("6015ac2c-3f1c-4072-a8c4-bc40110152f3")])
# [1] TRUE
# all(gene_counts_df[, c("bc7b5d4c-4668-4fba-99ab-88ce5ef041b7")] == gene_counts_df[, c("ea81d034-0f2d-46e6-90fa-769c62b91ef8")])
# [1] FALSE

### resetting row and col names for counts
## gene names have duplicates so we need to deal with these. We will keep the most expressed across patients
message("Working on count assay")
gene_counts_df_mod <- cbind(gene_name = gene_name_lookup$gene_name, row_sum = rowSums(gene_counts_df), gene_counts_df)
gene_counts_df_mod_sorted <- gene_counts_df_mod[order(gene_counts_df_mod$gene_name, -abs(gene_counts_df_mod$row_sum) ), ]
gene_counts_df_mod_sorted_dedup <- gene_counts_df_mod_sorted[ !duplicated(gene_counts_df_mod_sorted$gene_name), ]
rownames(gene_counts_df_mod_sorted_dedup) <- gene_counts_df_mod_sorted_dedup$gene_name
gene_counts_df_mod_sorted_dedup <- subset(gene_counts_df_mod_sorted_dedup, select = -c(row_sum))

colnames(gene_counts_df_mod_sorted_dedup) <- c("gene_name", sample_name_lookup$sampleID)
gene_counts_df_mod_sorted_dedup <- gene_counts_df_mod_sorted_dedup[, !duplicated(colnames(gene_counts_df_mod_sorted_dedup))]
message("Dimensions of raw table")
message(paste(dim(gene_counts_df), collapse = ' & '))
message("Dimensions of processed table")
message(paste(dim(gene_counts_df_mod_sorted_dedup), collapse = ' & '))

## same again for TPM table
## sample names have duplicates too. Looking into this I see same IDs with completely duplicated columns and some with differences?? Will just take the first!
message("Working on TPM assay")
gene_tpm_df_mod <- cbind(gene_name = gene_name_lookup$gene_name, row_sum = rowSums(gene_counts_df), gene_tpm_df)
gene_tpm_df_mod_sorted <- gene_tpm_df_mod[order(gene_tpm_df_mod$gene_name, -abs(gene_tpm_df_mod$row_sum) ), ]
gene_tpm_df_mod_sorted_dedup <- gene_tpm_df_mod_sorted[ !duplicated(gene_tpm_df_mod_sorted$gene_name), ]
rownames(gene_tpm_df_mod_sorted_dedup) <- gene_tpm_df_mod_sorted_dedup$gene_name
gene_tpm_df_mod_sorted_dedup <- subset(gene_tpm_df_mod_sorted_dedup, select = -c(row_sum))
gene_tpm_df_mod_sorted_dedup[,-1] <- round(gene_tpm_df_mod_sorted_dedup[,-1], 2)

colnames(gene_tpm_df_mod_sorted_dedup) <- c("gene_name", sample_name_lookup$sampleID)
gene_tpm_df_mod_sorted_dedup <- gene_tpm_df_mod_sorted_dedup[, !duplicated(colnames(gene_tpm_df_mod_sorted_dedup))]
message("Dimensions of raw table")
message(paste(dim(gene_tpm_df), collapse = ' & '))
message("Dimensions of processed table")
message(paste(dim(gene_tpm_df_mod_sorted_dedup), collapse = ' & '))

# ### removing 7SK from tables
# row_names_to_remove <- c("7SK")
# gene_counts_df_mod_sorted_dedup <- gene_counts_df_mod_sorted_dedup[!(row.names(gene_counts_df_mod_sorted_dedup) %in% row_names_to_remove),]
# gene_tpm_df_mod_sorted_dedup <- gene_tpm_df_mod_sorted_dedup[!(row.names(gene_tpm_df_mod_sorted_dedup) %in% row_names_to_remove),]
message("Finished creating formatted copies of data")

### write tables to file
message("Writing data to file")
write.table(gene_counts_df_mod_sorted_dedup, paste(file_out_suffix, "_count-mRNA.tsv", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_tpm_df_mod_sorted_dedup, paste(file_out_suffix, "_tpm-mRNA.tsv", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")

### Misc. summary stats
# ### some summary stuff
# temp <- t(gene_counts_df_mod_sorted_dedup[, -1])
# temp_summ <- summary(temp)
# write.table(t(temp_summ), "temp_sum.tsv")

message(paste("Finished RECOUNT3 data retrieval for", pID))

