#!/usr/bin/env Rscript

###*****************************************************************************
### DO DE analysis on RECODE3 random data
### JToubia - January 2023
###*****************************************************************************

###*****************************************************************************
### Import libraries ####
###*****************************************************************************
suppressMessages(if (!require("pacman")) install.packages("pacman"))
p_load(edgeR, Glimma, dplyr, tidyr, data.table, tibble, 
       ggplot2, ggforce, calibrate, gplots, RColorBrewer,
       rlogging)

SetLogFile()

###*****************************************************************************
#### Hard coded system variables ####
###*****************************************************************************


###*****************************************************************************
#### Read in Args ####
###*****************************************************************************
args <- commandArgs(trailingOnly = TRUE)

n_genes <- args[1]
message(paste("Random trials using ", n_genes, " per trial", sep = ""))

GOI <- "RANDOM"
rna_species <- "mRNA"
percentile <- as.integer(5)
switch <- "false"
ref_files_folder <- "/data_b/tools/SCA/REF_FILES"
# gsea_exe <- paste("/data_b/tools/SCA/bin/GSEA_Linux_4.2.3/gsea-cli.sh", "GSEAPreranked", sep =" ")
gsea_exe <- "false"
gene_fpkm_filename <- "/data_b/tools/SCA/data/RECOUNT3/BRCA/R3-BRCA-2024-03-07_tpm-mRNA.tsv"
gene_counts_filename <- "/data_b/tools/SCA/data/RECOUNT3/BRCA/R3-BRCA-2024-03-07_count-mRNA.tsv"
isomir_rpm_filename <- ""
normals <- ""

###*****************************************************************************
### Setup for analyses ####
###*****************************************************************************
main_dir <- getwd() 
dir.create(file.path(main_dir, "DE_Analysis"))
setwd(file.path(main_dir, "DE_Analysis"))

###*****************************************************************************
### Setup for random selection ####
###*****************************************************************************
message("Beginning random section of individuals from cohort")

### Check what species of RNA we are analysing and set/read appropriate variables
if (rna_species == "mRNA"){
  exp_data <- fread(file=gene_fpkm_filename, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
  exp_data <- exp_data %>% column_to_rownames('gene_name')
  exp_type <- "TPM"
} else if (rna_species == "miRNA"){
  mRNA_ids <- t(fread(file=gene_fpkm_filename, header = FALSE, nrows = 1))
  mRNA_ids <- mRNA_ids[-1]
  exp_data <- fread(file=isomir_rpm_filename, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
  exp_data <- exp_data %>% column_to_rownames('miR_name')
  exp_type <- "RPM"
} else {
  stop("Can only stratify by mRNA or miRNA", call=FALSE)
}
message(paste("Normalised expression data loaded for", GOI, ":", ncol(exp_data), "patients and", nrow(exp_data), "genes"))

TCGA_IDs <- colnames(exp_data)
random_low <- sample(TCGA_IDs, n_genes)
TCGA_IDs <- TCGA_IDs[!TCGA_IDs %in% random_low] 
random_high <- sample(TCGA_IDs, n_genes)

### now create design matrix
low_percentile_df <- as.data.frame(random_low)
matrix_colnames <- c("Sample")
colnames(low_percentile_df) <- matrix_colnames
low_percentile_df$lo <- 1
low_percentile_df$hi <- 0
message(paste(nrow(low_percentile_df), "patients stratified into low group for", GOI))

high_percentile_df <- as.data.frame(random_high)
colnames(high_percentile_df) <- matrix_colnames
high_percentile_df$lo <- 0
high_percentile_df$hi <- 1
message(paste(nrow(high_percentile_df), "patients stratified into high group for", GOI))

design <- rbind(low_percentile_df, high_percentile_df)
rownames(design) <- design$Sample
design <- design[c("lo", "hi")]

### check and strip normals out if requested
if (normals == "") {
  message(paste("Excluding non-diseased samples"))
  design$cancertype <- substr(rownames(design), 14, 15)
  design <- design[design$cancertype != 11, ]
  design <- subset(design, select=-cancertype)
}

write.table(design, "design.tsv", col.names=NA)

### all done here, sign off stratification
message("Finished random section of individuals from cohort")

###*****************************************************************************
### Setup for DE ####
###*****************************************************************************

message(paste("Beginning DE Analysis for", GOI))

hi <- 'hi'
lo <- 'lo'

### Read in data
raw_data <- fread(file=gene_counts_filename, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
colnames(raw_data)[1] <- 'gene_name'
message(paste("Raw count expression data loaded for", GOI, ":", ncol(raw_data)-1, "patients and", nrow(raw_data), "genes"))

GOI_label <- GOI

### Subset data
# dim(raw_data)
gene_counts <- raw_data[row.names(design)]
gene_counts[is.na(gene_counts)] <- 0
gene_names <- raw_data['gene_name']

### Setup groups
group <- vector()
for (row in 1:nrow(design)) {
	row_condition <- NaN
	for (condition in names(design)) {
		if (design[row, condition] != 0) {
			if (!is.na(row_condition)) {
				stop(paste("row ", row, " contains more than 1 condition"));
			} else {
				row_condition <- condition;
			}
		}
	}
	group <- append(group, row_condition);
}
group <- factor(group)
group <- relevel(group, "lo")

### Create DGEList object
y <- DGEList(counts=gene_counts, genes=gene_names, group=group)
# dim(y)

### CPM table
cpm_table <- round(cpm(y$counts),2)
write.csv(cbind(y$genes[rownames(cpm_table),1],cpm_table), "gene_raw_expression_data_cpm.csv", row.names=FALSE)

### Filtering and normalisation
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$gene_name)
y <- y[!d,]
# dim(y)

### filter lowly expressed genes
smallest_replicates_group_size <- min(table(group))
keep <- rowSums(cpm(y$counts)>5) >= smallest_replicates_group_size
y <- y[keep, , keep.lib.sizes=FALSE]
# dim(y)
message(paste("Filtering lowly expressed genes using a cpm cutoff of 5 in at least", smallest_replicates_group_size,
              "samples :", nrow(y), "kept,", nrow(raw_data)-nrow(y), "lost"))

### recompute library size
y$samples$lib.size <- colSums(y$counts)

### TMM normalisation
y <- calcNormFactors(y)

### Estimating the dispersion
y <- estimateDisp(y, design)

### Plot MDS, BCV and mean_var
glMDSPlot(y,labels=colnames(y$counts), groups=group, 
          main=paste(GOI_label,"MDS-Plot",sep="-",collapse=""), 
          html=paste(GOI_label,"MDS-Plot",sep="-",collapse=""), launch=FALSE)
# png("multi_dimensional_scaling_plot.png")
# plotMDS(y)
# invisible(dev.off())

png("mean_var.png");
plotMeanVar(y, show.tagwise.vars = TRUE, NBline = TRUE)
invisible(dev.off())

png("bcv.png");
plotBCV(y)
invisible(dev.off())

### Normalized counts raw
norm_raw <- y$counts
o <- order(rownames(norm_raw))
norm_raw <- norm_raw[o,]
write.csv(cbind(y$genes[rownames(norm_raw),1],norm_raw), "gene_normalised_expression_data_raw.csv", row.names=FALSE)
### Additionally also create files for updated GSEA analysis
temp_geneids <- as.data.frame(y$genes[rownames(norm_raw),1])
colnames(temp_geneids) <- "NAME"
temp_geneids$DESCRIPTION <- "na"
rev_norm_raw <- norm_raw[,order(ncol(norm_raw):1)]
write.table(cbind(temp_geneids, norm_raw), "gene_normalised_expression_data_raw_gsea.txt", sep='\t', row.names=FALSE, quote=FALSE)
writeLines(c(paste(length(design$hi), "2", "1", sep=" "), paste("#", "low", "high", sep=" "), paste0(design$hi, collapse=" ")), "gene_normalised_expression_data_raw_gsea.cls")

### Normalized CPM
norm_cpm <- cpm(y, normalized.lib.sizes=TRUE)
o <- order(rownames(norm_cpm))
norm_cpm <- norm_cpm[o,]
write.csv(cbind(y$genes[rownames(norm_cpm),1],norm_cpm), "gene_normalised_expression_data_cpm.csv", row.names=FALSE)

### log2 Normalized CPM
log_norm_cpm <- cpm(y, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1)
o <- order(rownames(log_norm_cpm))
log_norm_cpm <- log_norm_cpm[o,]
write.csv(cbind(y$genes[rownames(log_norm_cpm),1],log_norm_cpm), "gene_normalised_expression_data_logcpm.csv", row.names=FALSE)

###*****************************************************************************
#### Do DE test (exact) ####
###*****************************************************************************
### first check if switch comparison requested 
if (switch == "true") {
  comp <- c(hi, lo)
} else {
  comp <- c(lo, hi)
}

de <- exactTest(y, pair=comp)
dt <- decideTests(de,lfc=1)
message(paste("Differential expression using exact test performed", ":", 
              summary(dt)[[1]], "genes down regulated and", 
              summary(dt)[[3]], "genes up regulated"))

### save toptags and create .rnk file for GSEA
tt <- topTags(de, n=nrow(y))
write.csv(tt, "top_tags.csv", row.names=FALSE)
tt_rank <- tt$table
tt_rank["Rank"] <- (sign(tt_rank$logFC)) * (-log(tt_rank$FDR,10))
tt_rank <- tt_rank[c("gene_name", "Rank")]
tt_rank <- tt_rank[order(tt_rank$Rank), ]
tt_rank$Rank[tt_rank$Rank == Inf] = 320
tt_rank$Rank[tt_rank$Rank == -Inf] = -320
write.table(tt_rank, "top_tags_ranked.rnk", sep='\t', row.names=FALSE, quote=FALSE)

### make a deg_sigFC_table containing all the genes that are significantly up or downregulated
sorted_table <- tt$table[order(tt$table[,2]),]
rn <- rownames(tt$table)
deg_sigFC <- rn[(tt$table$FDR <= .05) & ((tt$table[,2] <= -1) | (tt$table[,2] >= 1))]
deg_sigFC_table <- sorted_table[deg_sigFC, ]

### Do the volcano plots - html versions with Glimma
results_df <- tt$table
results_df <- results_df[rownames(de$genes),]

results_df$PValue[results_df$PValue == 0.000000e+00] <- 1.0e-322
results_df$FDR[results_df$FDR == 0.000000e+00] <- 1.0e-320

if((summary(dt)[1,1] == 0) & (summary(dt)[3,1] == 0)){
	DE <- c("notDE")[as.factor(dt)]
	col_list = c("black")
} else if ((summary(dt)[1,1] > 0) & (summary(dt)[3,1] > 0)){
	DE <- c("downregulated", "notDE", "upregulated")[as.factor(dt)]
	col_list = c("blue","black","red")
} else if ((summary(dt)[1,1] > 0) & (summary(dt)[3,1] == 0)){
	DE <- c("downregulated", "notDE")[as.factor(dt)]
	col_list = c("blue","black")
} else if ((summary(dt)[1,1] == 0) & (summary(dt)[3,1] > 0)){
	DE <- c("notDE", "upregulated")[as.factor(dt)]
	col_list = c("black","red")
}
anno <- as.data.frame(cbind(dt,de$genes,DE))
anno <- cbind(rownames(anno), anno)
colnames(anno)<-c("EntrezID","dt", "GeneName", "DE")

with(results_df, glXYPlot(logFC, -log10(PValue), counts=cpm(y$counts,log=TRUE), 
      groups=group, samples=rownames(y$samples), main="High_Vs_Low", xlab="log2FC", 
      ylab="neg.log10.pValue", pch=20, cex=0.25, side.main="GeneName", 
      display.columns=c("EntrezID","GeneName","DE"), status=dt,anno=anno, 
      cols=col_list, html=paste(GOI_label,"High_Vs_Low-Volcano", sep="-",collapse=""), 
      launch=FALSE))

### Do the volcano plots - svg versions
svg(paste("High_Vs_Low", GOI_label,"volcano.svg",sep="_",collapse=""))
par(mar=c(5,6,4,2))
with(results_df, plot(logFC, -log10(PValue), pch=20, cex=0.25, col="grey", main="Volcano plot", cex.main=2.5, cex.lab=2.5, cex.axis=2.0))

### Add colored points: red if FDR<0.05, orange of log2FC>1, green if both)
with(subset(results_df, FDR < 0.05), points(logFC, -log10(PValue), pch=20, cex=0.25, col="green"))
with(subset(results_df, abs(logFC) >1 ), points(logFC, -log10(PValue), pch=20, cex=0.25, col="orange"))
with(subset(results_df, FDR < 0.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, cex=0.5, col="green"))
with(subset(results_df, FDR < 0.05 & logFC > 1), points(logFC, -log10(PValue), pch=20, cex=0.5, col="red"))
with(subset(results_df, FDR < 0.05 & logFC < -1), points(logFC, -log10(PValue), pch=20, cex=0.5, col="blue"))

### Label points with the textxy function from the calibrate plot
with(subset(results_df, -log10(results_df$PValue)>200 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_name, cex=.7,offset=0.5))

### Adding cut-off lines
FDR_sig_df <- subset(results_df, FDR < 0.05 & abs(logFC)>1)
write.csv(FDR_sig_df, file=paste("High_Vs_Low", GOI_label,"de_sigFC.csv",sep="_",collapse=""), row.names=FALSE)
yaxis_cuttoff <- -log10(FDR_sig_df[which.max(FDR_sig_df[,"FDR"]),"PValue"])
segments(-10, yaxis_cuttoff, -1, yaxis_cuttoff, col="grey", lty=3)
segments(1, yaxis_cuttoff, 10, yaxis_cuttoff, col="grey", lty=3)
segments(-1, yaxis_cuttoff, -1, 200, col="grey", lty=3)
segments(1, yaxis_cuttoff, 1, 200, col="grey", lty=3)

### Adding counts for sig up and down
total_reg_count <- nrow(FDR_sig_df)
up_reg_count <- nrow(subset(FDR_sig_df, logFC > 1))
down_reg_count <- nrow(subset(FDR_sig_df, logFC < -1))
mtext_string <- paste("# Sig DE:", total_reg_count, " (Up:", up_reg_count, ", Down:", down_reg_count, ")")
mtext(mtext_string, cex=1.0)
invisible(dev.off())

### Do Heatmap
logcounts_cpm <- cpm(y, normalized.lib.sizes = TRUE, log = TRUE)
logcounts_cpm_meansubtracted <- logcounts_cpm - rowMeans(logcounts_cpm)

d_table<-deg_sigFC_table

### Only keep columns relevant to the contrast being considered
cols<-c(1:ncol(logcounts_cpm))
DE_normalized_counts<-cbind(d_table,logcounts_cpm_meansubtracted[match(rownames(d_table),rownames(logcounts_cpm_meansubtracted)),cols])
DE_logcounts<-logcounts_cpm[match(rownames(d_table),rownames(logcounts_cpm)),cols]
DE_logcounts<-cbind(y$genes[rownames(d_table),],DE_logcounts)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# svg(filename=paste("High_Vs_Low", GOI, "heatmap.svg", sep="_"), width=10, height=10)
png(filename=paste("High_Vs_Low", GOI_label, "heatmap.png", sep="_"))
h<-heatmap.2(sapply(DE_normalized_counts[,(6:ncol(DE_normalized_counts))], as.numeric),
             hclustfun = function(x) hclust(x, method="single"),
             distfun=function(x) dist(x,method ='manhattan'),
             col=rev(morecols(50)), dendrogram="row", scale="row", trace="none",
             Colv = FALSE, margins=c(10,10), labRow = FALSE, labCol = FALSE, 
             ColSideColors = c(rep("blue", sum(design$lo)), rep("red", sum(design$hi))),
             key.title = NA,  
             main = paste("High Vs Low: ", GOI))
invisible(dev.off())

### This matrix should tally with the heatmap - the values in it are rowmean subtracted log counts.
write.csv(DE_logcounts[rev(h$rowInd),], file=paste("High_Vs_Low", GOI_label,"heatmap_scores.csv",sep="_"), row.names=FALSE)

### All done for DE so sign off
message(paste("Finished DE Analaysis for", GOI))

###*****************************************************************************
### Do GSEA / WebGestalt ####
###*****************************************************************************
### check that enrichment analyses requested
if (gsea_exe != "false"){
  gsea_exe <- paste(gsea_exe, "GSEAPreranked", sep =" ")
  
  message(paste("Starting GSEA Analysis for", GOI))
  
  main_dir <- getwd()
  dir.create(file.path(main_dir, "GSEA"))
  outdir <- paste(main_dir, "GSEA", sep="/")
  
  if (rna_species == "mRNA"){
    ### Hallmark
    message(paste("Passing DE results to GSEA using MSigDB Hallmark collection:", GOI))
    gsea_gmx <- "--gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v2023.1.Hs.symbols.gmt"
    gsea_rnk <- "-rnk top_tags_ranked.rnk"
    # gsea_res <- "-res gene_normalised_expression_data_raw_gsea.txt" 
    # gsea_cls <- "-cls gene_normalised_expression_data_raw_gsea.cls"
    gsea_chip <- "-chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_HGNC_ID_MSigDB.v7.4.chip"
    gsea_label <- paste("-rpt_label ", GOI_label, "_Strat_Vs_Hallmark", sep="")
    gsea_params1 <- "-collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -scoring_scheme weighted" 
    gsea_params2 <- "-create_svgs true -include_only_symbols true -make_sets true -plot_top_x 20"
    gsea_params3 <- "-rnd_seed timestamp -set_max 2000 -set_min 15 -zip_report false"
    gsea_out <- paste("-out ", outdir, sep="")
    gsea_command <- paste(gsea_exe, gsea_gmx, gsea_rnk, gsea_chip, gsea_label, gsea_params1, gsea_params2, gsea_params3, gsea_out, sep=" ") 
    
    system(gsea_command, ignore.stdout=TRUE, ignore.stderr=TRUE, wait=TRUE)
    
    ### Curated
    message(paste("Passing DE results to GSEA using MSigDB Curated collection:", GOI))
    gsea_gmx <- "-gmx ftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v2023.1.Hs.symbols.gmt"
    gsea_label <- paste("-rpt_label ", GOI_label, "_Strat_Vs_Curated", sep="")
    gsea_command <- paste(gsea_exe, gsea_gmx, gsea_rnk, gsea_chip, gsea_label, gsea_params1, gsea_params2, gsea_params3, gsea_out, sep=" ")
    
    system(gsea_command, ignore.stdout=TRUE, ignore.stderr=TRUE, wait=TRUE)
  }
  
  ### Create summary plots
  Sys.sleep(10)
  ### get dir and file names and process iteratively
  GSEA_output_dirs <- list.dirs(outdir, full.names=TRUE, recursive=FALSE)
  for (GSEA_output_dir in GSEA_output_dirs){
    combined_gsea_sig_results <- data.frame()
    GSEA_output_report_files <- list.files(path=GSEA_output_dir, pattern="gsea_report.*tsv", full.names=TRUE)
    for (GSEA_output_report_file in GSEA_output_report_files){
      sample_name <- strsplit(basename(GSEA_output_dir),"[.]")[[1]][1]
      gsea_report <- read.csv(GSEA_output_report_file, sep="\t", row.names=1)
      gsea_report <- Filter(function(x)!all(is.na(x)), gsea_report)
      temp_gsea_sig_results <- gsea_report[gsea_report$FDR.q.val < 0.05, ] # & abs(gsea_report$NES) > 2
      temp_gsea_sig_results <- head(temp_gsea_sig_results, n=10)
      combined_gsea_sig_results <- rbind(combined_gsea_sig_results, temp_gsea_sig_results)
      combined_gsea_sig_results <- combined_gsea_sig_results[order(combined_gsea_sig_results$NES),]
      if (GSEA_output_report_file == GSEA_output_report_files[length(GSEA_output_report_files)]){
        ###### Changed to detect zero length results
        if (nrow(combined_gsea_sig_results) > 0){
          ### all reports merged so create the plot
          combined_gsea_sig_results$NES <- as.numeric(combined_gsea_sig_results$NES)
          combined_gsea_sig_results$bar_color <- ifelse(combined_gsea_sig_results$NES > 0, "red", "blue")
          pdf(file=paste(outdir, "/",sample_name, ".pdf", sep=""),width=15,height=10)
          par(las=2) ; par(mar=c(10,50,10,5))
          barplot(combined_gsea_sig_results$NES, main=sample_name, horiz=TRUE, names.arg=row.names(combined_gsea_sig_results), xlab="Normalise Enrichment Score", cex.main=2, cex.lab=1.5, cex.axis=1.0, col=combined_gsea_sig_results$bar_color)
          invisible(dev.off())
          write.csv(combined_gsea_sig_results[,1:(length(combined_gsea_sig_results)-2)], 
                    paste(outdir, "/",sample_name, ".csv", sep=""), row.names=FALSE)
        } else {
          message("GSEA returned no significant hits. No summary plots produced")
        }  
        #########
      }
    }
  }
  
  ### All done for GSEA
  message(paste("Finished GSEA Analysis for", GOI))
  
  ### Do WebGestalt
  message(paste("Starting WebGestalt Analysis for", GOI))
  p_load(WebGestaltR)
  
  main_dir <- getwd()
  dir.create(file.path(main_dir, "WebGestalt"))
  WG_outdir <- paste(main_dir, "WebGestalt", sep="/")
  
  ### setup gene lists
  # subset all sig-up and sig-down regulated genes. If more than 500, take only top 500 
  FDR_sig_df_up <- FDR_sig_df[FDR_sig_df$logFC > 0,]
  FDR_sig_df_down <- FDR_sig_df[FDR_sig_df$logFC < 0,]
  FDR_sig_df_ord_up <- FDR_sig_df_up[order(FDR_sig_df_up$logFC, decreasing = TRUE), ]$gene_name
  FDR_sig_df_ord_down <- FDR_sig_df_down[order(FDR_sig_df_down$logFC), ]$gene_name
  if (length(FDR_sig_df_ord_up) > 500) { FDR_sig_df_ord_up <- FDR_sig_df_ord_up[1:500] }
  if (length(FDR_sig_df_ord_down) > 500) { FDR_sig_df_ord_down <- FDR_sig_df_ord_down[1:500] }
  
  ### setup reference list
  reference_geneset <- tt$table$gene_name
  
  ### setup databases required
  enrichDatabaseGO <- c("geneontology_Biological_Process_noRedundant","geneontology_Molecular_Function_noRedundant")
  enrichDatabasePW <- c("pathway_KEGG","pathway_Reactome")
  enrichDatabaseDS <- c("disease_Disgenet","disease_GLAD4U","disease_OMIM")
  enrichDBs <- list(enrichDatabaseGO = enrichDatabaseGO, enrichDatabasePW = enrichDatabasePW, enrichDatabaseDS = enrichDatabaseDS)
  
  ### Run each analysis
  if (length(FDR_sig_df_ord_up) > 50) {
    message("Passing DE UP regulated geneset to WebGestalt ORA")
    for (i in 1:length(enrichDBs)) { 
      sink(nullfile()) # don't want console messages
      enrichResult <- WebGestaltR(enrichMethod="ORA", 
                                  organism="hsapiens", 
                                  enrichDatabase=enrichDBs[[i]], 
                                  interestGene=FDR_sig_df_ord_up, 
                                  interestGeneType="genesymbol", 
                                  #referenceGene=reference_geneset, 
                                  referenceSet="genome_protein-coding", 
                                  referenceGeneType="genesymbol", 
                                  minNum=5, maxNum=2000, sigMethod="top", reportNumr=40,  
                                  isOutput=TRUE, 
                                  outputDirectory=WG_outdir, 
                                  projectName=paste(names(enrichDBs[i]), "UP-Reg", sep = "_"))
      sink()
    }
  } else {
    message("Not enough genes up regulated to do ORA analysis. Requires > 50")
  }
  
  if (length(FDR_sig_df_ord_down) > 50) {
    message("Passing DE Down regulated geneset to WebGestalt ORA")
    for (i in 1:length(enrichDBs)) { 
      sink(nullfile())
      enrichResult <- WebGestaltR(enrichMethod="ORA", 
                                  organism="hsapiens", 
                                  enrichDatabase=enrichDBs[[i]], 
                                  interestGene=FDR_sig_df_ord_down, 
                                  interestGeneType="genesymbol", 
                                  #referenceGene=reference_geneset, 
                                  referenceSet="genome_protein-coding", 
                                  referenceGeneType="genesymbol", 
                                  isOutput=TRUE, 
                                  outputDirectory=WG_outdir, 
                                  projectName=paste(names(enrichDBs[i]), "DOWN-Reg", sep = "_"))
      sink(nullfile())
    }
  } else {
    message("Not enough genes down regulated to do ORA analysis. Requires > 50")
  }
  ### All done for WebGestalt. Sign off now
  message(paste("Finished WebGestalt Analysis for", GOI))
}

