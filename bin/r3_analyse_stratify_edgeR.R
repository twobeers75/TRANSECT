#!/usr/bin/env Rscript

###*****************************************************************************
### DO DE analysis on RECODE3 stratified data
### JToubia - January 2023
###*****************************************************************************

###*****************************************************************************
### Import libraries ####
###*****************************************************************************
suppressMessages(if (!require("pacman")) install.packages("pacman"))
p_load(edgeR, Glimma, dplyr, tidyr, data.table, tibble, 
       ggplot2, ggforce, calibrate, gplots, RColorBrewer,
       rlogging)

SetLogFile(base.file=NULL)

###*****************************************************************************
#### Hard coded system variables ####
###*****************************************************************************


###*****************************************************************************
#### Read in Args ####
###*****************************************************************************
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	stop("At least one argument must be supplied (GOI).n", call.=FALSE)
} else if (length(args)==1) {
	# default to
  args[2] = "mRNA"
  args[3] = 5
  args[4] = "false"
  args[5] = "SCA/REF_FILES"
  args[6] = "SCA/bin/GSEA/gsea-cli.sh GSEAPreranked" 
  args[7] = list.files("SCA/data/", "RECODE3/BREAST", "R3-BREAST-*_tpm-mRNA.tsv", full.names=TRUE)
  args[8] = list.files("SCA/data/", "RECODE3/BREAST", "R3-BREAST-*_count-mRNA.tsv", full.names=TRUE)
  args[9] = ""
  arg2[10] = ""
}

GOI <- args[1]
rna_species <- args[2]
percentile <- as.integer(args[3])
switchDE <- args[4]
ref_files_folder <- args[5]
gsea_exe <- args[6]
gene_fpkm_filename <- args[7]
gene_counts_filename <- args[8]
isomir_rpm_filename <- args[9]
normals <- ""

###*****************************************************************************
### Setup for analyses ####
###*****************************************************************************
main_dir <- getwd()
dir.create(file.path(main_dir, "DE_Analysis"))
setwd(file.path(main_dir, "DE_Analysis"))

multi_GOI_analysis <- FALSE
ratio_GOI_analysis <- FALSE
single_GOI_analysis <- FALSE

###*****************************************************************************
### Setup for stratification ####
###*****************************************************************************
message(paste("Beginning stratification of", GOI))

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

### check GOI, split if necessary and set flags
if (grepl( "+", GOI, fixed = TRUE)){
  GOI_label <- gsub("[+]", "-", GOI)
  ### check that there is no more than 5 GsOI
  GsOI_split <- unlist(strsplit(GOI, "+", fixed=TRUE))
  if (length(GsOI_split) > 5) {
    stop("Can only have upto 5 genes in this type of analysis", call=FALSE)
  }
  message(paste("Stratification based on the addition of multiple genes of interest:", 
                paste(GOI, collapse="+")))
  ### set flag for multiple GOIs
  multi_GOI_analysis <- TRUE
} else if (grepl( "%", GOI, fixed = TRUE)){
  GOI_label <- gsub("%", "-", GOI)
  ### check that there is exactly 2 GsOI
  GsOI_split <- unlist(strsplit(GOI, "%", fixed=TRUE))
  if (length(GsOI_split) != 2) {
    stop("Must have 2 and 2 only genes in this type of analysis", call=FALSE)
  }
  message(paste("Stratification based on a ratio from 2 genes of interest:",  
                paste(GOI, collapse="%")))
  ### set flag for ratio GOIs
  ratio_GOI_analysis <- TRUE
} else {
  ### must be a single gene strat analysis
  message(paste("Stratification based on a single gene of interest:", GOI))
  GOI_label <- GOI
  GsOI_split <- GOI
  single_GOI_analysis <- TRUE
}

### subset data for G/sOI only. We require a few copies for histogram and boxplots (including the original (OG))
GOI_exp_raw_OG <- as.data.frame(t(exp_data[GsOI_split, ]))
GOI_exp_raw_OG <- na.omit(GOI_exp_raw_OG)
GOI_exp_hist <- GOI_exp_raw_OG
GOI_exp_boxplot <- GOI_exp_raw_OG
### save the table here
write.table(round(GOI_exp_raw_OG,2), "GOI_exp_raw_OG.tsv", sep='\t', col.names=NA)

### check to see if multi/ratio flags set
if (multi_GOI_analysis){
  ### create a long copy of raw expression for boxplot
  GOI_exp_sep <- GOI_exp_raw_OG
  GOI_exp_sep$patientID <- rownames(GOI_exp_sep)
  GOI_exp_sep$patientID <- factor(GOI_exp_sep$patientID)
  rownames(GOI_exp_sep) <- NULL
  GOI_exp_sep_long <- gather(GOI_exp_sep, "geneID", "expression", all_of(GsOI_split), factor_key=TRUE)
  GOI_exp_sep_long$expression <- log2(GOI_exp_sep_long$expression)
  ### create scaled df
  GOI_exp_scaled <- data.frame(row.names = rownames(GOI_exp_raw_OG))
  for (gene in GsOI_split){
    col_id <- paste(gene,"rank", sep="_")
    # GOI_exp_scaled[col_id] <- ntile(GOI_exp_raw_OG[gene], nrow(GOI_exp_raw_OG))
    GOI_exp_scaled[col_id] <- min_rank(GOI_exp_raw_OG[gene])
  }
  ### average all ranks for ordering
  GOI_exp_scaled <- GOI_exp_scaled %>%
    transmute(average=rowMeans(across(where(is.numeric))))
  colnames(GOI_exp_scaled) <- "ranking_score"
  ### now sum all raw expression values for plotting
  GOI_exp_raw_sum <- GOI_exp_raw_OG %>%
    transmute(sum=rowSums(across(where(is.numeric))))
  colnames(GOI_exp_raw_sum) <- GOI
  ### now cbind on index
  GOI_exp_boxplot <- cbind(GOI_exp_raw_sum, GOI_exp_scaled)
  GOI_exp_hist <- cbind(GOI_exp_raw_OG, GOI_exp_raw_sum, GOI_exp_scaled)
}
if (ratio_GOI_analysis){
  GOI_exp_boxplot = GOI_exp_boxplot + 1
  ### calc log2 FC
  GOI_exp_boxplot$ranking_score = (log2(GOI_exp_boxplot[,GsOI_split[1]])) - (log2(GOI_exp_boxplot[,GsOI_split[2]]))
  ### now divide raw expression values for plotting
  GOI_exp_boxplot$div <- GOI_exp_boxplot[,GsOI_split[1]] / GOI_exp_boxplot[,GsOI_split[2]]
  ### reset GOI variable
  # GOI <- paste(GOI, collapse="%")
  GOI_exp_boxplot <- rename(GOI_exp_boxplot, !!GOI := div)
  GOI_exp_hist <- GOI_exp_boxplot
}

### check and strip normals out if requested
if( grepl(pattern = "^TCGA", rownames(GOI_exp_hist)[1]) ){
  if (normals == "") {
    message(paste("Excluding non-diseased samples"))
    GOI_exp_hist$cancertype <- substr(rownames(GOI_exp_hist), 14, 15)
    GOI_exp_hist <- GOI_exp_hist[GOI_exp_hist$cancertype != 11, ]
    GOI_exp_hist <- subset(GOI_exp_hist, select=-cancertype)
  }
}

if (rna_species == "miRNA"){
  GOI_exp_boxplot <- subset(GOI_exp_boxplot, rownames(GOI_exp_boxplot) %in% mRNA_ids)
  GOI_exp_hist <- GOI_exp_boxplot
}

if (multi_GOI_analysis | ratio_GOI_analysis){
  GOI_exp_hist <- GOI_exp_hist %>% arrange(ranking_score)
  GOI_exp_wstrat = mutate(GOI_exp_hist, quantile_rank = ntile(GOI_exp_hist["ranking_score"],4))
  GOI_exp_wstrat = mutate(GOI_exp_wstrat, percentile_rank = round(cume_dist(GOI_exp_hist["ranking_score"])*100,2))
} else {
  GOI_exp_hist <- GOI_exp_hist %>% arrange(!!as.symbol(GOI))
  GOI_exp_wstrat = mutate(GOI_exp_hist, quantile_rank = ntile(GOI_exp_hist[GOI],4))
  if (nrow(GOI_exp_hist) < 100){
    GOI_exp_wstrat <- mutate(GOI_exp_wstrat, percentile_rank = ntile(GOI_exp_hist[GOI],10))
    GOI_exp_wstrat$percentile_rank <- GOI_exp_wstrat$percentile_rank * 10
  } else {
    GOI_exp_wstrat <- mutate(GOI_exp_wstrat, percentile_rank = round(cume_dist(GOI_exp_hist[GOI])*100,2))
  }
}

### write results to file
rownames(GOI_exp_wstrat) <- rownames(GOI_exp_hist)
write.table(GOI_exp_wstrat, "GOI_with_strat.tsv", sep='\t', col.names=NA)

### now create design matrix
# first_quart_df <- data.frame(matrix(ncol=3, nrow=dim(GOI_exp_wstrat[GOI_exp_wstrat$quantile_rank == 1, ])[1]))
# first_quart_df$Sample <- rownames(GOI_exp_wstrat[GOI_exp_wstrat$quantile_rank == 1, ])
low_percentile_df <- data.frame(matrix(ncol=3, nrow=dim(GOI_exp_wstrat[GOI_exp_wstrat$percentile_rank <= percentile, ])[1]))
matrix_colnames <- c("Sample", "lo", "hi")
colnames(low_percentile_df) <- matrix_colnames
low_percentile_df$Sample <- rownames(GOI_exp_wstrat[GOI_exp_wstrat$percentile_rank <= percentile, ])
low_percentile_df$lo <- 1
low_percentile_df$hi <- 0
message(paste(nrow(low_percentile_df), "patients stratified into low group for", GOI))

# third_quart_df <- data.frame(matrix(ncol=3, nrow=dim(GOI_exp_wstrat[GOI_exp_wstrat$quantile_rank == 4, ])[1]))
# third_quart_df$Sample <- rownames(GOI_exp_wstrat[GOI_exp_wstrat$quantile_rank == 4, ])
high_percentile_df <- data.frame(matrix(ncol=3, nrow=dim(GOI_exp_wstrat[GOI_exp_wstrat$percentile_rank >= (100-percentile), ])[1]))
colnames(high_percentile_df) <- matrix_colnames
high_percentile_df$Sample <- rownames(GOI_exp_wstrat[GOI_exp_wstrat$percentile_rank >= (100-percentile), ])
high_percentile_df$lo <- 0
high_percentile_df$hi <- 1
message(paste(nrow(high_percentile_df), "patients stratified into high group for", GOI))

design <- rbind(low_percentile_df, high_percentile_df)
rownames(design) <- design$Sample
design <- design[c("lo", "hi")]
write.table(design, "design.tsv", col.names=NA)

### need to check that all GOIs are expressed at decent levels for DE analysis
## require that the median expression value > 3 AND
## require that the max value > 5
GsOI_medians <- GOI_exp_raw_OG %>% summarise_at(GsOI_split, median)
GsOI_maxs <- GOI_exp_raw_OG %>% summarise_at(GsOI_split, max)
if (min(GsOI_medians) < 2 || min(GsOI_maxs) < 5) {
  # set DE flag to false but continue on to produce descriptive plots 
  run_DE <- FALSE
} else {
  run_DE <- TRUE
}

###*****************************************************************************
### Do descriptive plots ####
###*****************************************************************************

### Do frequency histogram OR scatterplots
if (ratio_GOI_analysis){
  ### Do Scatter
  GOI_exp_scatter <- log2(GOI_exp_wstrat[,1:2])
  GOI_exp_scatter$rank <- 1:nrow(GOI_exp_scatter)
  GOI_exp_scatter$percentile_rank <- GOI_exp_wstrat$percentile_rank
  svg(paste(GOI_label, "TPM_Scatter.svg", sep="_"))
  par(mar=c(5,6,4,2))
  plot(GOI_exp_scatter[,3], GOI_exp_scatter[,1], main="Ordered expression ratio", 
       xlab="Expression Ratio", ylab="log2(TPM)", col=alpha("#F8766D", 0.5), pch=15, cex=0.8,
       cex.main=3, cex.lab=2.5, cex.axis=2.0, xaxt="none", ylim=c(0,15))
  points(GOI_exp_scatter[,2], col=alpha("#00BFC4", 0.5), pch=17, cex=0.8)
  abline(v=max(GOI_exp_scatter[GOI_exp_scatter$percentile_rank <= percentile, "rank"]), col="blue", lty=5, lwd=2)
  abline(v=min(GOI_exp_scatter[GOI_exp_scatter$percentile_rank >= (100-percentile), "rank"]), col="red", lty=5, lwd=2)
  legend("top", legend=c(colnames(GOI_exp_scatter[1]), colnames(GOI_exp_scatter[2])), 
         col=c("#F8766D", "#00BFC4"), pch=c(15,17), cex=1.5, bg="transparent", bty = "n")
  invisible(dev.off())
} 
if (multi_GOI_analysis){
  ### Do Boxplots with sina
  GOI_exp_sep_long$cols_column <- ifelse(GOI_exp_sep_long$patientID %in% rownames(design[design$lo == 1,]), "#00BFC4", 
     ifelse(GOI_exp_sep_long$patientID %in% rownames(design[design$hi == 1,]), "#F8766D", "grey"))
  GOI_exp_sep_long$size_column <- ifelse(GOI_exp_sep_long$patientID %in% rownames(design[design$lo == 1,]), 2, 
     ifelse(GOI_exp_sep_long$patientID %in% rownames(design[design$hi == 1,]), 2, 0.75))
  GOI_exp_sep_long$in_column <- GOI_exp_sep_long$patientID %in% rownames(design)
  GOI_exp_sep_long <- GOI_exp_sep_long[order(GOI_exp_sep_long$geneID, GOI_exp_sep_long$in_column), ]
  GOI_exp_sep_long <- GOI_exp_sep_long[!is.infinite(GOI_exp_sep_long$expression),]
  p <- ggplot(GOI_exp_sep_long, aes(geneID, expression))
  p + geom_sina(colour = GOI_exp_sep_long$cols_column, size = GOI_exp_sep_long$size_column, alpha = 0.7) + 
    theme_classic() + theme(text=element_text(size=30)) + ylab(paste("Log2(", exp_type, ")", sep=""))
  invisible(ggsave(paste(GOI_label, "TPM_Boxplot_Sina.svg", sep="_"), last_plot(), width = 10, height = 8))
}
if (single_GOI_analysis){
  ### Do Histogram
  svg(paste(GOI_label, "TPM_histogram.svg", sep="_"))
  par(mar=c(5,5.5,4,2))
  hist(log2(GOI_exp_hist[, GOI]), breaks=100, main=GOI, 
       xlab=paste("Log2(", exp_type, ")", sep=""), ylab="Frequency",
       cex.main=3, cex.lab=2.5, cex.axis=2.0)
  abline(v=mean(log2(GOI_exp_hist[, GOI])), col="purple")
  abline(v=median(log2(GOI_exp_hist[, GOI])), col="green")
  abline(v=max(log2(GOI_exp_wstrat[GOI_exp_wstrat$percentile_rank <= percentile, GOI])), col="blue")
  abline(v=min(log2(GOI_exp_wstrat[GOI_exp_wstrat$percentile_rank >= (100-percentile), GOI])), col="red")
  # abline(v=max(log2(GOI_exp_wstrat[GOI_exp_wstrat$quantile_rank == 1, GOI])), col="blue")
  # abline(v=min(log2(GOI_exp_wstrat[GOI_exp_wstrat$quantile_rank == 4, GOI])), col="red")
  invisible(dev.off())
}

### Do tumor - normal expression boxplot
if( grepl(pattern = "^TCGA", rownames(GOI_exp_boxplot)[1]) ){
  GOI_exp_boxplot$Normal <- grepl(pattern = "(-11A|-11B)$", rownames(GOI_exp_boxplot))
} else {
  GOI_exp_boxplot$Normal <- "TRUE" 
}
GOI_exp_boxplot$Normal <- as.factor(GOI_exp_boxplot$Normal)
GOI_exp_boxplot[, GOI] <- log2(GOI_exp_boxplot[, GOI])
svg(paste(GOI_label, "TPM_N-T_boxplot.svg", sep="_"))
p <- ggplot(GOI_exp_boxplot, aes(x=Normal, y=get(GOI), color=Normal)) + 
  geom_boxplot(notch=TRUE) + theme_classic() + 
  theme(text=element_text(size=30))
p + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  ylab(paste("Log2(", GOI, " ", exp_type, ")", sep=""))
invisible(dev.off())

### Do stratified expression boxplot
GOI_exp_low <- GOI_exp_wstrat[GOI_exp_wstrat$percentile_rank <= percentile, ]
GOI_exp_low$stratification <- "Low"
GOI_exp_high <- GOI_exp_wstrat[GOI_exp_wstrat$percentile_rank >= (100-percentile), ]
GOI_exp_high$stratification <- "High"
GOI_selected_strat <- rbind(GOI_exp_low, GOI_exp_high)
GOI_selected_strat[, GOI] <- log2(GOI_selected_strat[, GOI])
svg(paste(GOI_label, "TPM_strat_boxplot.svg", sep="_"))
p <- ggplot(GOI_selected_strat, aes(x=stratification, y=get(GOI), color=stratification)) + 
  geom_boxplot(notch=FALSE, na.rm=TRUE) + 
  ylab(paste("Log2(", GOI, " ", exp_type, ")", sep="")) + 
  theme_classic() + theme(text=element_text(size=30))
p + geom_jitter(na.rm=TRUE, shape=16, position=position_jitter(0.2)) + 
  scale_x_discrete(limits=c("Low", "High")) 
# p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth=0.1) + 
invisible(dev.off())

### all done here, sign off stratification
message(paste("Finished stratification of", GOI))

###*****************************************************************************
### Setup for DE ####
###*****************************************************************************
if (!run_DE) {
  warning(paste("Unfortunately, there is not enough observations in this dataset for", GOI))
  stop("Ending this process") 
}

message(paste("Beginning DE Analysis for", GOI))

hi <- 'hi'
lo <- 'lo'

### Read in data
raw_data <- fread(file=gene_counts_filename, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
colnames(raw_data)[1] <- 'gene_name'
message(paste("Raw count expression data loaded for", GOI, ":", ncol(raw_data)-1, "patients and", nrow(raw_data), "genes"))

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
if (switchDE == "true") {
  comp <- c(hi, lo)
} else {
  comp <- c(lo, hi)
}

de <- exactTest(y, pair=comp)
dt <- decideTests(de,adjust.method="fdr",p.value=1e-5,lfc=1)
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
deg_sigFC <- rn[(tt$table$FDR <= 1e-5) & ((tt$table[,2] <= -1) | (tt$table[,2] >= 1))]
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

### Add colored points
with(subset(results_df, FDR < 1e-5), points(logFC, -log10(PValue), pch=20, cex=0.25, col="green"))
with(subset(results_df, abs(logFC) >1 ), points(logFC, -log10(PValue), pch=20, cex=0.25, col="orange"))
with(subset(results_df, FDR < 1e-5 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, cex=0.5, col="green"))
with(subset(results_df, FDR < 1e-5 & logFC > 1), points(logFC, -log10(PValue), pch=20, cex=0.5, col="red"))
with(subset(results_df, FDR < 1e-5 & logFC < -1), points(logFC, -log10(PValue), pch=20, cex=0.5, col="blue"))

### Label points with the textxy function from the calibrate plot
with(subset(results_df, -log10(results_df$PValue)>200 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=gene_name, cex=.7,offset=0.5))

### Adding cut-off lines
FDR_sig_df <- subset(results_df, FDR < 1e-5 & abs(logFC)>1)
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
  
  if (rna_species == "miRNA"){
    message(paste("Passing DE results to GSEA using Custom TargetScan collection:", GOI))
    gmx_file <- paste(ref_files_folder, "CSCS_Hs_8mer_CSfiltered.gmx", sep="/")
    ### In-house Targetscan collection
    gsea_gmx <- paste("-gmx", gmx_file, sep=" ")
    gsea_rnk <- "-rnk top_tags_ranked.rnk" 
    gsea_chip <- "-chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_HGNC_ID_MSigDB.v7.4.chip"
    gsea_label <- paste("-rpt_label ", GOI_label, "_Strat_Vs_Inhouse-TargetScan", sep="")
    gsea_params1 <- "-collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -scoring_scheme weighted" 
    gsea_params2 <- "-create_svgs true -include_only_symbols true -make_sets true -plot_top_x 20"
    gsea_params3 <- "-rnd_seed timestamp -set_max 2000 -set_min 15 -zip_report false"
    gsea_out <- paste("-out ", outdir, sep="")
    gsea_command <- paste(gsea_exe, gsea_gmx, gsea_rnk, gsea_chip, gsea_label, gsea_params1, gsea_params2, gsea_params3, gsea_out, sep=" ") 
    
    system(gsea_command, ignore.stdout=TRUE, ignore.stderr=TRUE, wait=TRUE)
  }
  
  # C3 not in use: ftp.broadinstitute.org://pub/gsea/gene_sets/c3.all.v7.4.symbols.gmt
  
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
        ### all reports merged so create the plot
        combined_gsea_sig_results$NES <- as.numeric(combined_gsea_sig_results$NES)
        combined_gsea_sig_results$bar_color <- ifelse(combined_gsea_sig_results$NES > 0, "red", "blue")
        pdf(file=paste(outdir, "/",sample_name, ".pdf", sep=""),width=15,height=10)
        par(las=2) ; par(mar=c(10,50,10,5))
        barplot(combined_gsea_sig_results$NES, main=sample_name, horiz=TRUE, names.arg=row.names(combined_gsea_sig_results), xlab="Normalise Enrichment Score", cex.main=2, cex.lab=1.5, cex.axis=1.0, col=combined_gsea_sig_results$bar_color)
        invisible(dev.off())
        write.csv(combined_gsea_sig_results[,1:(length(combined_gsea_sig_results)-2)], 
                  paste(outdir, "/",sample_name, ".csv", sep=""), row.names=FALSE)
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
  message("Passing DE UP regulated geneset to WebGestalt ORA")
  sink(nullfile()) # don't want console messages
  for (i in 1:length(enrichDBs)) { 
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
  }
  message("Passing DE Down regulated geneset to WebGestalt ORA")
  for (i in 1:length(enrichDBs)) { 
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
  }
  sink()
  ### All done for WebGestalt. Sign off now
  message(paste("Finished WebGestalt Analysis for", GOI))
}
