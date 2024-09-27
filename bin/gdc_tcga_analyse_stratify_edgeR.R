#!/usr/bin/env Rscript

###*****************************************************************************
### DO DE analysis on TCGA stratified data
### JToubia - January 2021
###*****************************************************************************

###*****************************************************************************
# Import libraries ####
###*****************************************************************************
suppressMessages(if (!require("pacman")) install.packages("pacman"))
p_load(edgeR, Glimma, dplyr, tidyr, data.table, tibble, 
       ggplot2, ggforce, calibrate, gplots, RColorBrewer,
       rlogging, plotly, htmlwidgets)

SetLogFile(base.file=NULL)
options(warn=-1)

###*****************************************************************************
# Hard coded system variables ####
###*****************************************************************************


###*****************************************************************************
# Read in Args ####
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
  args[7] = list.files("SCA/data/GDC/", "BRCA/mRNA_expression_fpkm", "GDC_TCGA-BRCA.*FPKM-mRNA_toTPM_all.tsv", full.names=TRUE)
  args[8] = list.files("SCA/data/GDC/", "BRCA/mRNA_expression_counts", "GDC_TCGA-BRCA.*Count-mRNA_all.tsv", full.names=TRUE)
  args[9] = list.files("SCA/data/GDC/", "BRCA/isomiR_expression_rpm", "GDC_TCGA-BRCA.*RPM-miRNAisoform_all.tsv", full.names=TRUE)
  args[10] = ""
}

GOI <- args[1]
rna_species <- args[2]
percentile <- as.numeric(args[3])
switchDE <- args[4]
ref_files_folder <- args[5]
gsea_exe <- args[6]
gene_fpkm_filename <- args[7]
gene_counts_filename <- args[8]
isomir_rpm_filename <- args[9]
normals <- ""

###*****************************************************************************
# Setup ####
###*****************************************************************************
main_dir <- getwd()
dir.create(file.path(main_dir, "DE_Analysis"))
setwd(file.path(main_dir, "DE_Analysis"))

###*****************************************************************************
## Setup custom thresholds ####
###*****************************************************************************
# remove GOI(s) from DE analyses
remove_GsOI <- FALSE
# filter lowly expressed genes threshold - require CPM > threshold in half the samples
low_gene_thrs <- 5
# decideTests pvalue and lfc
dt_pvalue <- 1e-5
dt_lfc <- 1
# min/max sig genes to use for ORA analysis
ORA_min <- 20
ORA_max <- 500

strat_do_multi_GOI_analysis <- FALSE
strat_do_ratio_GOI_analysis <- FALSE
strat_do_single_GOI_analysis <- FALSE

###*****************************************************************************
# Do stratification ####
###*****************************************************************************
message(paste("Beginning stratification of", GOI))

### Check what species of RNA we are analysing and set/read appropriate variables
if (rna_species == "mRNA"){
  strat_exp_data <- fread(file=gene_fpkm_filename, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
  strat_exp_data <- strat_exp_data %>% column_to_rownames('gene_name')
  strat_exp_data <- strat_exp_data + 0.01 # adding 0.1 for log calculations
  strat_exp_type <- "TPM"
} else if (rna_species == "miRNA"){
  mRNA_ids <- t(fread(file=gene_fpkm_filename, header = FALSE, nrows = 1))
  mRNA_ids <- mRNA_ids[-1]
  strat_exp_data <- fread(file=isomir_rpm_filename, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
  strat_exp_data <- strat_exp_data %>% column_to_rownames('miR_name')
  strat_exp_type <- "RPM"
} else {
  stop("Can only stratify by mRNA or miRNA", call=FALSE)
}
message(paste("Normalised expression data loaded for", GOI, ":", ncol(strat_exp_data), "patients and", nrow(strat_exp_data), "genes"))

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
  strat_do_multi_GOI_analysis <- TRUE
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
  strat_do_ratio_GOI_analysis <- TRUE
} else {
  ### must be a single gene strat analysis
  message(paste("Stratification based on a single gene of interest:", GOI))
  GOI_label <- GOI
  GsOI_split <- GOI
  strat_do_single_GOI_analysis <- TRUE
}

### subset data for G/sOI only. We require a few copies for histogram and boxplots (including the original (OG))
tmp_col_num <- ncol(strat_exp_data)
strat_GOI_exp_raw_OG <- as.data.frame(t(strat_exp_data[GsOI_split, ][,1:tmp_col_num]))
strat_GOI_exp_raw_OG <- na.omit(strat_GOI_exp_raw_OG) # omit all rows with na
strat_GOI_exp_hist <- strat_GOI_exp_raw_OG
strat_GOI_exp_boxplot <- strat_GOI_exp_raw_OG
### save the table here
write.table(round(strat_GOI_exp_raw_OG,2), "GOI_exp_raw_OG.tsv", sep='\t', col.names=NA)

### check to see if multi/ratio flags set
if (strat_do_multi_GOI_analysis){
  ### create a long copy of raw expression for boxplot
  strat_GOI_exp_sep <- strat_GOI_exp_raw_OG
  strat_GOI_exp_sep$patientID <- rownames(strat_GOI_exp_sep)
  strat_GOI_exp_sep$patientID <- factor(strat_GOI_exp_sep$patientID)
  rownames(strat_GOI_exp_sep) <- NULL
  strat_GOI_exp_sep_long <- gather(strat_GOI_exp_sep, "geneID", "expression", all_of(GsOI_split), factor_key=TRUE)
  strat_GOI_exp_sep_long$expression <- log2(strat_GOI_exp_sep_long$expression)
  ### get sum rank positions
  strat_GOI_exp_scaled <- data.frame(row.names = rownames(strat_GOI_exp_raw_OG))
  for (gene in GsOI_split){
    col_id <- paste(gene,"rank", sep="_")
    # strat_GOI_exp_scaled[col_id] <- ntile(strat_GOI_exp_raw_OG[gene], nrow(strat_GOI_exp_raw_OG))
    strat_GOI_exp_scaled[col_id] <- min_rank(strat_GOI_exp_raw_OG[gene])
  }
  ### average all ranks for ordering
  strat_GOI_exp_scaled <- strat_GOI_exp_scaled %>%
    transmute(average=rowMeans(across(where(is.numeric))))
  colnames(strat_GOI_exp_scaled) <- "ranking_score"
  ### now sum all raw expression values for plotting
  strat_GOI_exp_raw_sum <- strat_GOI_exp_raw_OG %>%
    transmute(sum=rowSums(across(where(is.numeric))))
  colnames(strat_GOI_exp_raw_sum) <- GOI
  ### now cbind on index
  strat_GOI_exp_boxplot <- cbind(strat_GOI_exp_raw_sum, strat_GOI_exp_scaled)
  strat_GOI_exp_hist <- cbind(strat_GOI_exp_raw_OG, strat_GOI_exp_raw_sum, strat_GOI_exp_scaled)
}
if (strat_do_ratio_GOI_analysis){
  strat_GOI_exp_boxplot = strat_GOI_exp_boxplot + 1
  ### calc log2 FC
  strat_GOI_exp_boxplot$ranking_score = (log2(strat_GOI_exp_boxplot[,GsOI_split[1]])) - (log2(strat_GOI_exp_boxplot[,GsOI_split[2]]))
  ### now divide raw expression values for plotting
  strat_GOI_exp_boxplot$div <- strat_GOI_exp_boxplot[,GsOI_split[1]] / strat_GOI_exp_boxplot[,GsOI_split[2]]
  strat_GOI_exp_boxplot <- rename(strat_GOI_exp_boxplot, !!GOI := div)
  strat_GOI_exp_hist <- strat_GOI_exp_boxplot
}

if (rna_species == "miRNA"){
  strat_GOI_exp_boxplot <- subset(strat_GOI_exp_boxplot, rownames(strat_GOI_exp_boxplot) %in% mRNA_ids)
  strat_GOI_exp_hist <- strat_GOI_exp_boxplot
}

### check and strip normals out if requested
if (normals == "") {
  message(paste("Excluding non-diseased  samples"))
  strat_GOI_exp_hist$cancertype <- substr(rownames(strat_GOI_exp_hist), 14, 15)
  strat_GOI_exp_hist <- strat_GOI_exp_hist[strat_GOI_exp_hist$cancertype != 11, ]
  strat_GOI_exp_hist <- subset(strat_GOI_exp_hist, select=-cancertype)
}

if (strat_do_multi_GOI_analysis | strat_do_ratio_GOI_analysis){
  strat_GOI_exp_hist <- strat_GOI_exp_hist %>% arrange(ranking_score)
  strat_GOI_exp_wstrat <- mutate(strat_GOI_exp_hist, quantile_rank = ntile(strat_GOI_exp_hist["ranking_score"],4))
  strat_GOI_exp_wstrat <- mutate(strat_GOI_exp_wstrat, percentile_rank = round(cume_dist(strat_GOI_exp_hist["ranking_score"])*100,2))
} else {
  strat_GOI_exp_hist <- strat_GOI_exp_hist %>% arrange(!!as.symbol(GOI))
  strat_GOI_exp_wstrat <- mutate(strat_GOI_exp_hist, quantile_rank = ntile(strat_GOI_exp_hist[GOI],4))
  if (nrow(strat_GOI_exp_hist) < 100){
    strat_GOI_exp_wstrat <- mutate(strat_GOI_exp_wstrat, percentile_rank = ntile(strat_GOI_exp_hist[GOI],10))
    strat_GOI_exp_wstrat$percentile_rank <- strat_GOI_exp_wstrat$percentile_rank * 10
  } else {
    strat_GOI_exp_wstrat <- mutate(strat_GOI_exp_wstrat, percentile_rank = round(cume_dist(strat_GOI_exp_hist[GOI])*100,2))
  }
}

### write results to file
rownames(strat_GOI_exp_wstrat) <- rownames(strat_GOI_exp_hist)
write.table(strat_GOI_exp_wstrat, "GOI_exp_with_strat.tsv", sep='\t', col.names=NA)

### need to check that all GOIs are expressed at decent levels for DE analysis
## require that the median expression value > 1 AND
## require that the max value > 5
## ie. Expression of GOI >1 in at least half the cohort AND at least one individual with expression >5
GsOI_medians <- strat_GOI_exp_raw_OG %>% summarise_at(GsOI_split, median)
GsOI_maxs <- strat_GOI_exp_raw_OG %>% summarise_at(GsOI_split, max)
if (min(GsOI_medians) < 1 || min(GsOI_maxs) < 5) {
  # set DE flag to false but continue on to produce descriptive plots 
  run_DE <- FALSE
} else {
  run_DE <- TRUE
}

if (run_DE) {
  ### now create design matrix
  # first_quart_df <- data.frame(matrix(ncol=3, nrow=dim(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$quantile_rank == 1, ])[1]))
  # first_quart_df$Sample <- rownames(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$quantile_rank == 1, ])
  strat_low_percentile_df <- data.frame(matrix(ncol=3, nrow=dim(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank <= percentile, ])[1]))
  strat_matrix_colnames <- c("Sample", "lo", "hi")
  colnames(strat_low_percentile_df) <- strat_matrix_colnames
  strat_low_percentile_df$Sample <- rownames(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank <= percentile, ])
  strat_low_percentile_df$lo <- 1
  strat_low_percentile_df$hi <- 0
  message(paste(nrow(strat_low_percentile_df), "patients stratified into low group for", GOI))
  
  # third_quart_df <- data.frame(matrix(ncol=3, nrow=dim(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$quantile_rank == 4, ])[1]))
  # third_quart_df$Sample <- rownames(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$quantile_rank == 4, ])
  strat_high_percentile_df <- data.frame(matrix(ncol=3, nrow=dim(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank >= (100-percentile), ])[1]))
  colnames(strat_high_percentile_df) <- strat_matrix_colnames
  strat_high_percentile_df$Sample <- rownames(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank >= (100-percentile), ])
  strat_high_percentile_df$lo <- 0
  strat_high_percentile_df$hi <- 1
  message(paste(nrow(strat_high_percentile_df), "patients stratified into high group for", GOI))
  
  design <- rbind(strat_low_percentile_df, strat_high_percentile_df)
  rownames(design) <- design$Sample
  design <- design[c("lo", "hi")]
  write.table(design, "design.tsv", col.names=NA)
}

###*****************************************************************************
## Do descriptive plots ####
###*****************************************************************************

### Do frequency histogram OR scatterplots
if (strat_do_ratio_GOI_analysis){
  ### Do Scatter
  strat_GOI_exp_scatter <- log2(strat_GOI_exp_wstrat[,1:2])
  strat_GOI_exp_scatter$rank <- 1:nrow(strat_GOI_exp_scatter)
  strat_GOI_exp_scatter$percentile_rank <- strat_GOI_exp_wstrat$percentile_rank
  min_bound <- max(strat_GOI_exp_scatter[strat_GOI_exp_scatter$percentile_rank <= percentile, "rank"])
  max_bound <- min(strat_GOI_exp_scatter[strat_GOI_exp_scatter$percentile_rank >= (100-percentile), "rank"])
  
  strat_GOI_exp_scatter <- log2(strat_GOI_exp_wstrat[GsOI_split])
  strat_GOI_exp_scatter <- rownames_to_column(strat_GOI_exp_scatter, var = "ID")
  strat_GOI_exp_scatter$rank <- 1:nrow(strat_GOI_exp_scatter)
  strat_GOI_exp_scatter <- strat_GOI_exp_scatter %>% pivot_longer(!c(rank,ID), names_to = "gene", values_to = "expression")
  
  pscatter <- ggplot(strat_GOI_exp_scatter, aes(x=rank, y=expression, color=gene, shape=gene, labels=ID)) + 
    geom_point(alpha=0.5, size=2) +
    labs(title="Ordered expression ratio",x="Expression Ratio", y ="log2(TPM)") +
    theme_classic() + ylim(0,15) +
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=30), 
          axis.ticks.x=element_blank(), axis.text.x=element_blank(), 
          legend.title=element_blank(), legend.position.inside = c(0.5,0.95))
  pscatter <- pscatter + geom_vline(aes(xintercept=min_bound), color="blue", linewidth=0.5, linetype="dashed") + 
    geom_vline(aes(xintercept=max_bound), color="red", linewidth=0.5, linetype="dashed") 
  # invisible(ggsave(paste(GOI_label, "TPM_Scatter.svg", sep="_"), width = 7, height = 7))
  
  p_ly = ggplotly(pscatter)
  p_ly$x$data[[1]]$text <- sub(paste("gene: ",GsOI_split[1],"<br />", sep=""), "", p_ly$x$data[[1]]$text, fixed = TRUE)
  p_ly$x$data[[2]]$text <- sub(paste("gene: ",GsOI_split[2],"<br />", sep=""), "", p_ly$x$data[[2]]$text, fixed = TRUE)
  
  saveWidget(ggplotly(p_ly), file = paste(GOI_label, "TPM_Scatter.html", sep="_"))
  
  # strat_GOI_exp_scatter <- log2(strat_GOI_exp_wstrat[,1:2])
  # strat_GOI_exp_scatter$rank <- 1:nrow(strat_GOI_exp_scatter)
  # strat_GOI_exp_scatter$percentile_rank <- strat_GOI_exp_wstrat$percentile_rank
  # svg(paste(GOI_label, "TPM_Scatter.svg", sep="_"))
  # par(mar=c(5,6,4,2))
  # plot(strat_GOI_exp_scatter[,3], strat_GOI_exp_scatter[,1], main="Ordered expression ratio", 
  #      xlab="Expression Ratio", ylab="log2(TPM)", col=alpha("#F8766D", 0.5), pch=15, cex=0.8,
  #      cex.main=3, cex.lab=2.5, cex.axis=2.0, xaxt="none", ylim=c(0,15))
  # points(strat_GOI_exp_scatter[,2], col=alpha("#00BFC4", 0.5), pch=17, cex=0.8)
  # abline(v=max(strat_GOI_exp_scatter[strat_GOI_exp_scatter$percentile_rank <= percentile, "rank"]), col="blue", lty=5, lwd=2)
  # abline(v=min(strat_GOI_exp_scatter[strat_GOI_exp_scatter$percentile_rank >= (100-percentile), "rank"]), col="red", lty=5, lwd=2)
  # legend("top", legend=c(colnames(strat_GOI_exp_scatter[1]), colnames(strat_GOI_exp_scatter[2])), 
  #        col=c("#F8766D", "#00BFC4"), pch=c(15,17), cex=1.5, bg="transparent", bty = "n")
  # invisible(dev.off())
}
if (strat_do_multi_GOI_analysis){
### Do Boxplots with sina
  strat_GOI_exp_sep_long$cols_column <- ifelse(strat_GOI_exp_sep_long$patientID %in% rownames(design[design$lo == 1,]), "#00BFC4", 
    ifelse(strat_GOI_exp_sep_long$patientID %in% rownames(design[design$hi == 1,]), "#F8766D", "grey"))
  strat_GOI_exp_sep_long$size_column <- ifelse(strat_GOI_exp_sep_long$patientID %in% rownames(design[design$lo == 1,]), 2, 
    ifelse(strat_GOI_exp_sep_long$patientID %in% rownames(design[design$hi == 1,]), 2, 0.75))
  strat_GOI_exp_sep_long$in_column <- strat_GOI_exp_sep_long$patientID %in% rownames(design)
  strat_GOI_exp_sep_long <- strat_GOI_exp_sep_long[order(strat_GOI_exp_sep_long$geneID, strat_GOI_exp_sep_long$in_column), ]
  strat_GOI_exp_sep_long <- strat_GOI_exp_sep_long[!is.infinite(strat_GOI_exp_sep_long$expression),]

  psina <- ggplot(strat_GOI_exp_sep_long, aes(geneID, expression))
  psina <- psina + geom_sina(colour = strat_GOI_exp_sep_long$cols_column, size = strat_GOI_exp_sep_long$size_column, alpha = 0.7) + 
    theme_classic() + theme(text=element_text(size=30)) + ylab(paste("Log2(", strat_exp_type, ")", sep=""))
  # invisible(ggsave(paste(GOI_label, "TPM_Boxplot_Sina.svg", sep="_"), last_plot(), width = 10, height = 8))
  
  p_ly = ggplotly(psina) %>% 
    style(text = paste("<b>Patient ID:</b>", strat_GOI_exp_sep_long$patientID, "<br><b>Expression:</b>", strat_GOI_exp_sep_long$expression))
  
  saveWidget(ggplotly(p_ly), file = paste(GOI_label, "TPM_Boxplot_Sina.html", sep="_"))
}
if (strat_do_single_GOI_analysis){
  ### Do Histogram
  phist <- ggplot(log2(strat_GOI_exp_hist), aes(x=get(GOI))) + 
    geom_histogram(color="black", fill="lightgrey", bins=100) + 
    labs(title=GOI,x=paste("Log2(", strat_exp_type, ")", sep=""), y = "Frequency") +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5), text=element_text(size=30))
  phist <- phist + geom_vline(aes(xintercept=mean(log2(strat_GOI_exp_hist[, GOI]))), color="purple", linewidth=0.5) + 
    geom_vline(aes(xintercept=median(log2(strat_GOI_exp_hist[, GOI]))), color="green", linewidth=0.5) +
    geom_vline(aes(xintercept=max(log2(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank <= percentile, GOI]))), color="blue", linewidth=0.5) +
    geom_vline(aes(xintercept=min(log2(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank >= (100-percentile), GOI]))), color="red", linewidth=0.5)
  # invisible(ggsave(paste(GOI_label, "TPM_histogram.svg", sep="_"), width = 7, height = 7))
  
  p_ly <- ggplotly(phist)
  p_ly$x$data[[1]]$text <- gsub("get(GOI)", GOI, p_ly$x$data[[1]]$text, fixed = TRUE)
  
  saveWidget(ggplotly(p_ly), file = paste(GOI_label, "TPM_histogram.html", sep="_"))
  
  # svg(paste(GOI_label, "TPM_histogram.svg", sep="_"))
  # par(mar=c(5,5.5,4,2))
  # hist(log2(strat_GOI_exp_hist[, GOI]), breaks=100, main=GOI, 
  #      xlab=paste("Log2(", strat_exp_type, ")", sep=""), ylab="Frequency",
  #      cex.main=3, cex.lab=2.5, cex.axis=2.0)
  # abline(v=mean(log2(strat_GOI_exp_hist[, GOI])), col="purple")
  # abline(v=median(log2(strat_GOI_exp_hist[, GOI])), col="green")
  # abline(v=max(log2(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank <= percentile, GOI])), col="blue")
  # abline(v=min(log2(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank >= (100-percentile), GOI])), col="red")
  # # abline(v=max(log2(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$quantile_rank == 1, GOI])), col="blue")
  # # abline(v=min(log2(strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$quantile_rank == 4, GOI])), col="red")
  # invisible(dev.off())
}

### Do tumor - normal expression boxplot
strat_GOI_exp_boxplot$Normal <- grepl(pattern = "(-11A|-11B)$", rownames(strat_GOI_exp_boxplot))
if( nrow(strat_GOI_exp_boxplot[strat_GOI_exp_boxplot$Normal == TRUE,]) > 0 ){
  TandN_samples <- TRUE
} else {
  TandN_samples <- FALSE
}
strat_GOI_exp_boxplot$Normal <- as.factor(strat_GOI_exp_boxplot$Normal)
strat_GOI_exp_boxplot[, GOI] <- log2(strat_GOI_exp_boxplot[, GOI])
strat_GOI_exp_boxplot <- rownames_to_column(strat_GOI_exp_boxplot, var = "ID")

suppressWarnings(p <- ggplot(strat_GOI_exp_boxplot, aes(x=Normal, y=get(GOI), color=Normal, label=ID)) + 
  geom_boxplot(notch=TRUE, outliers=FALSE) + theme_classic() + 
  theme(text=element_text(size=30)))
p <- p + geom_jitter(shape=16, position=position_jitter(0.2)) + 
  ylab(paste("Log2(", GOI, " ", strat_exp_type, ")", sep=""))
# invisible(ggsave(paste(GOI_label, "TPM_N-T_boxplot.svg", sep="_"), width = 7, height = 7))

p_ly = suppressWarnings(ggplotly(p))
if( TandN_samples ) {
  p_ly$x$data[[3]]$text <- gsub("get(GOI)", GOI, p_ly$x$data[[3]]$text, fixed = TRUE)
  p_ly$x$data[[4]]$text <- gsub("get(GOI)", GOI, p_ly$x$data[[4]]$text, fixed = TRUE)
  p_ly$x$data[[3]]$text <- gsub("Normal: FALSE<br />", "", p_ly$x$data[[3]]$text, fixed = TRUE)
  p_ly$x$data[[4]]$text <- gsub("Normal: TRUE<br />", "", p_ly$x$data[[4]]$text, fixed = TRUE)
} else {
  p_ly$x$data[[2]]$text <- gsub("get(GOI)", GOI, p_ly$x$data[[2]]$text, fixed = TRUE)
  p_ly$x$data[[2]]$text <- sub("Normal: TRUE<br />", "", p_ly$x$data[[2]]$text, fixed = TRUE)
}

saveWidget(ggplotly(p_ly), file = paste(GOI_label, "TPM_N-T_boxplot.html", sep="_"))

# svg(paste(GOI_label, "TPM_N-T_boxplot.svg", sep="_"))
# p <- ggplot(strat_GOI_exp_boxplot, aes(x=Normal, y=get(GOI), color=Normal)) + 
#   geom_boxplot(notch=TRUE) + theme_classic() + 
#   theme(text=element_text(size=30))
# p + geom_jitter(shape=16, position=position_jitter(0.2)) + 
#   ylab(paste("Log2(", GOI, " ", strat_exp_type, ")", sep=""))
# invisible(dev.off())

if (run_DE) {
  ### Do stratified expression boxplot
  strat_GOI_exp_low <- strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank <= percentile, ]
  strat_GOI_exp_low$stratification <- "Low"
  strat_GOI_exp_high <- strat_GOI_exp_wstrat[strat_GOI_exp_wstrat$percentile_rank >= (100-percentile), ]
  strat_GOI_exp_high$stratification <- "High"
  strat_GOI_selected_strat <- rbind(strat_GOI_exp_low, strat_GOI_exp_high)
  strat_GOI_selected_strat[, GOI] <- log2(strat_GOI_selected_strat[, GOI])
  strat_GOI_selected_strat <- rownames_to_column(strat_GOI_selected_strat, var = "ID")
  
  p <- ggplot(strat_GOI_selected_strat, aes(x=stratification, y=get(GOI), color=stratification, labels=ID)) + 
    geom_boxplot(notch=FALSE, na.rm=TRUE) + 
    ylab(paste("Log2(", GOI, " ", strat_exp_type, ")", sep="")) + 
    theme_classic() + theme(text=element_text(size=30))
  p <- p + geom_jitter(na.rm=TRUE, shape=16, position=position_jitter(0.2)) + 
    scale_x_discrete(limits=c("Low", "High"))
  
  p_ly = ggplotly(p)
  p_ly$x$data[[3]]$text <- gsub("get(GOI)", GOI, p_ly$x$data[[3]]$text, fixed = TRUE)
  p_ly$x$data[[4]]$text <- gsub("get(GOI)", GOI, p_ly$x$data[[4]]$text, fixed = TRUE)
  p_ly$x$data[[3]]$text <- sub("stratification: High<br />", "", p_ly$x$data[[3]]$text, fixed = TRUE)
  p_ly$x$data[[4]]$text <- sub("stratification: Low<br />", "", p_ly$x$data[[4]]$text, fixed = TRUE)
  
  saveWidget(ggplotly(p_ly), file = paste(GOI_label, "TPM_strat_boxplot.html", sep="_"))
  
  # svg(paste(GOI_label, "TPM_strat_boxplot.svg", sep="_"))
  # p <- ggplot(strat_GOI_selected_strat, aes(x=stratification, y=get(GOI), color=stratification)) + 
  #   geom_boxplot(notch=FALSE, na.rm=TRUE) + 
  #   ylab(paste("Log2(", GOI, " ", strat_exp_type, ")", sep="")) + 
  #   theme_classic() + theme(text=element_text(size=30))
  # p + geom_jitter(na.rm=TRUE, shape=16, position=position_jitter(0.2)) + 
  #   scale_x_discrete(limits=c("Low", "High")) 
  # # p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth=0.1) + 
  # invisible(dev.off())
}

### clean up and ready for DE TODO
rm(list=ls(pattern="strat_*"))

### all done here, sign off stratification
message(paste("Finished stratification of", GOI))

###*****************************************************************************
# Setup for DE ####
###*****************************************************************************
if (!run_DE) {
  warning(paste("Unfortunately, there is not enough observations in this dataset for", GOI))
  warning("Try inspecting the distribution plots to identify the cause")
  stop("Ending this process")
}

message(paste("Beginning DE Analysis for", GOI))

### setup group info
hi <- 'hi'
lo <- 'lo'

### Read in data
raw_data <- fread(file=gene_counts_filename, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
message(paste("Raw count expression data loaded for", GOI, ":", ncol(raw_data)-1, "patients and", nrow(raw_data), "genes"))

### If requested, remove GOI from datatable to negate the influence of these genes on the MDS
### Also, removes these points from volcano and heatmap and associated DE tables and others
if(remove_GsOI) {
  raw_data <- raw_data[!raw_data$gene_name == GsOI_split,]
}

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
keep <- rowSums(cpm(y$counts)>low_gene_thrs) >= smallest_replicates_group_size
y <- y[keep, , keep.lib.sizes=FALSE]
# dim(y)
message(paste("Filtering lowly expressed genes using a cpm cutoff of", low_gene_thrs, "in at least", smallest_replicates_group_size,
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
suppressWarnings(plotMeanVar(y, show.tagwise.vars = TRUE, NBline = TRUE))
invisible(dev.off())

png("bcv.png");
suppressWarnings(plotBCV(y))
invisible(dev.off())

### Normalized counts raw
norm_raw <- y$counts
o <- order(rownames(norm_raw))
norm_raw <- norm_raw[o,]
# write.csv(cbind(y$genes[rownames(norm_raw),1],norm_raw), "gene_normalised_expression_data.csv", row.names=FALSE)
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
## Do DE test (exact) ####
###*****************************************************************************
### first check if switch comparison requested 
if (switchDE == "true") {
  comp <- c(hi, lo)
} else {
  comp <- c(lo, hi)
}

de <- exactTest(y, pair=comp)
dt <- decideTests(de,adjust.method="fdr",p.value=dt_pvalue,lfc=dt_lfc)
message(paste("Differential expression using exact test performed", ":", 
              summary(dt)[[1]], "genes down regulated and", 
              summary(dt)[[3]], "genes up regulated"))

### save toptags and create .rnk file for GSEA
tt <- topTags(de, n=nrow(y))
write.csv(tt, "top_tags.csv", row.names=FALSE)
### create ranked file for GSEA (NOTE: not using this but will still produce it so that it can be used manually if required)
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
      display.columns=c("ID","GeneName","DE"), status=dt,anno=anno, 
      cols=col_list, html=paste(GOI_label,"High_Vs_Low-Volcano", sep="-",collapse=""), 
      launch=FALSE))

### Do the volcano plots - svg versions
# svg(paste("High_Vs_Low", GOI_label,"volcano.svg",sep="_", collapse=""))
png(paste("High_Vs_Low", GOI_label, "volcano.png", sep="_", collapse=""))
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
message(paste("Finished DE Analysis for", GOI))

###*****************************************************************************
## Do WebGestalt ####
###*****************************************************************************
### Do WebGestalt
message(paste("Starting WebGestalt Analysis for", GOI))
p_load(WebGestaltR)

main_dir <- getwd()
dir.create(file.path(main_dir, "WebGestalt"))
WG_outdir <- paste(main_dir, "WebGestalt", sep="/")

### setup gene lists
# subset all sig-up and sig-down regulated genes. If more than 500, take only top 500. If less than 20, don't run.
FDR_sig_df_up <- FDR_sig_df[FDR_sig_df$logFC > 0,]
FDR_sig_df_down <- FDR_sig_df[FDR_sig_df$logFC < 0,]
FDR_sig_df_ord_up <- FDR_sig_df_up[order(FDR_sig_df_up$logFC, decreasing = TRUE), ]$gene_name
FDR_sig_df_ord_down <- FDR_sig_df_down[order(FDR_sig_df_down$logFC), ]$gene_name
FDR_sig_df_ord_up_test <- ifelse(length(FDR_sig_df_ord_up) > ORA_min, TRUE, FALSE)
if (length(FDR_sig_df_ord_up) > ORA_max) { FDR_sig_df_ord_up <- FDR_sig_df_ord_up[1:ORA_max] }
FDR_sig_df_ord_down_test <- ifelse(length(FDR_sig_df_ord_down) > ORA_min, TRUE, FALSE)
if (length(FDR_sig_df_ord_down) > ORA_max) { FDR_sig_df_ord_down <- FDR_sig_df_ord_down[1:ORA_max] }

### setup reference list
reference_geneset <- tt$table$gene_name

### setup databases required
enrichDatabaseGO <- c("geneontology_Biological_Process_noRedundant","geneontology_Molecular_Function_noRedundant")
enrichDatabasePW <- c("pathway_KEGG","pathway_Reactome")
enrichDatabaseDS <- c("disease_Disgenet","disease_GLAD4U","disease_OMIM")
enrichDBs <- list(enrichDatabaseGO = enrichDatabaseGO, enrichDatabasePW = enrichDatabasePW, enrichDatabaseDS = enrichDatabaseDS)

### Run each analysis
sink(nullfile()) # don't want console messages
if (FDR_sig_df_ord_up_test){
  message("Passing DE UP regulated geneset to WebGestalt ORA")
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
} else {
  message("Not enough DE up regulated genes to execute WebGestalt")
}
if (FDR_sig_df_ord_down_test){
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
                                minNum=5, maxNum=2000, sigMethod="top", reportNumr=40,
                                isOutput=TRUE, 
                                outputDirectory=WG_outdir, 
                                projectName=paste(names(enrichDBs[i]), "DOWN-Reg", sep = "_"))
  }
} else {
  message("Not enough DE down regulated genes to execute WebGestalt")
}
sink()
### All done for WebGestalt
message(paste("Finished WebGestalt Analysis for", GOI))

###*****************************************************************************
## Do GSEA ####
###*****************************************************************************
### Do GSEA - If Requested
if (gsea_exe != "false"){
  gsea_exe <- paste(gsea_exe, "GSEAPreranked", sep =" ")
  message(paste("Starting GSEA Analysis for", GOI))
  
  main_dir <- getwd()
  dir.create(file.path(main_dir, "GSEA"))
  outdir <- paste(main_dir, "GSEA", sep="/")
  
  if (rna_species == "mRNA"){
    ### Hallmark
    message(paste("Passing DE results to GSEA using MSigDB Hallmark collection:", GOI))
    gsea_gmx <- "-gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v2023.1.Hs.symbols.gmt"
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
      combined_gsea_sig_results$NES <- as.numeric(combined_gsea_sig_results$NES)
      combined_gsea_sig_results <- combined_gsea_sig_results[order(combined_gsea_sig_results$NES),]
      if (GSEA_output_report_file == GSEA_output_report_files[length(GSEA_output_report_files)]){
        ### all reports merged so create the plot
        # combined_gsea_sig_results$bar_color <- ifelse(combined_gsea_sig_results$NES > 0, "red", "blue")
        # pdf(file=paste(outdir, "/",sample_name, ".pdf", sep=""),width=15,height=10)
        # par(las=2) ; par(mar=c(10,50,10,5))
        # barplot(combined_gsea_sig_results$NES, main=sample_name, horiz=TRUE, names.arg=row.names(combined_gsea_sig_results), xlab="Normalise Enrichment Score", cex.main=2, cex.lab=1.5, cex.axis=1.0, col=combined_gsea_sig_results$bar_color)
        # invisible(dev.off())
        if (nrow(combined_gsea_sig_results) > 0){
          combined_gsea_sig_results$cols_column <- ifelse(combined_gsea_sig_results$NES > 0, "#00BFC4", "#F8766D")
          
          gseaplot <- ggplot(data=combined_gsea_sig_results, aes(x=NES, y=GS.br..follow.link.to.MSigDB, fill=cols_column)) + 
            geom_bar(stat="identity") + scale_y_discrete(limits=combined_gsea_sig_results$GS.br..follow.link.to.MSigDB) + 
            theme_minimal() + xlab("Normalised Enrichment Score") +
            theme(legend.position="none", axis.title.y=element_blank(), 
                  axis.title.x=element_text(size=20), axis.text.x=element_text(size=18))
          
          setwd(outdir) # to bypass html files directory creation when saving into a folder
          saveWidget(ggplotly(gseaplot), file = paste(sample_name, ".html", sep=""), selfcontained=TRUE)
          setwd(main_dir)
          
          write.csv(combined_gsea_sig_results[,1:(length(combined_gsea_sig_results)-2)],
                    paste(outdir, "/",sample_name, ".csv", sep=""), row.names=FALSE)
        } else {
          message("No significant GSEA results to plot")
        } 
      }
    }
  }
  
  ### All done for GSEA
  message(paste("Finished GSEA Analysis for", GOI))
} 

