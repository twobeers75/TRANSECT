#!/usr/bin/env Rscript

###*****************************************************************************
### Create correlation table and plots from RECODE3 data
### JToubia - January 2023
###*****************************************************************************

###*****************************************************************************
### Import libraries ####
###*****************************************************************************
suppressMessages(if (!require("pacman")) install.packages("pacman"))
p_load(data.table, rlogging)
#suppressMessages(library(data.table))

SetLogFile(base.file=NULL)

###*****************************************************************************
### Read in Args ####
###*****************************************************************************
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (GOI).n", call.=FALSE)
} else if (length(args)==1) {
  # default to
  args[2] = list.files("/home/jtoubia/Desktop/Projects/SRt/RECOUNT3/BRCA/", "R3-BRCA-*_tpm-mRNA.tsv", full.names=TRUE)
  args[3] = "/home/jtoubia/Desktop/Projects/SRt/REF_FILES"
}

GOI <- args[1]
r3_z2n <- args[2]
ref_files_folder <- args[3]

###*****************************************************************************
### Setup for analyses ####
###*****************************************************************************
main_dir <- getwd()
gencode_lookup <- "gencode.v26.annotation.lookup"

###*****************************************************************************
### Start ####
###*****************************************************************************
message(paste("Beginning Corr Analysis for", GOI))

###*****************************************************************************
### Read in RECODE3 data ####
###*****************************************************************************
# df_r3_data <- read.delim(r3_z2n, row.names = 1)
exp_data <- fread(file=r3_z2n, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, data.table=FALSE)
exp_data <- exp_data[order(exp_data[,'gene_name']),]
exp_data <- exp_data[!duplicated(exp_data$gene_name),]
rownames(exp_data) <- exp_data$gene_name
exp_data <- subset(exp_data, select=-c(gene_name))
df_r3_data <- t(exp_data)

###*****************************************************************************
### Check there is enough data for the GOI ####
###*****************************************************************************
gene1 <- GOI
gene1 <- toString(GOI)

if (gene1 %in% colnames(df_r3_data)) {
  gene1_exp_value <- as.numeric(unlist(df_r3_data[,gene1]))
  if (sum(is.finite(gene1_exp_value)) < 50 || sum(gene1_exp_value > 1, na.rm=T) < 10) {
    message(paste("Unfortunately, there is not enough observations in this dataset for", GOI))
    stop("Ending this process") 
  }
}

dir.create(file.path(main_dir, "Corr_Analysis"))
dir.create(file.path(main_dir, "Corr_Analysis", "plots"))
setwd(file.path(main_dir, "Corr_Analysis"))

###*****************************************************************************
### Create Pairs ####
###*****************************************************************************
df_target_pairs <- read.delim(paste(ref_files_folder, gencode_lookup, 
                    sep="/"), header=FALSE, col.names=c("gene1_id", "gene2_id"))
df_target_pairs$gene1_id <- gsub('-', '.', GOI)
df_target_pairs$gene2_id <- gsub('-', '.', df_target_pairs$gene2_id)
df_target_pairs[,c("Cor","Pvalue","logExp_Cor","logExp_Pvalue")] <- ""

###*****************************************************************************
### Corr Analysis ####
###*****************************************************************************
for (row in 1:nrow(df_target_pairs)) {
# for (row in 1:10) { 
  gene2 <- df_target_pairs[row, "gene2_id"]
  gene2 <- toString(gene2)
  
  if (gene2 %in% colnames(df_r3_data)) {
    gene2_exp_value <- as.numeric(unlist(df_r3_data[,gene2]))
    if (sum(is.finite(gene2_exp_value)) < 50 || sum(gene2_exp_value > 1, na.rm=T) < 10) {
      df_target_pairs[row, 3:6] <- "low data"
    } else {
      cor_stats <- cor.test(gene1_exp_value,gene2_exp_value) #,use="pairwise.complete.obs")
      df_target_pairs[row, "Cor"] <- signif(as.numeric(cor_stats$estimate),3)
      df_target_pairs[row, "Pvalue"] <- signif(cor_stats$p.value,3)
      log_cor_stats <- cor.test(log2(gene1_exp_value),log2(gene2_exp_value)) #,use="pairwise.complete.obs")
      df_target_pairs[row, "logExp_Cor"] <- signif(as.numeric(log_cor_stats$estimate),3)
      df_target_pairs[row, "logExp_Pvalue"] <- signif(log_cor_stats$p.value,3)
      
      if(!is.na(log_cor_stats$estimate)){
        if (abs(signif(as.numeric(log_cor_stats$estimate))) > 0.8) {
          png(paste("plots/", gene1, "_", gene2, ".png", sep=""))
          cor_1 <- log_cor_stats
          lm1 <- lm(log2(gene1_exp_value)~log2(gene2_exp_value))
          par(mar=c(6,5.5,4,2))
          plot(log2(gene2_exp_value), log2(gene1_exp_value), pch=20, xlab=paste("log2(", gene2, ")", sep=""),
               ylab=paste("log2(", gene1, ")", sep=""), main="Expression Scatterplot", 
               cex.main=3, cex.lab=2.5, cex.axis=2.0)
          abline(lm1, col="red")
          legend("bottomleft", legend=c(paste("R =", signif(as.numeric(cor_1$estimate), 3)), 
                                        paste("p =", signif(cor_1$p.value, 3))))
          invisible(dev.off())
        }
      }
    }
  } else {
    df_target_pairs[row, 3:6] <- "no data"
  }
}

###*****************************************************************************
### Filter/Sort/Padjust/Write to file and plot then return to main working dir
###*****************************************************************************
### remove unwanted columns (NOTE. not using raw expression for correlations, only logged expression)
cols2keep <- c("gene1_id", "gene2_id", "logExp_Cor", "logExp_Pvalue")
df_target_pairs <- df_target_pairs[cols2keep]
### get rid of NAs, "no data" & "low data" rows (only here for historical reasons)
df_target_pairs <- na.omit(df_target_pairs) 
df_target_pairs <- df_target_pairs[!(df_target_pairs$logExp_Cor=="low data" | df_target_pairs$logExp_Cor=="no data"),]
### calculate adjusted P
df_target_pairs$logExp_Bonferroni <- signif(p.adjust(df_target_pairs$logExp_Pvalue, "bonferroni"),3)
df_target_pairs$logExp_FDR <- signif(p.adjust(df_target_pairs$logExp_Pvalue, "BH"),3)
### sort
df_target_pairs <- df_target_pairs[order(df_target_pairs$logExp_FDR, na.last = NA),]
### write out
write.table(df_target_pairs, paste(GOI, "corr.tsv", sep="_"), sep='\t', row.names=FALSE)

### Plot
df_target_pairs$logExp_FDR[df_target_pairs$logExp_FDR == 0] = 1E-320
df_target_pairs <- transform(df_target_pairs, logExp_Cor = as.numeric(logExp_Cor))

png(paste(GOI, "corr_volcano.png", sep="_"))
par(mar=c(5,6,4,2))
with(df_target_pairs, plot(logExp_Cor, -log10(logExp_FDR), pch=20, cex=0.25, col="grey", main="Pearson's Correlations", 
              xlab="R", ylab="-log10(FDR)", xlim=c(-1,1), ylim=c(0,310), cex.main=2.5, cex.lab=2.5, cex.axis=2.0))

### Add colored points:
# with(subset(df_target_pairs, logExp_Cor > 0.5 ), points(logExp_Cor, -log10(logExp_FDR), pch=20, cex=0.25, col="orange"))
# with(subset(df_target_pairs, logExp_Cor < -0.5 ), points(logExp_Cor, -log10(logExp_FDR), pch=20, cex=0.25, col="orange"))
with(subset(df_target_pairs, logExp_FDR < 1E-50), points(logExp_Cor, -log10(logExp_FDR), pch=20, cex=0.25, col="green"))
with(subset(df_target_pairs, logExp_FDR < 1E-150 & logExp_Cor > 0.8), points(logExp_Cor, -log10(logExp_FDR), pch=20, cex=0.75, col="red"))
with(subset(df_target_pairs, logExp_FDR < 1E-150 & logExp_Cor < -0.8), points(logExp_Cor, -log10(logExp_FDR), pch=20, cex=0.75, col="blue"))
invisible(dev.off())

setwd(main_dir)

### all done, sign off
message(paste("Finished Corr Analysis for", GOI))

