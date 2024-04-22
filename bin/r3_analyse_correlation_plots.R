#!/usr/bin/env Rscript

.libPaths( c( "/data/RLibs" , .libPaths() ) )

###*****************************************************************************
### Create correlation table and plots from RECODE3 data
### JToubia - January 2023
###*****************************************************************************

###*****************************************************************************
### Import libraries ####
###*****************************************************************************
suppressMessages(if (!require("pacman")) install.packages("pacman"))
p_load(data.table, plotly, htmlwidgets)
#suppressMessages(library(data.table))

###*****************************************************************************
### Read in Args ####
###*****************************************************************************
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (GOI).n", call.=FALSE)
} else if (length(args)==1) {
  # default to
  args[3] = list.files("/home/jtoubia/Desktop/Projects/SRt/RECOUNT3/BRCA/", "R3-BRCA-*_tpm-mRNA.tsv", full.names=TRUE)
  args[4] = "/home/jtoubia/Desktop/Projects/SRt/REF_FILES"
}

GOI <- args[1]
outdir <- args[2]
r3_z2n <- args[3]
ref_files_folder <- args[4]

###*****************************************************************************
### Setup for analyses ####
###*****************************************************************************
#main_dir <- getwd()
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
    message("Ending this process")
    quit("no", 2, FALSE)
  }
}

dir.create(file.path(outdir, "Corr_Analysis"))
dir.create(file.path(outdir, "Corr_Analysis", "plots"))
setwd(file.path(outdir, "Corr_Analysis"))

###*****************************************************************************
### Create Pairs ####
###*****************************************************************************
df_target_pairs <- read.delim(paste(ref_files_folder, gencode_lookup, 
                    sep="/"), header=FALSE, col.names=c("gene1_id", "gene2_id"))
df_target_pairs$gene1_id <- gsub('-', '.', GOI)
df_target_pairs$gene2_id <- gsub('-', '.', df_target_pairs$gene2_id)
df_target_pairs[,c("Cor","Pvalue","logExp_Cor","logExp_Pvalue")] <- ""

high_expr_df = data.frame(Names=df_r3_data[,1])

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
          # adding exprs to new df, to be outputted in the outputs
          high_expr_df[gene2] <- df_r3_data[gene2]

          #p <- ggplot(aes(x = log2(gene2_exp_value), y = log2(gene1_exp_value))) +
          #  geom_point(shape=20) + 
          #  labs(x = paste("log2(", gene2, ")", sep=""), y = paste("log2(", gene1, ")", sep=""), title = "Expression Scatterplot") + 
          #  theme_classic() +
          #  theme(text=element_text(size=30), 
          #        legend.position = "bottomleft", 
          #        legend.text=c(paste("R =", signif(as.numeric(cor_1$estimate), 3)), 
          #                      paste("p =", signif(cor_1$p.value, 3)))) +
          #  geom_smooth(formula=log2(gene1_exp_value)~log2(gene2_exp_value), method=lm, se=FALSE, color = "red")
          #p_ly = ggplotly(p) %>% style(text = paste("<b>Patient ID:</b>", high_expr_df$Names, 
          #                                  "<br><b>log2(", gene1, "):</b>", log2(gene1_exp_value),
          #                                  "<br><b>log2(", gene2, "):</b>", log2(gene2_exp_value),
          #                                  ))
          #saveWidget(ggplotly(p_ly), file = paste("plots/", gene1, "_", gene2, ".html", sep=""))
        }
      }
    }
  } else {
    df_target_pairs[row, 3:6] <- "no data"
  }
}

write.table(high_expr_df, paste(GOI, "most_correlated_gene_exprs.tsv", sep="_"), sep='\t', row.names=FALSE)

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

p <- ggplot(df_target_pairs, aes(x = logExp_Cor, y = -log10(logExp_FDR))) +
  geom_point(shape = 20, size = 0.25, color = "grey") + # adding additoinal points layer with the subset of data we want
  geom_point(data = subset(df_target_pairs, logExp_FDR < 1E-50), shape=20, size=0.25, color="green") +
  geom_point(data = subset(df_target_pairs, logExp_FDR < 1E-150 & logExp_Cor > 0.7), shape=20, size=0.75, color="red") +
  geom_point(data = subset(df_target_pairs, logExp_FDR < 1E-150 & logExp_Cor < -0.7), shape=20, size=0.75, color="blue") +
  labs(title = "Pearson's Correlations", x = "R", y = "-log10(FDR)") +
  theme_classic() +
  theme(text=element_text(size=30)) +
  xlim(-1, 1) +
  ylim(0, 310)

p_ly = ggplotly(p) %>% style(text = paste("<b>Gene 2:</b>", df_target_pairs$gene2_id, 
                                          "<br><b>logExp_Cor:</b>", df_target_pairs$logExp_Cor,
                                          "<br><b>logExp_Pvalue:</b>", df_target_pairs$logExp_Pvalue,
                                          "<br><b>logExp_Bonferroni:</b>", df_target_pairs$logExp_Bonferroni,
                                          "<br><b>logExp_FDR:</b>", df_target_pairs$logExp_FDR
                                          ))
saveWidget(ggplotly(p_ly), file = paste(GOI, "corr_volcano.html", sep="_"))

setwd(outdir)

### all done, sign off
message(paste("Finished Corr Analysis for", GOI))

