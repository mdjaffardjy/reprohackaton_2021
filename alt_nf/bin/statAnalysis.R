
#!/usr/bin/env Rscript 
set.seed(271120)
# This script is using the following librairies, already downloaded in the Docker image :
# DESeq2 : makes RNA-seq data analysis
# FactoMineR : makes PCA
# factoextra : diplays PCA graphics
# pheatmap : displays heat map
# RColorBrewer : contains colors for heat map

library("DESeq2")
library("FactoMineR")
library("factoextra")
library("pheatmap")
library("RColorBrewer")

#------------------------
# Getting data from bash script
#------------------------

args = commandArgs(TRUE) # list with two entries : meta-data and count data
metaData = read.table(args[1], header = TRUE, sep = ",", skip = 2, stringsAsFactors = TRUE)
countData = read.csv(args[2], header = TRUE, sep = "", skip = 1, stringsAsFactors = TRUE)

#------------------------
# Cleaning data
#------------------------

countData = countData[,-c(2,3,4,5,6)] # removing unused columns
                                      # colnames are : gene SRRxxxxxx.bam ... SRRxxxxxx.bam
colnames(countData)[c(2:9)] = gsub('.{0,4}$', '', colnames(countData)[c(2:9)]) # removing .bam (4 last characters)

# removing ENSG and 00000 from Gene ID to shorten names
countData$Geneid = gsub("[a-zA-Z ]", "", countData$Geneid)
countData$Geneid = as.integer(countData$Geneid)
countData$Geneid = as.factor(countData$Geneid)

#------------------------
# Analyzing data
#------------------------

## Changing data format (countData)
# DESeqDataSetFromMatrix
# design : we are intersting in SF3B1 mutation status, the formule is : data ~ sf3b1_status
# countData : counting matrix
# colData : meta data
# tidy : TRUE if countData first column is first row from metaData
DESeqData = DESeqDataSetFromMatrix(countData=countData,
                                   colData=metaData,
                                   design=~sf3b1_status,
                                   tidy = TRUE)

## Calculate statistics
# DESeq : Differential expression analysis based on the Negative Binomial distribution
# DESeq contains three functions :
# 1- Estimation of size factors
# 2- Dispersion
# 3- Wald test and adaptation to negative binomial distribution

analysis = DESeq(DESeqData) # results are in DESeq format

#------------------------
# Preparing results and saving
#------------------------
# Extracting results from DESeq analysis, giving samples names, mean expression, log2 fold change
# standards error, statistic test values, p-value and adjusted p-value

res = results(analysis)
write.table(res, "DESeq2_results.txt", sep = " ", row.names = TRUE, col.names = TRUE)

#------------------------
# Thresholds and selected genes IDs
#------------------------

#Thresholds
padj_thresh = 0.01
FC_thresh = 2

# Keep 20 genes with the highest abs(log2FoldChange)
res_sig <- subset(res, padj<padj_thresh)
res_lfc <- subset(res_sig, abs(log2FoldChange) > FC_thresh)
selected_genes <- order(abs(res_lfc$log2FoldChange), decreasing=TRUE)
selected_genes <- selected_genes[1:min(length(selected_genes),20)]

#------------------------
# Heat Map
#------------------------

# Normalizing results
pheatmap_data = analysis
pheatmap_data_norm = vst(pheatmap_data) #variance stabilizing transformation : for normalization
pheatmap_data_norm_df = as.data.frame(colData(pheatmap_data)[,c("sf3b1_status")])

# Changing rows and columns names
colnames(pheatmap_data_norm_df) = c("sf3b1_status")
pheatmap_data_to_plot = assay(pheatmap_data_norm)[selected_genes,]
rownames(pheatmap_data_norm_df) = colnames(pheatmap_data_to_plot) # mandatory to remove the error :
                        # Error in check.length("fill") : 'gpar' element 'fill' must not be length 0

# Gathering by sf3b1_status value (WT of R625 mutation)
new_order = order(pheatmap_data_norm_df$sf3b1_status)
pheatmap_data_to_plot = pheatmap_data_to_plot[,new_order]

# Figure to save - automatically saved in pdf
heatmap_to_save = pheatmap(pheatmap_data_to_plot, annotation_col=pheatmap_data_norm_df,
                           cluster_rows=T, cluster_cols=T, show_rownames=T,
                           main="Counting data for the 20 most expressed genes with VST method")

# genes are ordered by abs(log2FoldChange)
# values are normalized counting data

#------------------------
# Volcano Plot
#------------------------

output_name = args[3]
title3 = paste(output_name, "_VP.pdf", sep = "")
pdf(title3, width = 10, height = 10)

title = paste("Volcano plot with padj_thresh =", as.character(padj_thresh),
                "and log2FC_tresh =", as.character(FC_thresh), sep = " ")
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=title, xlim=c(-10,10)))
with(subset(res, padj<padj_thresh),
     points(log2FoldChange, -log10(pvalue),
            pch=20, col="blue"))
with(subset(res, padj<padj_thresh & abs(log2FoldChange)>FC_thresh),
     points(log2FoldChange, -log10(pvalue),
            pch=20, col="red"))
dev.off()


#------------------------
# PCA
#------------------------

# # We build a PCA with counting data, but coloring with DESeq analysis results
# PCAdata = countData
# PCAdata = na.omit(PCAdata)
# rownames(PCAdata) = PCAdata[,1]
# PCAdata = PCAdata[,-1]
# PCAdata = t(PCAdata) # transposing in order to make variables with genes and individuals with samples
# PCAdata.acp = PCA(PCAdata, graph=FALSE) # PCA without graph
# 
# selected_genes = results(analysis, tidy=TRUE)
# selected_genes = selected_genes[abs(selected_genes$log2FoldChange) > FC_thresh& selected_genes$padj < padj_thresh,] # keeping red genes with low pvalue and high log2FC
# 
# # Figure to save - automatically saved in pdf
# pca_to_save = fviz_pca_biplot(PCAdata.acp, label=c("var","ind"), col.ind = as.factor(metaData$sf3b1_status), repel = TRUE,
#                               mean.point = FALSE, labelsize = 2, legend.title = "Statut", select.var = list(name = selected_genes$row))
#

# keeping genes with low pvalue and high log2FC
selected_genes = results(analysis, tidy=TRUE)
selected_genes = selected_genes[abs(selected_genes$log2FoldChange) > FC_thresh & 
                                  selected_genes$padj < padj_thresh,]


# PCA data : countData = count(analysis)
PCAdata = vst(counts(analysis)) # normalizing data with vst
PCAdata = t(PCAdata) # transposing in order to make variables with genes and individuals with samples
PCAdata.acp = PCA(PCAdata, graph=FALSE) # PCA without graph

# Figure to save
pca_to_save = fviz_pca_biplot(PCAdata.acp, label=c("var","ind"),
                              col.ind = as.factor(metaData$sf3b1_status), repel = TRUE,
                              addEllipses = TRUE, ellipse.type = "convex",
                              mean.point = FALSE, labelsize = 4, legend.title = "Statut",
                              select.var = list(name = selected_genes$row))

# Associated scree plot
scree_to_save = fviz_eig(PCAdata.acp, ncp = 7)

#------------------------
# Opening PDF to save figures
#------------------------

output_name = args[3]

title1 = paste(output_name, "_PCA.pdf", sep = "")
pdf(title1, width = 10, height = 10)
pca_to_save

title2 = paste(output_name, "_HM.pdf", sep = "")
pdf(title2, width = 10, height = 10)
heatmap_to_save

title4 = paste(output_name, "_SCREE.pdf", sep = "")
pdf(title4, width = 10, height = 10)
scree_to_save

#------------------------------------------
# Number of genes differentially expressed 
#------------------------------------------
