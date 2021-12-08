#!/usr/bin/env Rscript 


library("DESeq2")
library("FactoMineR")
library("xlsx")
library("factoextra")
library("pheatmap")
library("RColorBrewer")


args = commandArgs(TRUE)
countData = read.csv(args[1], header = TRUE, sep = "", skip = 1, stringsAsFactors = TRUE)
countData = countData[,-(1:6)] # Suppression de colonnes superflues
gene_id = countData[, 1] # Noms des genes pour renommer les colonnes
row.names(countData) = gene_id # Les lignes sont identifiees par leur gene_id 
colnames(countData)[c(2:9)] = gsub('.{0,4}$', '', colnames(countData)[c(2:9)]) # removing .bam (4 last characters)

conds <- factor(c("mutant", "wt", "wt", "wt", "wt", "wt", "mutant", "mutant")) # Phenotype des echantillons
pheno_sample = c("mtn_89", "wt_88", "wt_87", "wt_86", "wt_85", "wt_84", "mtn_83", "mtn_82")

coldata <- data.frame(condition=conds) # etiqueter chaque echantillon par son type (wt ou mut)
row.names(coldata) = pheno_sample

ddsObjet <- DESeqDataSetFromMatrix(countData=count, colData=coldata, formula(~phen), tidy= TRUE) # On compare selon le phenotype

ddsObjet$condition <- relevel(ddsObjet$condition, ref = "wt")
ddsObjet <- estimateSizeFactors(ddsObjet) 
sizeFactors(ddsObjet) # harbour et al 2013

filtered_row <- rowSums(countData(ddsObjet)) >= 1 # On filtre les lignes sans read
dds <- ddsObjet[filtered_row,]

ddsEstim = DESeq(dds)
res <- results(ddsEstim, alpha=0.05, altHypothesis="greaterAbs") # padj < 0.05
write.table(res, "DESeq2_results.txt", sep = " ", row.names = TRUE, col.names = TRUE)

#Thresholds
padj_thresh = 0.01
FC_thresh = 2

# filter 20 genes with the highest FC
res_sig <- subset(res, padj<padj_thresh)
res_lfc <- subset(res_sig, abs(log2FoldChange) > FC_thresh)
selected_genes <- order(abs(res_lfc$log2FoldChange), decreasing=TRUE)
selected_genes <- selected_genes[1:min(length(selected_genes),20)]


# Volcano plot

pdf("analysis_VP.pdf", width = 10, height = 10)

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

v_plot_cols <- densCols(res_lfc0$log2FoldChange, -log10(res_lfc0$pvalue))
v_plot <- plot(res_lfc0$log2FoldChange, -log10(res_lfc0$padj), col=cols, panel.first=grid(), 
     main="Volcano plot", xlab="log2(fold-change)", ylab="-log10(p-value ajustee)",
     pch=20, cex=0.6, xlim=c(-4,7), ylim=c(-0.05,4))

abline(v=c(-1,1), col="brown") # Seuils
abline(h=-log10(0.05), col="brown")


# ACP

df_sample = as.data.frame(t(countData)) # Echantillons ligne, genes col
row.names(df_sample) = pheno_sample 
df_sample = cbind(df_sample, conds) # label (mutant, wt) 
acp = PCA(df_sample, scale.unit=TRUE, ncp=5, graph = FALSE, quali.sup=c(60613)) 
acp_plot = plot.PCA(acp, invisible=c('quali', 'ind.sup'), habillage=60613, title="ACP des échantillons", col.quali='#C724C7', label=c('ind'))


pdf("analyse_PCA.pdf", width = 10, height = 10)
acp_plot


# Comparaison des genes differentiellement exprimes avec Furney et al

rnaseq = read.csv("mutant_WT_lfc0_results.csv", header=TRUE)
furney = read.xlsx("Furney_genes_diff_exprimes.xlsx", 1, header=TRUE)

vec_colnames = furney[5,]
furney = furney[-c(1:5),] 
colnames(furney) = vec_colnames
furney = furney[, -c(13, 14)] 

vec_rnaseq = rnaseq[,1] 
vec_furney = furney[,12]

table(is.element(vec_rnaseq, vec_furney)) # genes communs entre les 2 analyses
rnaseq[which(is.element(vec_rnaseq, vec_furney)),] # genes detectes comme diff exprimes ds les 2 analyses

vec_harbour = c("ENSG00000132824", "ENSG00000104872", "ENSG00000131503", "ENSG00000166136", "ENSG00000134759", "ENSG00000151882", "ENSG00000168385", "ENSG00000146085", "ENSG00000078618", "ENSG00000116641")
table(is.element(vec_rnaseq, vec_harbour)) # genes communs avec Harbour et al


# Heat Map

heatmap_data = res
heatmap_data = vst(pheatmap_data) #normalisation
heatmap_data_df = as.data.frame(colData(heatmap_data)[,c("phenotype")])

colnames(heatmap_data_df) = c("phenotype")
heatmap_data_plot = assay(heatmap_data)[selected_genes,]
rownames(heatmap_data_df) = colnames(heatmap_data_plot) 

# classement par phenotype
new_order = order(heatmap_data_df$phenotype)
heatmap_data_to_plot = heatmap_data_plot[,new_order]


heatmap = pheatmap(heatmap_data_plot, annotation_col=pheatmap_data_df,
                           cluster_rows=T, cluster_cols=T, show_rownames=T,
                           main="comptes pour les 20 genes les plus exprimés")

pdf("analyse_HM.pdf", width = 10, height = 10)
heatmap





















