dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(Sys.getenv("R_LIBS_USER"))
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library("DESeq2")
library("org.Hs.eg.db")
library('RColorBrewer')

#reading and reformating data
featurecounts <- read.table('counts.txt', header = T, row.names = 1) #60671 obs
attach(featurecounts)
featurecounts$Chr <- NULL
featurecounts$Start <- NULL
featurecounts$End <- NULL
featurecounts$Strand <- NULL
featurecounts$Length <- NULL
featurecounts <- featurecounts[, c(4, 5, 6, 7, 1, 2, 3, 8)] #reoder the columns so that WT samples are followed by mutated samples


#building DESeqDataSet
condition <- c('WT', 'WT', 'WT', 'WT', 'mut', 'mut', 'mut', 'mut')
colData <- data.frame(row.names=colnames(featurecounts), condition=factor(condition, levels=c('WT','mut')))

dataset <- DESeqDataSetFromMatrix(countData = featurecounts,
                                  colData = colData,
                                  design = ~condition)
dds <- DESeq(dataset)

#performing differential expression analysis
result <- results(dds, contrast=c('condition','WT','mut'))
summary(result)
result <- result[complete.cases(result),]  #remove any rows with NA
resOrdered <- result[order(result$padj),]  #order by adjusted p-value
write.csv(resOrdered, "WT_vs_mutated.csv") #plain result table

#retrieving annotations for a user-friendly result
resulst_annot <- as.data.frame(resOrdered) #19490 obs
annot <- select(org.Hs.eg.db, keys = rownames(resulst_annot), columns = c("SYMBOL", "GENENAME"), keytype = "ENSEMBL") #19618 obs
#'select()' returned 1:many mapping between keys and columns -> need to get rid of duplicates
resulst_annot <- merge(resulst_annot, annot, by.x=0, by.y="ENSEMBL") #19618 obs
resulst_annot <- resulst_annot[!duplicated(resulst_annot$Row.names), ] #19490 obs
write.csv(resulst_annot, "WT_vs_mutated_annotated.csv") #result table with annotations

#PCA plot
rld <- rlogTransformation(dds, blind=TRUE)
plotPCA(rld, returnData=TRUE)

#plot log fold change vs mean expression for all the genes
#genes with p-value < 0.1 colored red
dev.new()
plotMA(result, main='DESeq2: WT vs. mutated', ylim=c(-7,7))


#heatmaps for significant (padj < 0.01) up-regulated and down-regulated genes with FC > 2
up <- resOrdered[ resOrdered[,'log2FoldChange'] > 1, ]
up <- up[ up[,'padj'] < 0.01, ]
#dim(up)  -- 20 genes

down <- resOrdered[ resOrdered[,'log2FoldChange'] < 1, ]
down <- down[ down[,'padj'] < 0.01, ]
#dim(down) -- 47 genes


nCounts <- counts(dds, normalized=TRUE)
hmcol <- rev(brewer.pal(11,'RdYlBu')) #red-high expression, blue - low expression

#substitute EnsebleIDs for gene symbols
a <- as.matrix(nCounts[ row.names(up), ])
names_a <- select(org.Hs.eg.db, keys = rownames(a), columns = "SYMBOL", keytype = "ENSEMBL")
a <- merge(a, names_a, by.x=0, by.y="ENSEMBL")
a <- a[!duplicated(a$Row.names), ]
row.names(a) <- a$Row.names
a$Row.names <- NULL

b <- as.matrix(nCounts[ row.names(down), ])
names_b <- select(org.Hs.eg.db, keys = rownames(b), columns = "SYMBOL", keytype = "ENSEMBL")
b <- merge(b, names_b, by.x=0, by.y="ENSEMBL")
b <- b[!duplicated(b$Row.names), ]
row.names(b) <- b$Row.names
b$Row.names <- NULL

# heatmaps
dev.new()
heatmap(as.matrix(a[,1:8]), col = hmcol, mar = c(8,5), Colv = NA, labRow = a$SYMBOL,
        cexRow = 0.9, cexCol = 1, main = "DEGs with padj < 0.01 and log2fold > 1")
# heatmap down
dev.new()
heatmap(as.matrix(b[,1:8]), col = hmcol, mar = c(8,5), Colv = NA, labRow = b$SYMBOL,
        cexRow = 0.7, cexCol = 1, main = "DEGs with padj < 0.01 and log2fold < 1" )

#volcano plot
alpha <- 0.01 # Threshold on the adjusted p-value
cols <- densCols(result$log2FoldChange, -log10(result$pvalue))
dev.new()
plot(result$log2FoldChange, -log10(result$padj), col=cols, panel.first=grid(),
     main="WT vs mutated", xlab="log2 fold change", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.8)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
#writing the names for the best DEGs
gn.selected <- abs(result$log2FoldChange) > 2 & result$padj < alpha/2
labels <- select(org.Hs.eg.db, keys = rownames(result)[gn.selected ], columns = "SYMBOL", keytype = "ENSEMBL")
text(result$log2FoldChange[gn.selected],
     -log10(result$padj)[gn.selected],
     lab=labels$SYMBOL, cex=0.4)
