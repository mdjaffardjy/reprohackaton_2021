# Charger les librairies :
library("DESeq2")


########################################################################################################
# La fonction DESeqDataSetFromMatrix prends en entree :                                                #
#                                                                                                      #
#     - une matrice de nombres qui sont les counts                                                     #
#     - le nom des colonnes au format data.frame, pour pouvoir faire les differences selon le nom      #
#     - un shema d'analyse                                                                             #
########################################################################################################

###### Recuperer le chemin (où s'exécute le script et où sont les fichiers dont on a besoin)
dir = getwd()
setwd(dir)

###### Recuperer pour chaque sample les counts :


test1 = read.table(snakemake@input[[1]])
test2 = read.table(snakemake@input[[2]])
test3 = read.table(snakemake@input[[3]])
test4 = read.table(snakemake@input[[4]])
test5 = read.table(snakemake@input[[5]])
test6 = read.table(snakemake@input[[6]])
test7 = read.table(snakemake@input[[7]])
test8 = read.table(snakemake@input[[8]])



##############################################################################################################################################
# les infos sur les counts sont dans la derniere colonne, donc on les reunit dans une seule matrice apres les avoir transformes en numeric : #
##############################################################################################################################################

###### assemblage des derniÃ¨eres lignes de chaque fichier en les tranformant en numeric :
essai = data.frame(as.numeric(test1[-1,7]), as.numeric(test2[-1,7]), as.numeric(test3[-1,7])
                   , as.numeric(test4[-1,7]), as.numeric(test5[-1,7]), as.numeric(test6[-1,7]), 
                   as.numeric(test7[-1,7]), as.numeric(test8[-1,7]))

###### passage au format matrix:
essai = as.matrix(essai)

###### recuperation au format data.frame des labels :
labels = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
mutations = c("M", "M", "WT", "WT", "WT", "WT", "WT", "M")
labels2 = data.frame(labels, mutations)
rowlabels = test1[-1,1]

###### ajout des labels des colonnes et lignes de la matrice :
rownames(essai)=rowlabels
colnames(essai) = labels

###### regarder a quoi cela ressemble :
head(essai)


###### Etude 
dds = DESeqDataSetFromMatrix(countData = essai, colData = labels2, ~mutations)
dds2 = DESeq(dds)
resultats = results(dds2)

head(resultats)

############# Ancien graphique, qu'on a remplacé par l'utilisation de EnhencedVolcano :

#par(mfrow=c(1,1))
#with(resultats, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-8,8)))
#with(subset(resultats, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
#with(subset(resultats, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

summary(resultats)

save(dds, dds2, resultats, essai, mutations, file=snakemake@output[[1]])

