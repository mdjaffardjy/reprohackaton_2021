
## Installer le package FactoMineR et factoextra pour l'ACP
install.packages("FactoMineR")
install.packages("factoextra")

## Installer le package EnhancedVolcano pour l'affichage des résultats de l'expression différentielle :
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

# Récupérer le chemin du répertoire de travail :
chemin = getwd()
setwd(chemin)


## Charger l'environnement Rdata sauvegardé lors de l'analyse statistique effectuée depuis Snakemake :
load("fichier.Rdata")
# Si vous avez une erreur d'ouverture du fichier, mettre le chemin absolu du fichier dans load (.../fichier.Rdata)

#########################################
# Expression différentielle des gènes : #
#########################################

# On enlève les NA du tableau
resultats <- resultats[complete.cases(resultats[,6]),]

## Comme on ne veut pas afficher de noms car cela rend le graphique moins lisible
# On veut lui dire de n'afficher les noms que des genes etant tres a gauche ou tres a droite du graphe
keyvals <- ifelse(resultats$log2FoldChange < -10, 'petit',
                  ifelse(resultats$log2FoldChange > 10, 'grand','normal'))

library(EnhancedVolcano)

## Création du plot pour l'expression différentielle :

png(filename = "enhancedVolcano.png")
EnhancedVolcano(resultats,
                lab = rownames(resultats),
                x = 'log2FoldChange', selectLab = rownames(resultats)[which(keyvals %in% c('petit', 'grand'))],
                y = 'padj', 
                title = 'WT VS muté avec p-value ajustée', xlim = c(-10,10), ylim = c(-1,20),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0, drawConnectors = TRUE,widthConnectors = 0.75,
                labSize = 6.0)
dev.off()


########################################
# Analyse en composantes principales : #
########################################


library("FactoMineR")
library("factoextra")

# Mettre les données en forme pour l'ACP :
essai2 = t(essai)
res_pca2 = PCA(X = essai2, graph = FALSE)
# --> regarder le PCA graph of individuals : on a bien 5 groupés et 3 pas groupés ! Et c'est les bons (mutés) qui sont pas groupés!

## Création du plot pour l'expression différentielle :

png(filename = "graph_individuals_pca.png")
fviz_pca_ind (res_pca2,
              repel = TRUE, # Évite le chevauchement de texte, 
              col.ind = mutations, # colorer by groups
              palette = c("#00AFBB", "#E7B800"),
              legend.title = "Groups",
              mean.point = FALSE
)
dev.off()
