# Projet Repro-Hackathon
### Auteur : Marine Djaffardjy

Ce projet cherche à reproduire en partie les analyses RNA-seq de deux papiers :<br>
  - https://pubmed.ncbi.nlm.nih.gov/23313955/<br>
  - https://pubmed.ncbi.nlm.nih.gov/23861464/<br>

Les analyses effectuées par ces papiers, sur lesquelles nous nous penchons dans ce projet, sont des échantillons RNA-seq obtenus dans le cadre de l'étude de mélanome uvéal, notamment l'expression de SF3B1, dont les mutations sont étudiées. Les données sont disponibles à cette adresse :  http://www.ncbi.nlm.nih.gov/sra?term=SRA062359

Ce projet contient deux workflows, un en snakemake et un en nextflow, qui effectuent la même analyse.
Ces workflows ont pour but d'analyser spécifiquement les données citées ci-dessus (et va ainsi comporter une étape intégrée de téléchargement de ces données). Ils contiennent également une étape d'analyse des données traitées, ainsi les différentes constantes (valeurs de seuil, par exemple), sont fixées pour ces données.

## Pré-requis 

Pour faire tourner les workflows, il faut 8 cpus et 40GB de RAM.
Les logiciels à avoir sont également : 
	nextflow (version 21.10.4)
	docker (version )
	snakemake (version )
	
## Faire tourner les workflows
Pour faire tourner les workflows, notamment sur le cloud IFB, il faut lancer, une fois dans le dossier du dépôt :
```bash
conda activate
nextflow run src/main.nf
```
Pour générer un rapport du workflow et une flowchart qui témoigne du data flow du workflow :

```bash
nextflow run src/main.nf -with-report reports/report.txt -with-dag results/flowchart.png
```

Pour faire tourner le workflow snakemake, il faut lancer :

## Résultats

Les résultats sont visibles dans la partie
