##### README écrit en se basant sur les instructions [IONOS](https://www.ionos.fr/digitalguide/sites-internet/developpement-web/fichier-readme/)
# Projet Repro-Hackathon - Goupe 5 (2020)
### Auteurs : Thomas El Khilali, Audrey Onfroy, Leila Outemzabet

## Etat du projet 
Le projet est fini !

## Description rapide
Le projet Repro-Hackathon s'inscrit dans le master M2 AMI2B. Dans le contexte de la bio-informatique, la reproductibilité est un élément essentiel pour permettre à qui le veuille de reproduire des études. La maîtrise des méthodes appliquées sur des données d'entrée doit toujours pemettre de reproduire les mêmes résultats en sortie.

Notre objectif est de reproduire la méthode appliquée dans les deux articles [Furney *et al.*, Cancer Discovery (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321577/) et [Harbour *et al.*, Nature Genetics (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789378/). Ils disposent des mêmes données de RNA-seq de tissus d'individus atteints ou non du mélanome uvéal, une tumeur de l'œil, les individus atteints du mélanome uvéal étant sujets à un mutation du gène SF3B1. Mais, à l'issue de l'analyse des données RNA-seq, Furney *et al.* observe une expression différentielle de certains gènes entre les sujets porteurs ou non de la mutation tandis que Harbour *et al.* n'ont pas cette conclusion. Afin d'explorer les raisons pour lesquelles les conclusions divergent entre les deux articles, nous avons construit un nouveau workflow, en Nextflow. Nous avons montré que 79 sont différentiellement exprimés telle que le log2 de *Fold Change* est supérieur à 2 en valeur absolue, et que la p-value ajustée (Benjamini-Hochberg) du test de Wald associé rend l'expression différentielle significative au seuil 1%.

## Organisation des fichiers
- **R_container_script** contient le fichier Docker qui nous a permis de créer l'image Docker pour l'analyse stastique des résultats ;
- **bin** contient les exécutables pour le workflow final, à savoir le script R permettant l'analyse statistique des résultats;
- **nextflow.config** contient les informations de configuration du workflow : utilisation de Docker et chargements des conteneurs ;
- **main.nf** est le workflow final en cours de construction ;
- **meta_Deta.txt** contient les méta-informations concernant les 8 échantillons dont nous analysons les données RNA-seq ;
- **statAnalysis.R** est le script R inclus dans le dernier processus du workflow. Il réalise l'analyse statistique des résultats de comptage issus des données RNA-seq ;
- **workflow_info** contient les informations sur l'exécution des différents processus, comme le temps requis ou la mémoire, ainsi que les deux tableaux finaux de données, à savoir la matrice de comptage (Feature Counts) et la matrice des fold change et p-value (DESeq2).

## Outils utilisés
Les outils sont issus des conteneurs du même nom. Hormis le dernier, ils proviennent tous du profil [evolbioinfo](https://hub.docker.com/u/evolbioinfo) sur Docker. Le dernier a été créé spécifiquement pour ce projet et est disponible sur le profil [leilaoutem](https://hub.docker.com/u/leilaoutem), toujours sur Docker.
* [samtools](https://github.com/samtools/samtools) : Version 1.11
* [sratoolkit](https://hpc.nih.gov/apps/sratoolkit.html) : Version 2.10.8
* [STAR](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) : Version 2.7.6a
* [subread](https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf) : Version 2.0.1
* Packages de R : [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) (version 1.28.1), [FactoMineR](http://factominer.free.fr/index.html) (version 2.3) et pour la visualisation : factoextra (version 1.0.7), pheatmap (version 1.0.12) et RColorBrewer (version 1.1-2)

## Exécution rapide
L'exécution du workflow (main.nf) nécessite Nextflow (version *20.10.0*), et se déroule bien dans une machine virtuelle avec **920** Go d'espace disponible et **64** Go de RAM. Pour l'exécution, il faut récupérer *nextflow.config*, *main.nf*, *metaData.txt* et *statAnalysis.R* en local. Ouvrir un terminal puis exécuter les commandes :

```
$ cd ../path/to/the/file
$ conda activate
$ nextflow main.nf -c nextflow.config -with-trace -with-timeline
```

Une fois l'exécution terminée, il s'agit de récupérer le PDF final si le workflow n'a pas été exécuté en local. Le fichier se trouve dans le répertoire results/statResults. Pour le récupérer localement, il suffit d'exécuter la commande ci-dessous dans un terminal ouvert en local, avec $XXX.XXX.XXX.XXX$ désignant l'adresse IP de la machine virtuelle : 

```
$ scp ubuntu@XXX.XXX.XXX.XXX:/path/to/PROJET_REPROCKATHON/results/statResults/output.pdf .
```


## Étapes du workflow
1) **Téléchargement des données** : chromosomes humains et mitochondrial, annotation du génome humain et les 8 données RNA-seq des 8 individus, issus du projet [SRP017413](https://www.ncbi.nlm.nih.gov/sra/?term=SRP017413)
2) **Construction** et **indexation** du génome humain avec STAR
3) **Conversion** des données RNA-seq en *.fastq* avec sratoolkit
4) **Alignement** et **Tri** des données RNA-seq sur le génome humain avec STAR. Les sorties sont des fichiers *.bam* triés.
5) **Indexation** des fichiers *.bam*. en *.bai* avec samtools
5) **Comptage** des séquences exprimées pour chaque patient avec la fonction featureCount de subread
6) **Analyse statistique** des résultats avec DESeq2 et FactoMineR

Les sorties du workflow sont des images vectorielles (format *.pdf*) présentant une Heat Map, un Volcano Plot et deux graphiques d'Analyse en Composantes Principales (ACP), le Scree Plot et le Biplot. Si le workflow a été exécuté sur machine virtuelle, il est possible de récupérer les graphiques avec la commande suivante :
```
$ scp ubuntu@XXX.XXX.XXX.XXX:/path/to/PROJET_REPROHACKATHON/results/statResults/* ./path/to/save/
```

#### Informations pratiques concernant l'installation de Docker et la publication de l'image sur le *hub* Docker
D'après : https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-20-04-fr
```
$ # Installation et activation de docker
$ sudo apt update # update existing packages list
$ sudo apt install apt-transport-https ca-certificates curl software-properties-common
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add - # add GPG key
$ sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu focal stable"
$ sudo apt update # update new packages
$ sudo apt install docker-ce # install docker
$ sudo systemctl status docker # verify everything is activated : Ctrl+C to close
$ sudo service docker start
$
$ cd /path/to/PROJET_REPROCKATHON/
$
$ # Pour créer une image docker, il faut les permissions d'accéder à /var/run/docker.sock
$ ll /var/run/docker.sock
$ sudo chmod 777 -R /var/run/docker.sock
$ docker build ./R_container_script/ -t rtools_stat_image
$
$ # Pour utiliser les images construites localement, il faut les partager sur Docker Hub :
$ docker login
$ docker tag rtools_stat_image leilaoutem/rtools_stat:v2020.11.20
$ docker push leilaoutem/rtools_stat:v2020.11.20 # docker push identifiant/repository:image
```
