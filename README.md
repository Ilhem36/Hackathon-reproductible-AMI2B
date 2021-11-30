# Projet Hackathon 
# Auteurs : Weng Lorraine, Riahi Ilhem, Ammar Fatma, Cherifi Sana

#Description:
Le but de ce projet est de reproduire des parties de l'analyse décrite dans ces articles : https://pubmed.ncbi.nlm.nih.gov/23313955/ https://pubmed.ncbi.nlm.nih.gov/23861464/

Ils ont effectué un RNA-Seq dans des échantillons de patients atteints de mélanome de l'uvée. Certains échantillons sont mutés en SF3B1 qui est un facteur d'épissage. Pour en savoir plus sur SF3B1(https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000115524;r=2:197388515-197435079). nf-uveal melanoma/rnaseq est un pipeline bioinformatique qui analyse ces données afin de trouver des gènes exprimés de manière différentielle. En d'autres termes, ce pipeline vise à analyser des gènes plus (ou moins) exprimés dans une condition (échantillons mutés SF3B1) par rapport à une autre (échantillons non mutés SF3B1).

Le pipeline est construit à l'aide de Nextflow, un outil de workflow pour exécuter des tâches sur plusieurs infrastructures de calcul de manière très portable. Il utilise des conteneurs Docker/Singularity rendant l'installation triviale et les résultats hautement reproductibles.

# Etat du projet= finalisé

# Outils utilisés
Les outils ont tous la même provenance, c'est-à-dire qu'ils sont tous issus du même profil hormis le dernier outil. Il s'agit du profil [evolbioinfo] sur Docker. Le dernier outil en question  a été spécifiquement créé dans le cadre de ce projet et est  disponible sur le profil [fatmaammar] toujours sur le Docker.

[samtools](https://github.com/samtools/samtools) : Version 1.11
[sratoolkit](https://hpc.nih.gov/apps/sratoolkit.html) : Version 2.10.8
[STAR](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf) : Version 2.7.6
[Subread](https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf) : Version 2.0.1
- Packages de R avec r-base:3.6.3 : [BiocManager] pour installer quelques packages [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) (version 1.28.1), [FactoMineR](http://factominer.free.fr/index.html) (version 2.3) matrixStats (version 0.61.0) et pour la visualisation : ggplot2 (version 3.3.5), factoextra (version 1.0.7), pheatmap (version 1.0.12), RColorBrewer (version 1.1-2), EnhancedVolcano(version 1.12.0)
Le package tidyverse (version 1.3.1), nous n'avons pas pu l'installer dans le container donc nous avons ajouté une commande dans notre code R qui va nous le faire: install.packages("tidyverse")

# Etapes du workflow

1) "DownloadSRR", "fastqZipping" et "chromosomedownloading": Il s'agit du téléchargement de données des chromosomes humains et mitochondriaux ainsi que la compression du fichier. Il s'agit également de l'annotation du génome humain et des 8 données RNA-seq des 8 individus, issus du projet [SRP017413].
2) "CompleteGenome" et "genomeIndexation" : C'est la construction ainsi que l'indexation du génome humain avec STAR.
3) "HumanGenomeAnnotation": Conversion des données RNA-seq en .fastq avec sratoolkit.
4) "mappingfq": C'est la réalisation de l'alignement et du tri des données RNA-seq sur le génome humain avec STAR. Les sorties sont des fichiers .bam triés.
5) "generating_bam_files" : Il s'agit de l'indexation des fichiers .bam. en .bai avec samtools
5) "counting_Reads_Matrix" : On réalise le comptage des séquences exprimées pour chaque patient avec la fonction featureCount de subread
6) "stat_analysis": C'est l'analyse statistique des résultats avec DESeq2 et FactoMineR.

#Exécution rapide
Afin de pouvoir exécuter le workflow (main.nf), il est nécessaire d'utiliser Nextflow (version 20.10.0). De plus, cette exécution se déroule bien dans une machine virtuelle avec Go de RAM. Pour exécuter le workflow, il faut donc récupérer nextflow.config, main.nf, metaData.txt et stat_dernier.R en local, ouvrir un terminal puis exécuter les commandes :

$ cd ../path/to/the/file
$ conda activate nextflow
$ nextflow run main.nf 
En cas d'interruption d'un process:
$ nextflow run main.nf -resume


#Les différents fichiers utilisés lors de ce projets

- "main.nf" est le workflow final en cours de construction ;
- "Dockerfile" contient le fichier Docker dont l'utilisation nous a permis de créer l'image Docker afin de réaliser l'analyse statistque des résultats ;
- "stat" un dossier contenant le script R qui permet l'excecution du process Stat_analysis
- "nextflow.config" contient les informations de configuration du workflow : utilisation de Docker et chargements des conteneurs ;
- "metaData.txt" contient les méta-informations concernant les 8 échantillons dont nous analysons les données RNA-seq ;
- "stat_dernier.R" il s'agit du script R que l'on retrouve dans le dernier processus du workflow. Il permet de réaliser l'analyse statistique des résultats de comptage qui sont issus des données RNA-seq ;
- "resultats": ce ficher contient tous les résultats obtenus.

#Informations concernant l'installation de l'image du Docker et sa publication sur le hub Docker

Installation et activation de docker image
Pour créer une image Docker, tout d'abord nous avons installé le Docker Desktop et nous avons créé un compte. Ensuite, nous avons créé un Dockerfile contenant nos packages à installer, puis nous avons utilisé la commande: 
docker build -t fatmaammar/test14
Pour utiliser l'image construite localement, nous l'avons partagée sur Docker Hub en faisant:
Push to hub

