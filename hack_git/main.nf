#!/usr/bin/env nextflow
// ====================================================== Main workflow for RNA-seq data analysis ======================================================= // 


File_SRP = Channel.from('SRR628589', 'SRR628588', 'SRR628587', 'SRR628586', 'SRR628585', 'SRR628584', 'SRR628583', 'SRR628582') 
/* Ce process sert à télécharger les 8 fichiers  SRR du site NCBI SRA en utilisant le NCBI_API_KEY mentionné dans le fichier config.nextflow */ 
process DownloadSRR {

/*
Input : 8  identifiants d'échantillons d'ARNm de patients atteints de mélanome
Output: 8 fichiers SRRxxxxxx.sra stockés dans le dossier SRR
Function : Télécharger les fichiers SRR
*/ 
    publishDir "results/SRR/"

    input:
    val srr_id from File_SRP

    output:
    tuple val(srr_id), file("${srr_id}.sra") into file_srr 
    
    script:
    """
    wget -O ${srr_id}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${srr_id}/${srr_id}.1
    """
}

/* Ce process sert à convertir   les 8 fichiers  téléchargés dans le process précedent en format .gz 
NB: chaque échantillon SRR fonctionne en paire de reads ce qui explique 1.fastq.gz et 2.fastq.gz dans le output de ce process
*/ 
process fastqZipping {
/* 	Input: les valeurs des identifiants SRR récupéres depuis file_srr  du process DownloadSRR
	Output: 8 échantillons au format fastq ( 2 paires pour chaque échantillon)
	Fonction : compresser les fichier SRR
*/

    publishDir "results/FastQ/"
    
    input:
    tuple val(sraid), file(srr) from file_srr
    
    output:
    tuple val(sraid), file("*1.fastq.gz"), file("*2.fastq.gz") into fastqFiles
    
    script:
    """    
    fastq-dump --gzip --split-files ${srr}
    """
}
/* Ce process sert à télécharger lec chromosomes (= génome humain de référence ) */ 

chromosome= Channel.from ("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT","X","Y")

process  chromosomedownloading{ 
/* 	Input: valeurs de chaque chromosome récupéré de la liste "chromosome" ci dessus.
	Output: 25 fichiers sous format n.fa.gz avec n indique la valeur de chaque chromosome. 
	Fonction: Télécharger des fichiers  fasta depuis le lien ci dessous
*/

	publishDir "results/Chrs/"
	input:
	val chromosome from chromosome
	
	output:
	file "${chromosome}.fa.gz" into chromosomefagz
	"""
	wget -O ${chromosome}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chromosome}.fa.gz
	"""
}

/* Ce process sert à fusionner et décompresser les fichiers de chromosome en un fichier ref pour former le génome humain */ 
process CompleteGenome {
/* 	Input: Les 25 chromosomes téléchargés dans le process chromosomedownloading
	Output: 1 fichier HumanGenome.fa regroupant tous les chromosomes avec la fonction collect 
	Fonction: Fusionner  et Décompresser  les 25 fichiers de chromosome pour former le génome humain
*/

    publishDir "results/Genome/"

    input:
    file CompleteChromosomes from chromosomefagz.collect()
		
    output: 
    file "HumanGenome.fa" into human
    
    script:
    """
    gunzip -c ${CompleteChromosomes} > HumanGenome.fa
    """
}
 
/* Ce process sert à créer des index du génome humain */ 
process genomeIndexation {
  /* 	Input:  Le fichier HumanGenome.fa
	Output:  Un fichier ref  contenant un index du génome humain
	Fonction: Indexer le génome humain à partir d'un fichier fasta
*/
    publishDir "results/IndexGenome/"
    
    input: 
    file (Genome) from human.collect()

    output:
    path "ref" into genomeIndex

    script:
    """
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ${Genome}
    """
}

/* Ce process sert à télécharger le fichier annotation pour pouvoir annoter le génome humain et générér un fichier compressé  */ 
process HumanGenomeAnnotation{
  /* 	Input: Ce process ne nécéssite pas un Input 
	Output:  Fichier sous format .gtf(
	Fonction: Annotation et compression du génome humain 
*/
   
   
    publishDir "results/HGA"
    
    output:
    file "HumanGenome.gtf" into gff

    script:
    """
    wget -O HumanGenome.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip -c HumanGenome.gtf.gz > HumanGenome.gtf 
    """
}
/* Ce process sert à aligner les paires de reads sur le génome humain annoté */ 
process mappingfq {
  /* 	Input: les paires de reads récupéres dans le process fastqZipping et le génome humain annoté
	Output:  Fichier sous l'extension.bam
	Function: Alignement des paires de read sur le génome humain annoté en utilisant la méthode STAR
*/
   

    publishDir "results/Mapping/"
    
    input:
    tuple val(sraid), file(SRR_1), file(SRR_2) from fastqFiles
    path IndexGenome from genomeIndex
    
    output:
    file "${sraid}.bam" into bam1, bam2
    
    script:
    """
    STAR --outSAMstrandField intronMotif \
    --outFilterMismatchNmax 4 \
    --outFilterMultimapNmax 10 \
    --genomeDir ${IndexGenome} \
    --readFilesIn <(gunzip -c ${SRR_1}) <(gunzip -c ${SRR_2}) \
    --runThreadN ${task.cpus} \
    --outSAMunmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --outStd BAM_SortedByCoordinate \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM ${task.memory.toBytes()} \
    > ${sraid}.bam
    """
}
/* Ce process sert à indexer les fichiers bam  généré par le process mapping */ 
process generating_bam_files {
    /*
    Input : les fichiers Bam  generés par le process précédent
    Output : Fichiers Bam indexés 
    Function : Indexer les fichiers bam en utilisant le container samtools 
    */
    publishDir "results/bam_files/"
    
    input: 
    file bam from bam1
    
    output:
    file "${bam}.bai" into bam_index_files
    
    script:
    """
    samtools index $bam
    """
}
/* Ce process sert à compter le nombre de read aligné sur le génome humain pour chaque géne et dans chaque échantillon  */ 
process counting_Reads_Matrix {
/*
    Input : collecte des Fichiers bam, des fichiers bam indexés et le génome humain annoté 
    Output : Une matrice de comptage contenant le nombre de read pour chaque géne et chaque échantillons
    Function : Compte le nombre de read alignés sur chaque géne du génome humain
    */
   publishDir "results/featureCounts/"
    
    input:
    file bam_files from bam_index_files.collect()
    file bam from bam2.collect()
    val gtf from gff
    
    output:
    file "featureCounts_Matrix.txt" into Counting_reads
    
    script: 
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a $gtf -o featureCounts_Matrix.txt ${bam}
    """
}
/* Ce process sert à détécter les génes différentiellement exprimés en utilisant le package DESeq2 */ 

meta_Data_channel = Channel.fromPath('./metaData.txt')

process stat_analysis {
/*
    Input : La matrice de comptage et le fichier metadata qui contient des informations sur les 8 échantillons
    Output : Un fichier text contenant les fichiers différentiellement exprimés et les graphes
    Function : Exécution du package DESeq2 sur la matrice de comptage pour trouver les génes différentiellement exprimés entre les patients WT  et mutés. 
    */

    publishDir "results/stat_analysis/"
    
    input:
    file count_Data from Counting_reads
    file meta_Data from meta_Data_channel
    
    output:
    tuple file("DE_gene.txt"), file("Plots*") into channel_end
    
    script:
    """
    Rscript ${workflow.projectDir}/bin/stat_dernier.R ${meta_Data} ${count_Data} "Plots" 
    """
}


workflow.onComplete = {
    println "Sucess"
} 

 
workflow.onError = {
    println "Pipeline Error is : ${workflow.errorMessage}"
}
 

