#!/usr/bin/env nextflow
// ====================================================== Main workflow for RNA-seq data analysis ======================================================= // 


File_SRP = Channel.from('SRR628589', 'SRR628588', 'SRR628587', 'SRR628586', 'SRR628585', 'SRR628584', 'SRR628583', 'SRR628582') 
process DownloadSRR {

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


process fastqZipping {

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


chromosome= Channel.from ("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT","X","Y")

process  chromosomedownloading{ 

	publishDir "results/Chrs/"
	input:
	val chromosome from chromosome
	
	output:
	file "${chromosome}.fa.gz" into chromosomefagz
	"""
	wget -O ${chromosome}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chromosome}.fa.gz
	"""
}


process CompleteGenome {

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
 
process genomeIndexation {
  
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


process HumanGenomeAnnotation{
   
    publishDir "results/HGA"
    
    output:
    file "HumanGenome.gtf" into gff

    script:
    """
    wget -O HumanGenome.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip -c HumanGenome.gtf.gz > HumanGenome.gtf 
    """
}
/*
process mappingfq {

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
*/