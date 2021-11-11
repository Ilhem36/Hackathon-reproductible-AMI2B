Introduction 

The goal of this project is to reproduce parts of the analysis described in these papers :
https://pubmed.ncbi.nlm.nih.gov/23313955/
https://pubmed.ncbi.nlm.nih.gov/23861464/

They performed RNA-Seq in samples from patients with uveal melanoma. Some samples are mutated in SF3B1 which is a splicing factor. To learn more about SF3B1(https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000115524;r=2:197388515-197435079).
nf-uveal melanoma/rnaseq is a bioinformatics pipeline that analyzes this data in order to find diffentially expressed genes. In other words, this pipeline aims to analyze genes that are more (or less ) expressed in one condition(SF3B1 mutated samples)compared to another (SF3B1 non mutated samples). 

 The pipeline in built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It uses Docker/Singularity containers making installation trivial and results highly reproducible. 

