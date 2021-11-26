##################Input Data #####################
library(DESeq2)
sampleData<-read.table("SraRunTable.txt", header = TRUE, sep = ",",row.names = 1)
head(sampleData)
rawCount<-read.table("featureCounts_Matrix.txt", header = TRUE, sep = "",row.names = 1)
head(rawCount)
###################Cleaning Data#################################
new_rawCount= rawCount[,-c(1:5)]
#Remove the extension.bam
colnames(new_rawCount)[c(1:8)] = gsub('.{0,4}$', '', colnames(new_rawCount)[c(1:8)])
#cheking if the order of columns of new_rawcount is the same in sample data (rownames) 
all(rownames(sampleData) == colnames(new_rawCount))
#order it : 
new_rawCount<-new_rawCount[,rownames(sampleData)]
all(rownames(sampleData) == colnames(new_rawCount))
all(rownames(sampleData) == colnames(new_rawCount))
######order sample data#####
labels = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
mutations = c("M", "M", "WT", "WT", "WT", "WT", "WT", "M")
labels2 = data.frame(labels, mutations)

# Order the tissue types so that it is sensible and make sure the control sample is first: WT -> M:muté
labels2$mutations<-factor(labels2$mutations,levels=c("WT","M"))
#Create the DESeqDataSet object: 
deseq2Data <- DESeqDataSetFromMatrix(countData=new_rawCount, colData=labels2, design= ~ mutations)
deseq2Data <- DESeq(deseq2Data)
View(counts(deseq2Data ))
#Pre-filtering of data: 
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) >0, ])
#[1] 25761     8 : un grand nombre quand méme donc je ne pense pas les enlever:
# Si dans l'autre il a eu 35179 parmi 65217 qui ont un  nombre de read >5
####Differential Expression Analysis #### 
#### Estimation of size factors ###
##Check the size factors : 
sizeFactors(deseq2Data)
estimateSizeFactors(deseq2Data)
colSums(counts(deseq2Data))
