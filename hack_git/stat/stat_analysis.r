## Setup
install.packages("tidyverse")
### Bioconductor and CRAN libraries used
library(RColorBrewer)
library(apeglm)
library(dplyr)
library(DESeq2)
library(pheatmap)
library("FactoMineR")
library("factoextra")
library(ggrepel)
library(EnhancedVolcano)
################### Loading  Data ################### 
#### Load in data ###
args = commandArgs(TRUE)
meta<-read.table(args[1], header = TRUE, sep = ",",row.names = 1)
data<-read.table(args[2], header = TRUE, sep = "",row.names = 1)
rownames(data) = gsub("[a-zA-Z ]", "", rownames(data))
rownames(data) = as.integer(rownames(data))
### Check classes of the data we just brought in
class(meta)
class(data)
################# Cleaning data #################"" 
coldata<-meta[,c('sf3b1_mutation_status','disease')]
new_data= data[,-c(1:5)]
colnames(new_data)[c(1:8)] = gsub('.{0,4}$', '', colnames(new_data)[c(1:8)])
### Check that sample names match in both files
all(colnames(new_data) %in% rownames(coldata))
all(colnames(new_data) == rownames(coldata))
###### Match sample names  files in  both files ########
new_data<-new_data[,rownames(coldata)]
all(colnames(new_data) == rownames(coldata))
######Order sample data#####
labels = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
mutations = c("M", "M", "WT", "WT", "WT", "WT", "WT", "M")
labels2 = data.frame(labels, mutations)
# Order the tissue types so that it is sensible and make sure the control sample is first: WT -> M:muté
labels2$mutations<-factor(labels2$mutations,levels=c("WT","M"))

###################### I- RNA-seq count distribution   #########################

ggplot(new_data) +geom_histogram(aes(x =SRR628582), stat = "bin", bins = 200) + xlim(-5, 500)  +xlab("Raw expression counts") + ylab("Number of genes")

############# Modeling count data  for Mutated samples ##################

mean_counts <- apply(new_data[c(1,2,8)], 1, mean)
variance_counts <- apply(new_data[c(1,2,8)], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) + geom_point(aes(x=mean_counts, y=variance_counts)) + geom_line(aes(x=mean_counts, y=mean_counts, color="red")) + scale_y_log10() + scale_x_log10()

#####Count Normalization of data using DESeq2#############
########### 1- Create DESEq2 Object for quality assessment###########
# Create DESeq2Dataset object

dds <- DESeqDataSetFromMatrix(countData = new_data, colData = labels2, design = ~ mutations)
View(counts(dds))
dds$mutation <- factor(dds$mutations, levels = c("WT","M"))
dds2 <- DESeq(dds)

############ Quality assessment and exploratory analysis using DESeq2################

############### Transform normalized counts using the rlog transformation #######

rld<-rlog(dds2,blind=TRUE)

############### PCA ###############

essai<-t(new_data)
res_pca<-PCA(essai)
png(filename = "PCA_plot.png")
fviz_pca_ind (res_pca,
              repel = TRUE, # Évite le chevauchement de texte, 
              col.ind = mutations, # colorer by groups
              palette = c("#8f2c4e", "#30bfbf"),
              legend.title = "Sample Type",
              addEllipses = TRUE,
              ellipse.level=0.95,
              mean.point = FALSE)
###############  Hierarchical Clustering ###############
# Extract the rlog matrix from the object
rld_mat <- assay(rld)

# Compute pairwise correlation values
rld_cor <- cor(rld_mat)

# Plot heatmap

pheatmap(rld_cor)

##################################### Differencial gene expression with DESeq2 ##################################
#Create DESeq object 
DESeq_object <- DESeqDataSetFromMatrix(countData = new_data, colData = labels2, design = ~ mutations)
#Run analysis /
#By re-assigning the results of the function back to the same variable name (dds), we can fill in the slots of our DESeqDataSet object.
DESeq_object<-DESeq(DESeq_object)
#Let's do a quick check of the output DESeq function #### 
## Check the size factors 
sizeFactors(DESeq_object)
#0 verifier cette courbe 
## Check the size factors
sizeFactors(DESeq_object)
## Total number of raw counts per sample
colSums(counts(DESeq_object))
## Total number of normalized counts per sample
colSums(counts(DESeq_object, normalized=T))

#####Model fitting and hypothesis testing #### 
####### Building the results table ########
## Define contrasts, extract results table, and shrink the log2 fold changes
res <- results(DESeq_object)
res<-results(DESeq_object,contrast = c("mutations","M","WT"))

class(res)
summary(res)
mcols(res, use.names=T)

###Summarizing results#####
summary(res)
#Remove genes with NA 
res_NO_NA<-res[complete.cases(res[,6]),]

####Extracting significant differential expressed genes: 
### Set thresholds:
p_adj<-0.05
LFC<-0.58
#order results by adj values : 

###Enhances volcano plot: 

lab_italics <- paste0("italic('", rownames(res_NO_NA), "')")
selectLab_italics = paste0(
  "italic('",
  c(""),
  "')")

EnhancedVolcano(res_NO_NA,
                lab = lab_italics,
                x = 'log2FoldChange',
                y = 'pvalue',
                title='WT vs Mutated with p_adj',
                selectLab = selectLab_italics,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = p_adj,
                FCcutoff = LFC,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'pink', 'purple', 'red3'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()
