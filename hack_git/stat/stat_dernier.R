##### Set Seed #####
set.seed(123)
##### install tidyverse
install.packages("tidyverse")
### Bioconductor and CRAN libraries used
library("DESeq2")
library("pheatmap")
library("FactoMineR")
library("factoextra")
library("EnhancedVolcano")
library("tidyverse")
library("RColorBrewer")

################### Loading  Data ################### 
args = commandArgs(TRUE)
meta<-read.table(args[1], header = TRUE, sep = ",",row.names = 1)
data<-read.table(args[2], header = TRUE, sep = "",row.names = 1)
rownames(data) = gsub("[a-zA-Z ]", "", rownames(data))
rownames(data) = as.integer(rownames(data))

### Check classes of the data we just brought in
class(meta)
class(data)

################# Cleaning data #################
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

Hist=ggplot(new_data) +geom_histogram(aes(x =SRR628582), stat = "bin", bins = 200) + xlim(-5, 500)  +xlab("Raw expression counts") + ylab("Number of genes")

############# Modeling count data  for Mutated samples ##################

mean_counts <- apply(new_data[c(1,2,8)], 1, mean)
variance_counts <- apply(new_data[c(1,2,8)], 1, var)
df <- data.frame(mean_counts, variance_counts)

Disp=ggplot(df) + geom_point(aes(x=mean_counts, y=variance_counts)) + geom_line(aes(x=mean_counts, y=mean_counts, color="red")) + scale_y_log10() + scale_x_log10()

############### ACP  ###############
tran_data<-t(new_data)
res_pca<-PCA(tran_data,graph = FALSE)

plot_pca=fviz_pca_ind (res_pca,
              repel = TRUE, # Évite le chevauchement de texte, 
              col.ind = mutations, # colorer par groupe
              palette = c("#8f2c4e", "#30bfbf"),
              legend.title = "Sample Type",
              addEllipses = TRUE,
              ellipse.level=0.95,
              mean.point = FALSE)
###############  Clustering ###############
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = labels2, design = ~ mutations)
dds <- DESeq(dds)
#Transform normalized counts using the rlog transformation #
rld<-rlog(dds,blind=TRUE)
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)
# Plot heatmap
heat.colors <- brewer.pal(6, "Blues")
Heatmap=pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)

##################################### Differencial gene expression with DESeq2 ##################################

#Create DESeq object 
DESeq_object <- DESeqDataSetFromMatrix(countData = new_data, colData = labels2, design = ~ mutations)
#Run analysis 
DESeq_object<-DESeq(DESeq_object)
#### Let's do a quick check of the output DESeq function #### 

## Check the size factors 
sizeFactors(DESeq_object)

## Total number of raw counts per sample
colSums(counts(DESeq_object))

#results DESeq function 
res<-results(DESeq_object)
summary(res)

#### Set thresholds to extract significant differencial expression genes  ###########
p_adj<-0.05
LFC<-0.5

#### Results with the fixed thersholds#####
new_res<-results(DESeq_object, alpha = p_adj, lfcThreshold = LFC)
summary(new_res)


#Remove genes with NA 
res_NO_NA<-new_res[complete.cases(new_res[,6]),]
summary(res_NO_NA)
write.table(res_NO_NA, "DE_gene.txt", sep = " ", row.names = TRUE, col.names = TRUE)
###Store results information in a dataframe object 
res_sig<-res_NO_NA%>%data.frame() %>%rownames_to_column(var="gene") %>% as_tibble()
### Filter results according to p_adj and LFC thresholds
sig_DE_gene <- res_sig %>%filter(padj < p_adj & abs(log2FoldChange) >LFC)
summary(sig_DE_gene)
res_NO_NA<-res_NO_NA[order(res_NO_NA$padj),]
###Extract top 10 signficant genes
top10_sigDE_genes <- sig_DE_gene %>% data.frame() %>%arrange(padj) %>% pull(gene) %>% head(n=10) 
### Extract informations about these top 10 genes
info_top_10_gene<- res_sig %>%filter(padj < p_adj & abs(log2FoldChange) >LFC)  %>%arrange(padj) %>% head(n=10)
#Store theses  informations in a data   frame object 
DF_top_10_gene<-as.data.frame(info_top_10_gene)

##################################### Results visualization ##################################

###Enhances volcano plot: 
lab_italics <- paste0("italic('", rownames(res_NO_NA), "')")
selectLab_italics = paste0("italic('",c(top10_sigDE_genes),"')")

volc=EnhancedVolcano(res_NO_NA,
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
                #parseLabels = TRUE,
                col = c('black', 'pink', 'purple', 'red3'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()
#####Save plot ####
output_name = args[3]

Fig1 = paste(output_name, "Hist_.pdf", sep = "")
pdf(Fig1, width = 10, height = 10)
Hist

Fig2 = paste(output_name, "Disp_.pdf", sep = "")
pdf(Fig2, width = 10, height = 10)
Disp

Fig3 = paste(output_name, "pca_.pdf", sep = "")
pdf(Fig3, width = 10, height = 10)
plot_pca

Fig4= paste(output_name, "Heat.pdf", sep = "")
pdf(Fig4, width = 10, height = 10)
Heatmap
Fig5= paste(output_name, "volc.pdf", sep = "")
pdf(Fig5, width = 10, height = 10)
volc
volcano
