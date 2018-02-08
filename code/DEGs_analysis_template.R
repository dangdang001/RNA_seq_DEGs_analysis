
################################################################################################
######## Differential Expression Analysis with DESeq2 (Based on the count table) ######
################################################################################################
### Created by Donglei Yin
### Date: 2018-01-29
### Project: KLF5-knockout RNA-seq (Mouse)
#####################################################################


setwd("C:/Users/dyin/Documents/RNA-seq project-20171213")


##############################################################
############     Load Required R Libraries   #################
##############################################################

list.of.packages <- c("devtools", "Seurat", "dplyr", "Matrix", "reshape", "ggplot2","stats","pheatmap","digest")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("AnnotationDbi")
biocLite("EnsDb.Mmusculus.v79")

library(devtools)
library(Matrix)
library(Seurat)
library(dplyr)
library(DESeq2)
library(reshape)
library(ggplot2)
library(stats)
library(pheatmap)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
library(digest)

##############################################################
##Prepare data file, 4 treatment vs 4 controls##
##############################################################
Ctrl_1=read.table("./data/Sample_Ctrl-1/analysis/Ctrl-1.counts.txt",header=FALSE,skip=2)
Ctrl_2=read.table("./data/Sample_Ctrl-2/analysis/Ctrl-2.counts.txt",header=FALSE,skip=2)
Ctrl_3=read.table("./data/Sample_Ctrl-3/analysis/Ctrl-3.counts.txt",header=FALSE,skip=2)
Ctrl_4=read.table("./data/Sample_Ctrl-4/analysis/Ctrl-4.counts.txt",header=FALSE,skip=2)
KO_1=read.table("./data/Sample_KO-1/analysis/KO-1.counts.txt",header=FALSE,skip=2)
KO_2=read.table("./data/Sample_KO-2/analysis/KO-2.counts.txt",header=FALSE,skip=2)
KO_3=read.table("./data/Sample_KO-3/analysis/KO-3.counts.txt",header=FALSE,skip=2)
KO_4=read.table("./data/Sample_KO-4/analysis/KO-4.counts.txt",header=FALSE,skip=2)

colnames(Ctrl_1)<-c("Geneid","Chr","Start","End","Strand","Length","Ctrl_1")
colnames(Ctrl_2)<-c("Geneid","Chr","Start","End","Strand","Length","Ctrl_2")
colnames(Ctrl_3)<-c("Geneid","Chr","Start","End","Strand","Length","Ctrl_3")
colnames(Ctrl_4)<-c("Geneid","Chr","Start","End","Strand","Length","Ctrl_4")
colnames(KO_1)<-c("Geneid","Chr","Start","End","Strand","Length","KO_1")
colnames(KO_2)<-c("Geneid","Chr","Start","End","Strand","Length","KO_2")
colnames(KO_3)<-c("Geneid","Chr","Start","End","Strand","Length","KO_3")
colnames(KO_4)<-c("Geneid","Chr","Start","End","Strand","Length","KO_4")

Ctrl_1<-Ctrl_1[,c("Geneid","Ctrl_1")];Ctrl_2<-Ctrl_2[,c("Geneid","Ctrl_2")];Ctrl_3<-Ctrl_3[,c("Geneid","Ctrl_3")];Ctrl_4<-Ctrl_4[,c("Geneid","Ctrl_4")];
KO_1<-KO_1[,c("Geneid","KO_1")];KO_2<-KO_2[,c("Geneid","KO_2")];KO_3<-KO_3[,c("Geneid","KO_3")];KO_4<-KO_4[,c("Geneid","KO_4")];

data.all<-Reduce(function(x,y) merge(x,y,by="Geneid",all=TRUE),list(Ctrl_1,Ctrl_2,Ctrl_3,Ctrl_4,KO_1,KO_2,KO_3,KO_4))

rownames(data.all)<-data.all$Geneid
data.all$Geneid<-NULL
rownames(data.all)<-sub('\\..*', '', rownames(data.all))

####################################################################################################################################
########################## The follwing part depends on which samples you want to compare############################################
####################################################################################################################################

################################################
# Samples used in this template: 

# Treatment:KO_1,KO_2,KO_3
# Control:Ctrl_1,Ctrl_2,Ctrl_3,Ctrl_4

# Contents:
    # 1. Generate DE gene lists
    # 2. PCA plot
    # 3. Heatmap (both samples and genes)
    # 4. Volcano Plot

################################################

# Before run the following steps, make sure you have created a folder called "Results_KOxxxx_vs_Ctrlxxxx" under your current working directory (step 1)

# For example, if you want to compare KO 1,2,4 with Ctrl 1,2,3,4, then name it as "Results_KO124_vs_Ctrl1234", so that corresposiding results will be saved under the same folder.

# Here is the example comparing KO 1,2,3 with Ctrl 1,2,3,4

setwd("./Results_KO123_vs_Ctrl1234")   ## Sample Parameter


#######################################1. Generate DE gene lists###########################################

############################################
##Convert count data into a DESeq2 object###
############################################


data.selected=data.all[c("Ctrl_1","Ctrl_2","Ctrl_3","Ctrl_4","KO_1","KO_2","KO_3")]  
## Sample Parameter, remove or add other samples interested

coldata<-data.frame(condition=c(rep("Control",4),rep("KLF5_knockout",3)))  
## Sample Parameter, change the two numbers(#samples in each of the two groups), we have 4 control samples and 3 KO samples for now

rownames(coldata)<-c("Ctrl_1","Ctrl_2","Ctrl_3","Ctrl_4","KO_1","KO_2","KO_3")  
## Sample Parameter, remove or add other samples interested

# make sure the orders are consistent, True is correct
all(rownames(coldata) == colnames(data.selected))

dds <- DESeqDataSetFromMatrix(countData = data.selected,
                              colData = coldata,
                              design = ~ condition)
dds #45232 

#####################################################################
##Perform preliminary filtering, trim out extremely low counts read##
#####################################################################

dds <- dds[ rowSums(counts(dds)) >= 3, ] 
dds
#25954

##############################################################
##DESeq differential expression analysis result object #######
##############################################################
dds <- DESeq(dds)
deResult <- results(dds) 
deResult

##############################################################
##Analyse DE analysis results#################################
##############################################################


# Annotating (Ensembl gene IDs to gene symbols)

library("AnnotationDbi")
library(EnsDb.Mmusculus.v79)



columns(EnsDb.Mmusculus.v79)
keytypes(EnsDb.Mmusculus.v79)

rownames(dds) <- mapIds(EnsDb.Mmusculus.v79,
                        keys=rownames(dds),
                        column="SYMBOL",
                        keytype="GENEID",
                        multiVals="first")


# Order by p values, signifiant level was set to be 0.05, you can change it to any other level interested by changing alpha

alpha = 0.05
deResultOrdered <- deResult[order(deResult$padj),]
# Select rows that p value < alpha
deResultSig <- deResult[(!is.na(deResult$padj)) & (deResult$padj < alpha),]
# Statistical summary of DE results
summary(deResult)
summary(deResultSig)
# Number of significant genes under alpha
sum(deResult$padj < alpha, na.rm=TRUE) #2940

deResultOrdered$symbol <- mapIds(EnsDb.Mmusculus.v79,
                                 keys=row.names(deResultOrdered),
                                 column="SYMBOL",
                                 keytype="GENEID",
                                 multiVals="first")



#########################################
## Export result to csv #################
#########################################
write.csv(as.data.frame(deResultOrdered),
          file="./markers_DE_0.05.csv")


#######################################2. PCA plot ###########################################

# log-transformed data
rld <- rlog(dds, blind = FALSE)

# PCA plot

source("./../code/pca.R") # this file includes an R function I wrote for this use

# Changing i,j(1<=i<j<=4) to view PCA plot of PCi vs PCj
# Might need to tune the parameters "xlimits" to make sure all the samples were displayed

png("./PCA_P1vsP2.png",height=6, width=8,units='in',res = 600)
plotPCA.dy(rld, intgroup=c("condition"),i=1,j=2,xlimits=c(-25,25))
dev.off()

png("./PCA_P2vsP3.png",height=6, width=8,units='in',res = 600)
plotPCA.dy(rld, intgroup=c("condition"),i=2,j=3,xlimits=c(-25,25))
dev.off()

png("./PCA_P3vsP4.png",height=6, width=8,units='in',res = 600)
plotPCA.dy(rld, intgroup=c("condition"),i=3,j=4,xlimits=c(-10,20))
dev.off()


#######################################3. Heatmap ##############################################

# 3.1 Sample distances heatmap

library("pheatmap")
library("RColorBrewer")

sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rownames(colData(rld)), sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png("./sample_distcance_heatmap.png",height=6, width=8,units='in',res = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

# 3.2 Gene (20 genes with the highest variance across samples) heatmap

# In the sample distance heatmap made previously, the dendrogram at the side shows us a hierarchical clustering of the samples.
# Such a clustering can also be performed for the genes. Since the clustering is only relevant for genes that actually carry 
# a signal, one usually would only cluster a subset of the most highly variable genes. Here, for demonstration, 
# let us select the 20 genes with the highest variance across samples. We will work with the rlog transformed counts:


rld <- rlog(dds, blind = FALSE)

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 30)

mat  <- assay(rld)[topVarGenes, ]
mat  <- mat - rowMeans(mat)

mat<-mat[!is.na(row.names(mat)),][1:20,]

anno <- data.frame(Condition=colData(rld)[, c("condition")])

rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno,cluster_rows=FALSE)

png("./heatmap_20markers.png",height=6, width=8,units='in',res = 600)
pheatmap(mat, cluster_rows=FALSE)
dev.off()

###########################################4. Volcano Plot ################################################

# Genes with FDR<0.05 were marked as red 
# Genes with padj<0.001,abs(log2FoldChange)>3 were labeled with gene symbols

library(ggrepel)

results<-data.frame(Gene=deResultOrdered$symbol,log2FoldChange=deResultOrdered$log2FoldChange, pvalue=deResultOrdered$pvalue,padj=deResultOrdered$padj)

results = results[!is.na(results$padj),]

results = mutate(results, Group=ifelse(results$padj<0.05, "FDR<0.05", "Not Significant"))

p = ggplot(results, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=Group)) +
  scale_color_manual(values=c("red", "black"))

p<-p+geom_text_repel(data=filter(results, padj<0.001,abs(log2FoldChange)>3), aes(label=Gene))

png("./volcano_plot.png",height=6, width=8,units='in',res = 600)
p
dev.off()

