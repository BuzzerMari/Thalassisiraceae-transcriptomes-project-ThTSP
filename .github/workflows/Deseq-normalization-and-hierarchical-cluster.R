library(DESeq2)
#read table with gene counts
x <- read.table("Gene2_count_all_ThTSP",row.names="Gene",sep = "\t",header=TRUE)
#COndition file contains the origin of samples
cond <- read.table("condition.txt",header = TRUE,sep = "\t" )
cond$condition <- factor(cond$condition)
dds <- DESeqDataSetFromMatrix(countData = x, colData = cond, design = ~ condition)
dds <- estimateSizeFactors(dds)

# DESeq2 normalization counts
y = counts(dds, normalized = TRUE)
head(y)
sizeFactors(dds)
#save normalized table output and size factors
write.table(y, file = "norm_deseq_ThalassiosiralestRANSgenes2", sep = "\t",
            row.names = TRUE)
size<-sizeFactors(dds)
write.table(size, file = "size_factor_deseq_ThalassiosiralestRANSgenes", sep = "\t",
            row.names = TRUE)
         
install.packages("ggfortify")
library(ggfortify)
#rlog normalization of dds normalized counts 
rld <- rlog(dds, blind=FALSE)
head(assay(rld),3)
ntd <- normTransform(dds)
BiocManager::install("vsn")
library("vsn")
meanSdPlot(assay(rld))
# calculate sample dists
sampleDists <- dist(t(assay(rld)))
# libraries forhierachical cluster visualization
library("RColorBrewer")
library("viridis")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$sample, rld$condition, sep="-")
colnames(sampleDistMatrix) <- paste(dds$sampler, rld$condition, sep="-")
colors= viridis(16)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         show_rownames=T, show_colnames=T,
         col=colors, cellwidth=24, cellheight = 24

         )
         
         
