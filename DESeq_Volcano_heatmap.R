library(DESeq2)
library (ggplot2)
library (EnhancedVolcano)
library(pheatmap)
library(viridis)

setwd("D:/Seq_Analysis/2023April_Mito_Yng_Old_Novogene/zUMI/UMI_Analysis/OLD")

x <- read.delim("CLK856_mito_exon_umicounts.txt",row.names="Gene")   #read count file

meta <- read.delim("metadata_CLK_mito.txt",row.names=1) #read metafile that has age and genotype for all samples 

dds <- DESeqDataSetFromMatrix(countData = x, colData = meta, design = ~genotype)   #read data into DESeq object   

# remove samples with less than 12 counts in total across all samples and run DESeq
dim(dds)
keep <- rowSums(counts(dds)) >= 6          
dds <- dds[keep,]
dim(dds)

dds = DESeq(dds)

res.opa1OLD <- results(dds, contrast=c("genotype","Opa1","Cas9"), independentFiltering=FALSE)  # compare 2 genotypes using DESeq2

#write.csv(as.data.frame(res.opa1OLD),  file="Deseq2_UMI_CLK856_w1118_OLDvsYoung.csv")    #write csv output

res.opa1OLD$log2FoldChange[res.opa1OLD$log2FoldChange > 5] <- 5
res.opa1OLD$log2FoldChange[res.opa1OLD$log2FoldChange < -5] <- -5
res.opa1OLD$padj[res.opa1OLD$padj < 10e-10] <- 10e-10

#make volcano plot
p <- EnhancedVolcano(res.opa1OLD,
                     lab = NA,
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'Opa1-mut vs Cas9 OLD',
                     pCutoff = 0.05,
                     FCcutoff = 0.4,
                     pointSize = 3.75,
                     col=c('gray50', 'gray30', 'gray40', 'hotpink4'),
                     colAlpha = 1,
                     ylim=c(-1,10),
                     axisLabSize = 46,
                     xlim=c(-6, 6))
p +   ggplot2::coord_cartesian(xlim=c(-6, 6)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 2))
#export as png height - 1200; width - 900

#----Heatmap----
glycolysis <- read.delim("glycolysis.txt")   #read list of genes for which heatmap is to be made
glycolysis <- glycolysis[,1]
ntd <- normTransform(dds)

pheatmap(assay(ntd)[glycolysis,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, scale="row", border_color=NA, color = inferno(100, alpha = 1, begin = 0, end = 1, direction = -1), fontsize=20)

#export as png height - 900; width - 600