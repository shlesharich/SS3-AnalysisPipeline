library(DESeq2)
library (ggplot2)

setwd("D:/Seq_Analysis/2023April_Mito_Yng_Old_Novogene/zUMI/UMI_Analysis/BOTH")

x <- read.delim("CLK856_mito_exon_umicounts.txt",row.names="Gene")   #read count file

meta <- read.delim("metadata_CLK_mito.txt",row.names=1) #read metafile that has age and genotype for all samples 

dds <- DESeqDataSetFromMatrix(countData = x, colData = meta, design = ~Group)   #read data into DESeq object   

# remove samples with less than 12 counts in total across all samples and run DESeq
dim(dds)
keep <- rowSums(counts(dds)) >= 12          
dds <- dds[keep,]
dim(dds)

dds = DESeq(dds)

age <- meta$age
genotype <- meta$genotype

#plot
d <- plotCounts(dds, gene='Ldh', intgroup='age', returnData=TRUE)
names(d)[1] <- "NormalizedCounts"
p <- ggplot(d, aes(x = genotype, y = NormalizedCounts, color=age))+ggtitle("Ldh")+ theme_bw()+theme(text = element_text(color='Grey25', size=18))
p + geom_point(position=position_jitter(w = 0.1,h = 0), size=4)


