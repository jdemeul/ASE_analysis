### reattemt to normalize log2fc values

## get log2-fold expression values of T-ALL samples
library(DESeq2)

countsdir <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_TCGA/htseq_counts/"
# setwd(countsdir)

sampleFiles <- list.files(path = countsdir, pattern = ".htseq.counts")
sampleName <- paste0("htseq_",substr(sampleFiles, 1, 8))
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = "T-ALL")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = countsdir,
                                       design= ~ 1)
# ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) >= 100, ]

dds <- estimateSizeFactors(ddsHTSeq)
dds_vst <- vst(dds, blind = T)
# head(assay(dds_vst), 3)
library(vsn)
meanSdPlot(assay(dds_vst))
# ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))

# hist(assay(dds_vst)[3,])

dds_vstmeans <- rowMeans(x = assay(dds_vst))
vst_fc <- assay(dds_vst) - dds_vstmeans
plot(dds_vstmeans, vst_fc[,1])
plot(log10(rowMeans(counts(dds, normalized=TRUE))+1), vst_fc[,1])

resdf <- as.data.frame(counts(dds, normalized=TRUE))

l2fcdf <- as.data.frame(vst_fc)

gene_ids_names <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20180525_HTSeqCount_gene_names.txt", as.is = T, header = F)
rownames(gene_ids_names) <- substr(gene_ids_names$V1, 1, 15)

l2fcdf$gene_name <- gene_ids_names[substr(rownames(l2fcdf), 1, 15), "V2"]
l2fcdf$mean_expression <- 2^rowMeans(x = log2(resdf+1))

resdf$gene_name <- gene_ids_names[substr(rownames(resdf), 1, 15), "V2"]

write.table(x = l2fcdf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_TCGA/20180731_RNAlog2fc_vst_AML_TCGA.txt", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(x = as.data.frame(assay(dds_vst)), file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_TCGA/20180731_RNAcounts_vst_AML_TCGA.txt", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(x = resdf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_TCGA/20180731_RNAcounts_normalised_AML_TCGA.txt", quote = F, sep = "\t", row.names = T, col.names = T)



### show distrib of TCGA expression data
alloccgeq2 <- unique(read.delim(file = "20180706_l2fc_vst_AIrecurrenceGreaterThan1_table_AML.txt", as.is = T)$gene_name)
l2fcgenesofinterest <- l2fcdf[l2fcdf$gene_name %in% c(alloccgeq2, "PRDM16"), ]

library(reshape2)
library(ggplot2)

l2fcgenesofinterestmelt <- melt(data = l2fcgenesofinterest[, !colnames(l2fcgenesofinterest) == "mean_expression"], id.vars = "gene_name")


p1 <- ggplot(data = l2fcgenesofinterestmelt, mapping = aes(x = gene_name, y = value))
# p1 <- p1 + geom_point(alpha = .5, shape = 16, size = 1.5, stroke = 0)
p1 <- p1 + geom_jitter(colour = "grey", alpha = .6)
p1 <- p1 + geom_violin(scale = "width", alpha =.5)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90)) + labs(y = "log2 fold change in TCGA", x = "gene name")
p1
ggsave(filename = "20180731_l2fc_vst_AIrecurrenceGreaterThan1_AML_TCGA.png", plot = p1, dpi = 300, width = 15, height = 6)
write.table(x = l2fcgenesofinterest, file = "20180731_l2fc_vst_AIrecurrenceGreaterThan1_table_AML_TCGA.txt", sep = "\t", col.names = T, row.names = F, quote = F)



#### GATA2, MECOM, PRDM16, NPM1 expression
l2fcgataprdm <- as.data.frame(t(l2fcdf[l2fcdf$gene_name %in% c("GATA2", "MECOM", "PRDM16", "NPM1"), 1:145]))
colnames(l2fcgataprdm) <- gene_ids_names[substr(colnames(l2fcgataprdm), 1, 15), "V2"]
# l2fcgataprdmmelt <- melt(data = l2fcgataprdm)

par(mfrow = c(1,3))
plot(l2fcgataprdm$PRDM16, l2fcgataprdm$GATA2, xlab = "PRDM16", ylab = "GATA2", main = "log2 fold change")
plot(l2fcgataprdm$PRDM16, l2fcgataprdm$MECOM, xlab = "PRDM16", ylab = "MECOM", main = "log2 fold change")
plot(l2fcgataprdm$GATA2, l2fcgataprdm$MECOM, xlab = "GATA2", ylab = "MECOM", main = "log2 fold change")

