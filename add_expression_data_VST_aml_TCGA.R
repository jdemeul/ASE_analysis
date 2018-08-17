### reattemt to normalize log2fc values

## get log2-fold expression values of T-ALL samples
library(DESeq2)

countsdir <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/AML_TCGA/htseq_counts/"
annotsfile <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/AML_TCGA/annotations/gdc_sample_sheet.2018-07-30-2.tsv"
# setwd(countsdir)

annots <- read.delim(file = annotsfile, as.is = T)

sampleFiles <- sub(pattern = ".gz$", replacement = "", x = annots$File.Name)
  # list.files(path = countsdir, pattern = ".htseq.counts")
sampleName <- annots$Sample.ID
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
# library(vsn)
# meanSdPlot(assay(dds_vst))
# ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))

# hist(assay(dds_vst)[3,])

dds_vstmeans <- rowMeans(x = assay(dds_vst))
vst_fc <- assay(dds_vst) - dds_vstmeans
plot(dds_vstmeans, vst_fc[,1])
plot(log10(rowMeans(counts(dds, normalized=TRUE))+1), vst_fc[,1])

resdf <- as.data.frame(counts(dds, normalized=TRUE))

l2fcdf <- as.data.frame(vst_fc)

gene_ids_names <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180525_HTSeqCount_gene_names.txt", as.is = T, header = F)
rownames(gene_ids_names) <- substr(gene_ids_names$V1, 1, 15)

l2fcdf$gene_name <- gene_ids_names[substr(rownames(l2fcdf), 1, 15), "V2"]
l2fcdf$mean_expression <- 2^rowMeans(x = log2(resdf+1))

resdf$gene_name <- gene_ids_names[substr(rownames(resdf), 1, 15), "V2"]

write.table(x = l2fcdf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/AML_TCGA/20180731_RNAlog2fc_vst_AML_TCGA.txt", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(x = as.data.frame(assay(dds_vst)), file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/AML_TCGA/20180731_RNAcounts_vst_AML_TCGA.txt", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(x = resdf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/AML_TCGA/20180731_RNAcounts_normalised_AML_TCGA.txt", quote = F, sep = "\t", row.names = T, col.names = T)



### show distrib of TCGA expression data
alloccgeq2 <- unique(read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180706_l2fc_vst_AIrecurrenceGreaterThan1_table_AML.txt", as.is = T)$gene_name)
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
ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180731_l2fc_vst_AIrecurrenceGreaterThan1_AML_TCGA.png", plot = p1, dpi = 300, width = 15, height = 6)
write.table(x = l2fcgenesofinterest, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180731_l2fc_vst_AIrecurrenceGreaterThan1_table_AML_TCGA.txt", sep = "\t", col.names = T, row.names = F, quote = F)



#### GATA2, MECOM, PRDM16, NPM1 expression
l2fcgataprdm <- as.data.frame(t(l2fcdf[l2fcdf$gene_name %in% c("GATA2", "MECOM", "PRDM16", "NPM1"), 1:145]))
colnames(l2fcgataprdm) <- gene_ids_names[substr(colnames(l2fcgataprdm), 1, 15), "V2"]
# l2fcgataprdmmelt <- melt(data = l2fcgataprdm)

par(mfrow = c(1,3))
plot(l2fcgataprdm$PRDM16, l2fcgataprdm$GATA2, xlab = "PRDM16", ylab = "GATA2", main = "log2 fold change")
plot(l2fcgataprdm$PRDM16, l2fcgataprdm$MECOM, xlab = "PRDM16", ylab = "MECOM", main = "log2 fold change")
plot(l2fcgataprdm$GATA2, l2fcgataprdm$MECOM, xlab = "GATA2", ylab = "MECOM", main = "log2 fold change")





###### actual use of the TCGA ASE output

## combine p-values (of powered SNP loci) per gene
library(GenomicFeatures)
library(ggplot2)

### functions

combine_pvals <- function(ase_pergene) {
  # ase_pergene <- ase_results_annot[4:14, ]
  outdf <- data.frame(contig = ase_pergene[1, "contig"], positions = paste0(unique(ase_pergene$position), collapse = ","), pcombined = 1, gene = ase_pergene[1, "gene"], stringsAsFactors = F, power = T)
  
  is_duplicated <- duplicated(ase_pergene$position)
  has_power <- ase_pergene$filter <= 0.01
  ase_pergene <- ase_pergene[!is_duplicated & has_power, ]
  
  if (nrow(ase_pergene) > 1) {
    outdf$pcombined <- fishersMethod(x = ase_pergene$pval)
  } else if (nrow(ase_pergene) == 1) {
    outdf$pcombined <- ase_pergene$pval
  } else {
    outdf$power <- F
  }
  return(outdf)
}


fishersMethod <- function(x) {
  pchisq(q = -2 * sum(log(x)), df = 2*length(x), lower.tail = F)
}


plot_imbalance_expression <- function(imbalancedf) {
  imbalancedf <-  imbalancedf[order(imbalancedf$mean_expression, decreasing = F), ]
  labeldf <- data.frame(pos = unlist(lapply(X = 10^(0:4), FUN = function(x) sum(imbalancedf$mean_expression < x))), expr = 10^(0:4), stringsAsFactors = F)
  
  outdf_bak <- imbalancedf
  imbalancedf <- imbalancedf[!grepl(pattern = "^HLA.*", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^IG[HLK].*", x = imbalancedf$gene_name, perl = T) &
                               !grepl(pattern = "^TR[ABDG][VCDJ].*", x = imbalancedf$gene_name, perl = T), ]
  
  imbalancedf$notes <- ifelse(imbalancedf$padj > 0.05, "nonsig", 
                              ifelse(imbalancedf$log2fc >= 1, "up",
                                     ifelse(imbalancedf$log2fc <= -.73, "down", "nonsig")))
  
  p1 <- ggplot(data = imbalancedf, mapping = aes(x = 1:nrow(imbalancedf), y = -sign(log2fc)*log10(pcombined)))
  p1 <- p1 + geom_point(mapping = aes(colour = notes, size = abs(log2fc)), 
                        show.legend = F, alpha = .4)
  p1 <- p1 + geom_hline(yintercept = c(-1,1)*-log10(max(imbalancedf[imbalancedf$padj < .05, "pcombined"])), linetype = "dashed", colour = "grey") +
    geom_text(data = imbalancedf[imbalancedf$notes != "nonsig", ], mapping = aes(x = which(imbalancedf$notes != "nonsig"), y = -sign(log2fc)*log10(pcombined), label = gene_name), size = 1.5, angle = 45, hjust = 0, nudge_x = nrow(imbalancedf)/250, nudge_y = 0.1, alpha = .5, show.legend = F)
  p1 <- p1 + scale_y_continuous(breaks = seq(-10,10,2), oob = scales::squish, limits = c(-10,10))
  p1 <- p1 + scale_x_continuous(breaks = labeldf$pos, labels = labeldf$expr, name = "mean expression (normalised)")
  # p1 <- p1 + scale_color_brewer(type = "div", palette = "RdBu", direction = -1)
  p1 <- p1 + scale_color_manual(values = c(nonsig = "#e0e0e0", up = "#ef8a62", down = "#67a9cf"))
  p1 <- p1 + scale_size_continuous(range = c(1,7.5))
  p1 <- p1 + theme_minimal() + theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = -90)) + labs(x = NULL)
  return(p1)  
}
### end functions





# generate library of all Hs exons (should be the covered regions).
gtffile <- "/srv/shared/vanloo/pipeline-files/human/references/annotation/GDC.h38/gencode.v22.annotation.gtf.gz"
hstxdb <- makeTxDbFromGFF(file = gtffile, organism = "Homo sapiens")
# seqlevels(hstxdb) <- sub(pattern = "chr", replacement = "", x = seqlevels(seqinfo(hstxdb)))
hsexondb <- exons(x = hstxdb, columns = c("gene_id"))
seqlevelsStyle(hsexondb) <- "Ensembl"

annotsfile <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/AML_TCGA/annotations/gdc_sample_sheet.2018-07-30-2.tsv"
annots <- read.delim(file = annotsfile, as.is = T)

# add in the log2-fold change data and actual gene names
l2fcfile <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/AML_TCGA/20180731_RNAlog2fc_vst_AML_TCGA.txt"
l2fcdf <-  read.delim(file = l2fcfile, as.is = T)

for (SAMPLEID in annots$Sample.ID) {
  # i <- 1
  # SAMPLEID <- sampledf[i, "Number"]
  
  ## read a results file
    # TWESID <- paste0("WES_", sampledf[i, "t_wes_id"])
  ase_resultsfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/TCGA_ASE/", SAMPLEID, "/", SAMPLEID, "_ase_out.txt")
  if (!file.exists(ase_resultsfile))
    next
  
  ase_results <- read.delim(file = ase_resultsfile, as.is = T)
  
  
  # make results into GRanges object, identify all exonic SNPs and create new df with all of these (contains duplicate SNPs)
  asegr <- GRanges(seqnames = ase_results$contig, ranges = IRanges(start = ase_results$position, end = ase_results$position))
  annothits <- findOverlaps(query = asegr, subject = hsexondb)
  # in one case, there were two genes using the same exon ... this just takes the first
  hitgenes <- sapply(mcols(hsexondb[subjectHits(annothits)])$gene_id, FUN = function(x) x[[1]])
  ase_results_annot <- data.frame(ase_results[queryHits(annothits), colnames(ase_results) != "gene"], gene = hitgenes, stringsAsFactors = F)
  
  
  # create output dataframe with combined p-value per gene + adjust for multiple testing
  outdf <- do.call(rbind, by(data = ase_results_annot, INDICES = ase_results_annot$gene, FUN = combine_pvals))
  outdf$padj <- 1
  outdf[outdf$power, "padj"] <- p.adjust(p = outdf[outdf$power, "pcombined"], method = "fdr")
  # outdf$gene_name <- gene_ids_names[outdf$gene, "V2"]
  
  # add in the log2-fold change data and actual gene names
  outdf[, c("log2fc", "mean_expression", "gene_name")] <- l2fcdf[outdf$gene, c(gsub(pattern = "-", replacement = ".", x = SAMPLEID), "mean_expression", "gene_name")]
  
  # format
  outdf$contig <- factor(outdf$contig, levels = c(1:22, "X"))
  outdf <- outdf[order(outdf$contig, as.integer(unlist(lapply(strsplit(outdf$positions, split = ","), FUN = function(x) x[1])))), ]
  outfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/TCGA_ASE/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.txt")
  write.table(x = outdf, file = outfile, quote = F, sep = "\t", row.names = F, col.names = T)
  
  
  # plot
  p1 <- plot_imbalance_expression(imbalancedf = outdf)
  plotfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/TCGA_ASE/", SAMPLEID, "/", SAMPLEID, "_imbalance_expression_vst.png")
  ggsave(filename = plotfile, plot = p1, dpi = 300, width = 15, height = 6)
}

