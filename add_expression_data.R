## get log2-fold expression values of T-ALL samples
library(DESeq2)

countsdir <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/countsRNAseq/"
# setwd(countsdir)

sampleFiles <- list.files(path = countsdir, pattern = ".alignments.bam.count.name")
sampleName <- sub(".alignments.bam.count.name", "", sampleFiles)
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = "T-ALL")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = countsdir,
                                       design= ~ 1)
# dds <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) >= 100, ]
dds <- estimateSizeFactors(ddsHTSeq)
resdf <- counts(dds, normalized=TRUE)
medsample <- rowMedians(x = resdf)

l2fcdf <- as.data.frame(log2((resdf+1)/(medsample+1)))

gene_ids_names <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/countsRNAseq/1-JE11333287.alignments.bam.count", as.is = T, header = F)[, c(1,2)]
rownames(gene_ids_names) <- gene_ids_names$V1

l2fcdf$gene_name <- gene_ids_names[rownames(l2fcdf), "V2"]
l2fcdf$median_expression <- medsample

write.table(x = l2fcdf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20171116_RNAlog2fc", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(x = as.data.frame(counts(ddsHTSeq)), file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20171116_RNAcounts", quote = F, sep = "\t", row.names = T, col.names = T)



# 
# library(GenomicFeatures)
# # library(BSgenome.Hsapiens.UCSC.hg19)
# 
# hstxdb <- makeTxDbFromGFF(file = "/srv/shared/vanloo/pipeline-files/human/references/annotation/GENCODE/gencode.v23lift37.annotation.gtf", organism = "Homo sapiens")
# seqlevels(hstxdb) <- sub(pattern = "chr", replacement = "", x = seqlevels(seqinfo(hstxdb)))
# hsexondb <- exons(x = hstxdb, columns = c("gene_id"))
# 
# 
# ase_annotate_local <- function(asedf, exondb) {
#   
#   geneids <- unlist(mcols(exondb)$gene_id)
#   asegr <- GRanges(seqnames = asedf$contig, ranges = IRanges(start = asedf$position, end = asedf$position),
#                    mcols = asedf[, -c(1,2)], seqinfo = seqinfo(bsgenome_hs37d5))
#   
#   annothits <- findOverlaps(query = asegr, subject = exondb)
#   
#   asedf$gene <- NA
#   asedf[unique(queryHits(annothits)), "gene"] <- do.call(rbind, by(data = geneids[subjectHits(annothits)], INDICES = queryHits(annothits),
#                                                                    FUN = function(x) data.frame(genes = paste0(unique(x), collapse = ","), stringsAsFactors = F)))
#   
#   return(asedf)
# }
# 
# ase_results_annot <- ase_annotate_local(asedf = ase_results, exondb = hsexondb)


# get_log2fc <- function(ased)
# ase_results_annotonly <- ase_results[!is.na(ase_results$gene), ]
# splitgenes <- strsplit(x = ase_results_annotonly$gene, split = ",")
# dfreps <- rep(1:nrow(ase_results_annotonly), lengths(splitgenes))
# 
# ase_results_annotonly <- ase_results_annotonly[dfreps, ]
# ase_results_annotonly$gene <- unlist(splitgenes)
# ase_results_annotonly[, c("log2fc", "median_expression")] <- l2fcdf[match(x = ase_results_annotonly$gene, table = l2fcdf$gene_name), c(grep(pattern = "SC540875", x = colnames(l2fcdf)), 53)]



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
  imbalancedf <-  imbalancedf[order(imbalancedf$median_expression, decreasing = F), ]
  labeldf <- data.frame(pos = unlist(lapply(X = 10^(0:4), FUN = function(x) sum(imbalancedf$median_expression < x))), expr = 10^(0:4), stringsAsFactors = F)
  
  outdf_bak <- imbalancedf
  imbalancedf <- imbalancedf[!grepl(pattern = "^HLA.*", x = imbalancedf$gene_name, perl = T), ]
  imbalancedf$notes <- ifelse(imbalancedf$padj > 0.05, "nonsig", 
                        ifelse(imbalancedf$log2fc >= 1, "up",
                               ifelse(imbalancedf$log2fc <= -.73, "down", "nonsig")))
  
  p1 <- ggplot(data = imbalancedf, mapping = aes(x = 1:nrow(imbalancedf), y = -sign(log2fc)*log10(pcombined)))
  p1 <- p1 + geom_point(mapping = aes(colour = notes, size = abs(log2fc)), 
                        show.legend = F, alpha = .4)
  p1 <- p1 + geom_hline(yintercept = c(-1,1)*-log10(max(imbalancedf[imbalancedf$padj < .05, "pcombined"])), linetype = "dashed", colour = "grey") +
    geom_text(data = imbalancedf[imbalancedf$notes != "nonsig", ], mapping = aes(x = which(imbalancedf$notes != "nonsig"), y = -sign(log2fc)*log10(pcombined), label = gene_name), size = 1.5, angle = 45, hjust = 0, nudge_x = nrow(imbalancedf)/250, nudge_y = 0.1, alpha = .5, show.legend = F)
  p1 <- p1 + scale_y_continuous(breaks = seq(-10,10,2), oob = scales::squish, limits = c(-10,10))
  p1 <- p1 + scale_x_continuous(breaks = labeldf$pos, labels = labeldf$expr, name = "median expression (normalised)")
  # p1 <- p1 + scale_color_brewer(type = "div", palette = "RdBu", direction = -1)
  p1 <- p1 + scale_color_manual(values = c(nonsig = "#e0e0e0", up = "#ef8a62", down = "#67a9cf"))
  p1 <- p1 + scale_size_continuous(range = c(1,7.5))
  p1 <- p1 + theme_minimal() + theme(panel.grid.minor.x = element_blank(), axis.text.x = element_text(angle = -90)) + labs(x = NULL)
  return(p1)  
}
### end functions

# generate library of all Hs exons (should be the covered regions).
gtffile <- "/srv/shared/vanloo/pipeline-files/human/references/annotation/GENCODE/gencode.v23lift37.annotation.gtf"
hstxdb <- makeTxDbFromGFF(file = gtffile, organism = "Homo sapiens")
seqlevels(hstxdb) <- sub(pattern = "chr", replacement = "", x = seqlevels(seqinfo(hstxdb)))
hsexondb <- exons(x = hstxdb, columns = c("gene_id"))

sampledf <- read.delim(file = "20171108_complete_samples", as.is = T)


for (i in 1:nrow(sampledf)) {
  # i <- 18
  SAMPLEID <- sampledf[i, "sampleid"]
  TWESID <- paste0("WES_", sampledf[i, "t_wes_id"])
  
  ## read a results file
  ase_resultsfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/nomatch/", TWESID, "/", TWESID, "_ase_out.txt")
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
  l2fcfile <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20171116_RNAlog2fc"
  l2fcdf <-  read.delim(file = l2fcfile, as.is = T)
  outdf[, c("log2fc", "median_expression", "gene_name")] <- l2fcdf[outdf$gene, c(grep(pattern = SAMPLEID, x = colnames(l2fcdf)), 53, 52)]
  
  # format
  outdf$contig <- factor(outdf$contig, levels = c(1:22, "X"))
  outdf <- outdf[order(outdf$contig, as.integer(unlist(lapply(strsplit(outdf$positions, split = ","), FUN = function(x) x[1])))), ]
  outfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/nomatch/", TWESID, "/", TWESID, "_imbalance_expression.txt")
  write.table(x = outdf, file = outfile, quote = F, sep = "\t", row.names = F, col.names = T)
  
  # plot
  p1 <- plot_imbalance_expression(imbalancedf = outdf)
  plotfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/nomatch/", TWESID, "/", TWESID, "_imbalance_expression.png")
  ggsave(filename = plotfile, plot = p1, dpi = 300, width = 15, height = 6)
  
}

