### PCA plot AML17 with annotation of the PRDM high samples

vsd <- read.delim(file = "20180706_RNAlog2fc_vst_AML.txt", as.is = T)

# vsd2 <- assay(dds_vst)
vsd_prcomp <- prcomp(x = t(vsd2), retx = T, center = T, scale. = T)

plotdata <- as.data.frame(vsd_prcomp$x)
plotdata$prdm_hi <- c(vsd[vsd$gene_name == "PRDM16", 1:16] > 1)
plotdata$mecom_hi <- c(vsd[vsd$gene_name == "MECOM", 1:16] > 1)

library(ggplot2)
p1 <- ggplot(data = plotdata, mapping = aes(x = PC1, y = PC2)) + geom_point(mapping = aes(color = interaction(prdm_hi, mecom_hi)))
p1

des2pca <- plotPCA(dds_vst, returnData = T)
des2pca$prdmhi <- c(vsd[vsd$gene_name == "PRDM16", 1:16] > 1)
des2pca$mecomhi <- c(vsd[vsd$gene_name == "MECOM", 1:16] > 1)
percentVar <- round(100 * attr(des2pca, "percentVar"))
ggplot(des2pca, aes(PC1, PC2, colour = interaction(prdmhi, mecomhi))) +
  geom_text(mapping = aes(label = name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

  
  
### check for chimeric reads deriving from PRDM16 or LINC00982

library(GenomicRanges)

chimreadsfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/AML/RNA/mapped/", pattern = "*_Chimeric.out.junction", full.names = T, recursive = T)
chimeric_reads <- lapply(X = chimreadsfiles, FUN = read.delim, as.is = T, header = F, col.names = c("chr1", "start1", "strand1", "chr2", "start2", "strand2", "junction", "replenleft", "replenright", "read_name", "base1seg1", "CIGARseg1", "base1seg2", "CIGARse2"))
chimeric_reads_oi <- lapply(chimeric_reads, FUN = function(x) x[(x$chr1 == "chr1" & x$start1 > 2950000 & x$start1 < 3370000) | 
                                                                  (x$chr2 == "chr1" & x$start2 > 2950000 & x$start2 < 3370000), ])

IDs <- sub(basename(chimreadsfiles), pattern = "_Chimeric.out.junction", replacement = "")
chimeric_reads_oi_df <- do.call(rbind, chimeric_reads_oi)
chimeric_reads_oi_df$sample <- rep(IDs, sapply(chimeric_reads_oi, nrow))
chimeric_reads_oi_df$chr1pos <- ifelse(chimeric_reads_oi_df$chr1 == "chr1" & chimeric_reads_oi_df$start1 > 2950000 & chimeric_reads_oi_df$start1 < 3370000, chimeric_reads_oi_df$start1, chimeric_reads_oi_df$start2)
chimeric_reads_oi_df$chrmate <- ifelse(chimeric_reads_oi_df$chr1 == "chr1" & chimeric_reads_oi_df$start1 > 2950000 & chimeric_reads_oi_df$start1 < 3370000, chimeric_reads_oi_df$chr2, chimeric_reads_oi_df$chr1)


library(ggplot2)

p1 <- ggplot(data = chimeric_reads_oi_df, mapping = aes(x = chr1pos, y = sample, colour = chrmate)) + geom_point(shape = "|", alpha = .8, size = 4, position = position_jitter(width = 1000, height = 0))
p1 <- p1 + theme_minimal()
p1

p1 <- ggplot(data = chimeric_reads_oi_df, mapping = aes(x = chr1pos)) + geom_histogram(aes(fill = sample), bins = 250)
p1




aseoutfile <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_ase/22/22_ase_out.txt"
aseoutfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_ase/", pattern = "*_ase_out.txt", full.names = T, recursive = T)
lapply(aseoutfiles, plotfile)
# plotfile(aseoutfile = aseoutfile)

plotfile <- function(aseoutfile) {
asesample <- sub(basename(aseoutfile), pattern = "_ase_out.txt", replacement = "")
aseout <- read.delim(file = aseoutfile, as.is = T)
aseout <- aseout[aseout$contig == "1" & aseout$position > 1000000 & aseout$position < 5000000,]
p1 <- ggplot(data = aseout, mapping = aes(x = position)) + geom_point(mapping = aes(y = altCountGenome/(altCountGenome+refCountGenome))) +
  geom_point(data = aseout[aseout$filter < .001,], mapping = aes(y = altCount/(altCount+refCount)), colour = "green") +
  theme_minimal() + labs(x = "chr1 position", y = "B-allele frequency") + ylim(c(0,1)) + xlim(c(1000000, 5000000))
ggsave(filename = paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_PRDM16_UPD/", asesample, "_BAFplot.png"))
}
