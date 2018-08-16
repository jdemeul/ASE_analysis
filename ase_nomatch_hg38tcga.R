## analysis pipeline
library(readr)
library(VariantAnnotation)
library(BSgenome.Hsapiens.NCBI.GRCh38)
bsgnom <- BSgenome.Hsapiens.NCBI.GRCh38
library(biomaRt)
library(ggplot2)
library(rslurm)
library(VGAM)

source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/1000Genomes_getAllelecounts.R")
source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/utils.R")
source(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/ase_nomatch_hg38tcga_functions.R")

## statics
ALLELECOUNTER <- "/srv/sw/eb/software/alleleCount/4.0.0-GCCcore-6.4.0/bin/alleleCounter"
ALLELESDIR <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_dbSNP151_GRCh38/"
JAVA <- "/srv/sw/eb/software/Java/1.8.0_162/bin/java"

EXOMEDIR <- "/srv/shared/vanloo/TCGA/TCGA-LAML/WXS/"
RNADIR <- "/srv/shared/vanloo/TCGA/TCGA-LAML/RNAseq/"
RNAREFGENOME <- "/srv/shared/vanloo/pipeline-files/human/references/alignment/GRCh38.d1.vd1/GRCh38.d1.vd1.fa"

minMapQ <- 35
minBaseQ <- 20

# for cluster
NCORES <- 1

chrs <- c(1:22, "X")
sampledf <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180812_tcgadata_formatted.tsv", as.is = T)



get_ase_tcga_laml <- function(sample_id) {

  TBAMFILE <- file.path(EXOMEDIR, sampledf[sampledf$sample_id == sample_id, "wxs_bam"])
  TRNABAMFILE <- file.path(RNADIR, sampledf[sampledf$sample_id == sample_id, "rna_bam"])

  TOUTDIR <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/TCGA_ASE/", sample_id)
  
  dir.create(path = TOUTDIR)
  
  ## get allelecounts for tumour exome
  mcmapply(FUN = alleleCount,
           locifile = paste0(ALLELESDIR, "1000genomeslocidbSNP151GRCh38_CHR", chrs, ".txt"),
           outfile = file.path(TOUTDIR, paste0(sample_id, "_alleleCounts_chr", chrs, ".txt")),
           MoreArgs = list(bam = TBAMFILE, min_baq = minBaseQ, min_maq = minMapQ),
           mc.cores = NCORES, mc.preschedule = T)
  
  mclapply(X = chrs, FUN = get_alleles_chr_nomatch, allelesdir = ALLELESDIR, countsdir = TOUTDIR, sample_id = sample_id, mc.cores = NCORES, mc.preschedule = T)
  combine_loci_nomatch(countsdir = TOUTDIR, sample_id = sample_id, bsgnom = bsgnom)
  
  ## get allelecounts for tumour transcriptome
  ASEReadCount(hetSNPvcf = file.path(TOUTDIR, paste0(sample_id, "_hetSNPs_nomatch.vcf")),
               bamfile = TRNABAMFILE,
               refgenome = RNAREFGENOME,
               outfile = file.path(TOUTDIR, paste0(sample_id, "_asereadcounts_nomatch.rtable")),
               minBaseQ = minBaseQ, minMapQ = minMapQ)
  
  ## compute statistics
  asedf <- compute_pvals_nomatch(toutdir = TOUTDIR, tsample = sample_id, exclude_bad_snps = F)
  
  ## Some QC and plotting
  # assess filtering
  p2 <- ggplot(data = asedf, mapping = aes(x = pval, fill = filter <= 0.01)) + geom_histogram(binwidth = 0.01) + scale_y_log10()
  ggsave(filename = file.path(TOUTDIR, paste0(sample_id, "_filter.png")), plot = p2, dpi = 300, width = 10, height = 7)
  
  # sum(asedf$padj < .05)
  # ggd.qqplot(asedf[asedf$pval > 0, "pval"])
  # ggd.qqplot(p.adjust(asedf[asedf$pval > 0, "pval"], method = "fdr"))
  p3 <- ggqq(asedf[asedf$pval > 0, "pval"])
  ggsave(filename = file.path(TOUTDIR, paste0(sample_id, "_QQ.png")), plot = p3, dpi = 300, width = 10, height = 7)
  
  ## Gene annotation
  asedf <- merge(x = asedf, y = ase_annotate(asedf = asedf[asedf$pval <= 0.01, ], bsgnom = bsgnom), all.x = T)
  asedf$contig <- factor(x = asedf$contig, levels = c(1:22, "X"))
  asedf <- asedf[order(asedf$contig, asedf$position), ]
  write_tsv(x = asedf, path = file.path(TOUTDIR, paste0(sample_id, "_ase_out.txt")))
  
  # final plot
  p4 <- plot_ase_manhattan(asedf = asedf)
  ggsave(filename = file.path(TOUTDIR, paste0(sample_id, "_manhattan.png")), plot = p4, dpi = 300, width = 10, height = 3)
  
  return(NULL)
}


# debug(get_ase_tcga_laml)
# get_ase_tcga_laml(sample_id = sampledf[1,"sample_id"])

amlasejob <- slurm_apply(f = get_ase_tcga_laml, params = sampledf[,"sample_id", drop = F], jobname = "ase_aml_tcga", nodes = 15, cpus_per_node = 2, add_objects = ls(),
                          pkgs = rev(.packages()), libPaths = .libPaths(), slurm_options = list(), submit = T)
print_job_status(amlasejob)
# cancel_slurm(amlasejob)
