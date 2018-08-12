## analysis pipeline
library(readr)
library(VariantAnnotation)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
bsgenome_hs37d5 <- BSgenome.Hsapiens.1000genomes.hs37d5
library(BSgenome.Hsapiens.UCSC.hg19)
bsgenome_hg19 <- BSgenome.Hsapiens.UCSC.hg19
library(biomaRt)
library(ggplot2)

source(file = "1000Genomes_getAllelecounts.R")
source(file = "utils.R")


## statics
ALLELECOUNTER <- "/home/jdemeul/bin/alleleCounter"
ALLELESDIR <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_20130501_v5b/"
JAVA <- "/srv/sw/eb/software/Java/1.8.0_141/bin/java"


sampledf <- read.delim(file = "20171108_complete_samples", as.is = T)
sampledf <- sampledf[!is.na(sampledf$n_wes_id), ]
# sampledf[4,]

# for (i in 4:11) {
i <- 1
  ## required input
  SAMPLEID <- sampledf[i, "sampleid"]
  TWESID <- paste0("WES_", sampledf[i, "t_wes_id"])
  MWESID <- paste0("WES_", sampledf[i, "n_wes_id"])
  
  print(SAMPLEID)


EXOMEDIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/Exome_Data/"
RNADIR <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/RNA-seq_T-ALL/"
MBAMFILE <- list.files(path = EXOMEDIR, pattern = paste0(MWESID, "_.*bam$"), recursive = T, full.names = T)
TBAMFILE <- list.files(path = EXOMEDIR, pattern = paste0(TWESID, "_.*bam$"), recursive = T, full.names = T)
TRNABAMFILE <- list.files(path = RNADIR, pattern = paste0(SAMPLEID, ".*bam$"), recursive = T, full.names = T)
RNAREFGENOME <- "/srv/shared/vanloo/pipeline-files/human/references/alignment/hg19/ucsc.hg19.fasta"

MOUTDIR <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/", MWESID)
TOUTDIR <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/", TWESID)

minMapQ <- 35
minBaseQ <- 20

NCORES <- 12

## get allelecounts for matched normal exome
chrs <- c(1:22, "X")

dir.create(path = MOUTDIR)
dir.create(path = TOUTDIR)

mcmapply(FUN = alleleCount,
         locifile = paste0(ALLELESDIR, "1000genomesloci2013v5b_chr", chrs, ".txt"),
         outfile = file.path(MOUTDIR, paste0(MWESID, "_alleleCounts_chr", chrs, ".txt")),
         MoreArgs = list(bam = MBAMFILE, min_baq = minBaseQ, min_maq = minMapQ),
         mc.cores = NCORES)


####### read in loci from the matched normal exome, keep only the clean heterozygous ones

## collect allelecounts for ref and alt, keep only those sites exceeding minimal depth requirement
get_alleles_chr <- function(chr, allelesdir, countsdir, outdir, mindepth = 10, sample_id) {
  # allelesfile <- "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_20130501_v5b/1000genomesAlleles2013v5b_chr22.txt"
  # countsfile <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/WES_1/WES_1_alleleCounts_chr22.txt"
  # mindepth <- 10
  allelesfile <- file.path(allelesdir, paste0("1000genomesAlleles2013v5b_chr", chr, ".txt"))
  countsfile <- file.path(countsdir, paste0(sample_id, "_alleleCounts_chr", chr, ".txt"))
  outfile <- file.path(outdir, paste0(sample_id, "_inform_alleles_chr", chr, ".txt"))
  
  alleles <- read_tsv(file = allelesfile, col_names = T, col_types = "cicc")
  counts <- read_tsv(file = countsfile, col_names = c("chr", "pos", "count_A", "count_C", "count_G", "count_T", "depth"), col_types = "ciiiiii", comment = "#")
  alleles$count_ref <- ifelse(alleles$ref == "A", counts$count_A,
                              ifelse(alleles$ref == "C", counts$count_C, 
                                     ifelse(alleles$ref == "G", counts$count_G, counts$count_T)))
  alleles$count_alt <- ifelse(alleles$alt == "A", counts$count_A,
                              ifelse(alleles$alt == "C", counts$count_C, 
                                     ifelse(alleles$alt == "G", counts$count_G, counts$count_T)))
  
  output <- alleles[rowSums(alleles[, c("count_ref", "count_alt")]) >= mindepth, ]
  write_tsv(x = output, path = outfile, col_names = F)
  return(NULL)
}


## further filter allelecounts to retain only confident heterozygous ones, write output into hetSNPs file
## also write loci to feed to AlleleCount for tumour exome, write vcf to feed to ASEReadCounter for tumour RNA
construct_new_loci <- function(countsdir, outdir, alpha = .1, lower = 1/3, mid = .5, upper = 2/3, sample_id, genome = bsgenome_hg19) {
  
  allelecounts_files <- list.files(path = countsdir, pattern = "_inform_alleles_chr", full.names = T)
  allelecounts <- do.call(rbind, lapply(X = allelecounts_files, FUN = read_tsv, col_names = c("chr", "pos", "ref", "alt", "count_ref", "count_alt"), col_types = "ciccii"))
  
  allelecounts <- allelecounts[allelecounts$count_ref > 0 & allelecounts$count_alt > 0, ]
  pconfint <- apply(X = allelecounts[, c("count_ref", "count_alt")], MARGIN = 1, FUN = function(x) qbeta(c(alpha/2, 1-alpha/2), shape1 = x[1]+1, shape2 = x[2]+1, lower.tail = TRUE, log.p = FALSE))
  is_het <- pconfint[1,] >= lower & pconfint[1,] < mid & pconfint[2,] <= upper & pconfint[2,] > mid
  
  outfile <- file.path(outdir, paste0(sample_id, "_hetSNPs.txt"))
  locifile <- file.path(outdir, paste0(sample_id, "_hetSNPs_loci.txt"))
  allelecounts <- allelecounts[is_het, ]
  write_tsv(x = allelecounts, path = outfile, col_names = T)
  write_tsv(x = allelecounts[ , c("chr", "pos")], path = locifile, col_names = F)
  
  locivcf <- file.path(outdir, paste0(sample_id, "_hetSNPs.vcf"))
  ## annoyingly, transcriptome is mapped to hg19, while exomes are mapped to 1000G_GRCh37d5 (positions of main chroms match, but names differ of course)
  allelecounts_vr <- VRanges(seqnames = paste0("chr", allelecounts$chr), ranges = IRanges(start = allelecounts$pos, end = allelecounts$pos), ref = allelecounts$ref, alt = allelecounts$alt, seqinfo = seqinfo(genome))
  allelecounts_vr <- sort(allelecounts_vr)
  sampleNames(allelecounts_vr) <- sample_id
  writeVcf(obj = allelecounts_vr, filename = locivcf)
  
  return(NULL)
}


mclapply(X = chrs, FUN = get_alleles_chr, allelesdir = ALLELESDIR, countsdir = MOUTDIR, outdir = MOUTDIR, sample_id = MWESID, mc.cores = NCORES)
# debug(construct_new_loci)
construct_new_loci(countsdir = MOUTDIR, outdir = MOUTDIR, sample_id = MWESID)


## run allelecount on the new hetloci for tumour
alleleCount(locifile = file.path(MOUTDIR, paste0(MWESID, "_hetSNPs_loci.txt")),
            bam = TBAMFILE,
            outfile = file.path(TOUTDIR, paste0(TWESID, "_hetSNPs_allelecounts.txt")),
            min_baq = minBaseQ, min_maq = minMapQ)



## Use ASEReadCounter to get allelecounts on the RNA

ASEReadCount <- function(hetSNPvcf, bamfile, refgenome, outfile, minBaseQ = 20, minMapQ = 35) {
  cmd <- paste0(JAVA, " -jar /srv/sw/eb/software/GATK/3.8-Java-1.8.0_141/GenomeAnalysisTK.jar",
               " -R ", refgenome,
               " -T ASEReadCounter",
               " -o ", outfile,
               " -I ", bamfile,
               " -sites ", hetSNPvcf,
               " -U ALLOW_N_CIGAR_READS",
               " -minDepth 1",
               " --minMappingQuality ", minMapQ,
               " --minBaseQuality ", minBaseQ,
               " --countOverlapReadsType COUNT_FRAGMENTS_REQUIRE_SAME_BASE")
  system(cmd, wait = T)
  # return(cmd)
}

# RNAREFGENOME <- "/srv/shared/vanloo/pipeline-files/human/references/alignment/hg19/ucsc.hg19.fasta"
# outfile <- file.path(TOUTDIR, paste0(TWESID, "_asereadcounts.rtable"))
# bamfile <- TRNABAMFILE
# hetSNPvcf <- file.path(MOUTDIR, paste0(MWESID, "_hetSNPs.vcf"))
# minMapQ <- 35
# minBaseQ <- 20

ASEReadCount(hetSNPvcf = file.path(MOUTDIR, paste0(MWESID, "_hetSNPs.vcf")),
             bamfile = TRNABAMFILE,
             refgenome = RNAREFGENOME,
             outfile = file.path(TOUTDIR, paste0(TWESID, "_asereadcounts.rtable")),
             minBaseQ = minBaseQ, minMapQ = minMapQ)


## read in loci from the tumour exome, get counts and construct Beta distrib (theta)
## then read in loci from the tumour transcriptome, compute Betabin(X>x|theta)
compute_pvals <- function(toutdir, tsample, filtercutoff = 0.01) {
  asecountsfile <- file.path(toutdir, paste0(tsample, "_asereadcounts.rtable"))
  genomecountsfile <- file.path(toutdir, paste0(tsample, "_hetSNPs_allelecounts.txt"))
  
  asecounts <- read_tsv(file = asecountsfile, col_names = T, col_types = "ciccciiiiiiii")
  genomecounts <- read_tsv(file = genomecountsfile, col_names = c("chr", "pos", "count_A", "count_C", "count_G", "count_T", "depth"), col_types = "ciiiiii", comment = "#")
  
  asecounts$contig <- sub(pattern = "chr", replacement = "", x = asecounts$contig)
  asedf <- merge(x = asecounts, y = genomecounts, by.x = c("contig", "position"), by.y = c("chr", "pos"))
  
  asedf$refCountGenome <- ifelse(asedf$refAllele == "A", asedf$count_A,
                              ifelse(asedf$refAllele == "C", asedf$count_C, 
                                     ifelse(asedf$refAllele == "G", asedf$count_G, asedf$count_T)))
  asedf$altCountGenome <- ifelse(asedf$altAllele == "A", asedf$count_A,
                              ifelse(asedf$altAllele == "C", asedf$count_C, 
                                     ifelse(asedf$altAllele == "G", asedf$count_G, asedf$count_T)))
  asedf <- asedf[ , c("contig", "position", "refAllele", "altAllele", "refCountGenome", "altCountGenome",
                      "refCount", "altCount")]
  
  asedf$filter <- apply(X = asedf[, c("refCountGenome", "altCountGenome", "refCount", "altCount")], MARGIN = 1,
                      FUN = function(x) min(betabinom.test.ab(q = 0, size = x["refCount"] + x["altCount"], shape1 = x["refCountGenome"] + 1, shape2 = x["altCountGenome"] + 1, alternative = "two.sided"),
                                            betabinom.test.ab(q = x["refCount"] + x["altCount"], size = x["refCount"] + x["altCount"], shape1 = x["refCountGenome"] + 1, shape2 = x["altCountGenome"] + 1, alternative = "two.sided")))
  
  asedf$pval <- apply(X = asedf[, c("refCountGenome", "altCountGenome", "refCount", "altCount")], MARGIN = 1,
                FUN = function(x) betabinom.test.ab(q = x["refCount"], size = x["refCount"] + x["altCount"], shape1 = x["refCountGenome"] + 1, shape2 = x["altCountGenome"] + 1, alternative = "two.sided"))
  
  asedf$padj <- 1
  # asedf$is_testworthy <- asedf$filter <= filtercutoff & asedf$refCountGenome > 0 & asedf$altCountGenome > 0
  # asedf$is_testworthy <- asedf$filter <= filtercutoff
  asedf[asedf$filter <= filtercutoff, "padj"] <- p.adjust(asedf[asedf$filter <= filtercutoff, "pval"], method = "fdr")
  
  return(asedf)
}

asedf <- compute_pvals(toutdir = TOUTDIR, tsample = TWESID)


## Some QC

# assess filtering
p2 <- ggplot(data = asedf, mapping = aes(x = pval, fill = filter <= 0.01)) + geom_histogram(binwidth = 0.01) + scale_y_log10()
ggsave(filename = file.path(TOUTDIR, paste0(TWESID, "_filter.png")), plot = p2, dpi = 300, width = 10, height = 7)

# sum(asedf$padj < .05)
# ggd.qqplot(asedf[asedf$pval > 0, "pval"])
# ggd.qqplot(p.adjust(asedf[asedf$pval > 0, "pval"], method = "fdr"))
p3 <- ggqq(asedf[asedf$pval > 0, "pval"])
ggsave(filename = file.path(TOUTDIR, paste0(TWESID, "_QQ.png")), plot = p3, dpi = 300, width = 10, height = 7)



## Gene annotations:

ase_annotate <- function(asedf) {
  chromregions <- paste(asedf$contig, asedf$position, asedf$position, sep = ":", collapse = ",")
  
  # identify and get data from ensembl74 biomart
  # listMarts(host="grch37.ensembl.org")
  # ensembl=useMart(host="www.ensembl.org",biomart = "ENSEMBL_MART_ENSEMBL")
  # listDatasets(ensembl)
  ensembl <- useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  # filters <- listFilters(ensembl)
  # attributes <- listAttributes(ensembl)
  annot <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
                 filters = c("chromosomal_region"), mart = ensembl, values = chromregions)
  
  asegr <- GRanges(seqnames = asedf$contig, ranges = IRanges(start = asedf$position, end = asedf$position),
                   mcols = asedf[, -c(1,2)], seqinfo = seqinfo(bsgenome_hs37d5))
  annotgr <- GRanges(seqnames = annot$chromosome_name, ranges = IRanges(start = annot$start_position, end = annot$end_position),
                     mcols = annot[, -c(1:3)], seqinfo = seqinfo(bsgenome_hs37d5))
  annothits <- findOverlaps(query = asegr, subject = annotgr)
  
  asedf$gene <- NA
  asedf[unique(queryHits(annothits)), "gene"] <- do.call(rbind, by(data = annot[subjectHits(annothits), "external_gene_name"], INDICES = queryHits(annothits),
                                                                   FUN = function(x) data.frame(genes = paste0(unique(x), collapse = ","), stringsAsFactors = F)))
  
  asedf$contig <- factor(x = asedf$contig, levels = c(1:22, "X"))
  asedf <- asedf[order(asedf$contig, asedf$position), ]
  
  return(asedf)
}


asedf <- ase_annotate(asedf = asedf)
write_tsv(x = asedf, path = file.path(TOUTDIR, paste0(TWESID, "_ase_out.txt")))


## Plotting

plot_ase_manhattan <- function(asedf) {
  labeldf <- data.frame(chr = c(1:22,"X"), pos = as.vector(by(data = 1:nrow(asedf), INDICES = asedf$contig, FUN = mean)), brks = as.vector(by(data = 1:nrow(asedf), INDICES = asedf$contig, FUN = max)))
  
  p1 <- ggplot(data = asedf, mapping = aes(x = 1:nrow(asedf), y = -log10(pval))) + geom_point(mapping = aes(colour = contig %in% as.character(seq(2,22,2))), show.legend = F, alpha = .3, size = .5) +
    geom_hline(yintercept = -log10(max(asedf[asedf$padj < .05, "pval"])), linetype = "dashed", colour = "grey") +
    geom_text(data = asedf[asedf$padj <= .05, ], mapping = aes(x = which(asedf$padj <= .05), y = -log10(pval), label = gene), size = 1.5, angle = 45, hjust = 0, nudge_x = nrow(asedf)/250, nudge_y = 0.1, alpha = .5, show.legend = F)
  p1 <- p1 + scale_x_continuous(breaks = labeldf$pos, labels = labeldf$chr, minor_breaks = labeldf$brks) + scale_y_continuous(breaks = seq(0,10,2), oob = scales::squish, limits = c(0,10))
  p1 <- p1 + theme_minimal() + theme(panel.grid.major.x = element_blank()) + labs(x = NULL)
  return(p1)
}

p4 <- plot_ase_manhattan(asedf = asedf)
ggsave(filename = file.path(TOUTDIR, paste0(TWESID, "_manhattan.png")), plot = p4, dpi = 300, width = 10, height = 3)


# }
