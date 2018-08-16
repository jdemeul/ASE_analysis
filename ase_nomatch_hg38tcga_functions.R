


## collect allelecounts for ref and alt, keep only those sites exceeding minimal depth requirement for each
get_alleles_chr_nomatch <- function(chr, allelesdir, countsdir, mindepth = 3, sample_id) {
  allelesfile <- file.path(allelesdir, paste0("1000genomesAllelesdbSNP151GRCh38_chr", chr, ".txt"))
  countsfile <- file.path(countsdir, paste0(sample_id, "_alleleCounts_chr", chr, ".txt"))
  outfile_allelecounts <- file.path(countsdir, paste0(sample_id, "_inform_alleles_nomatch_chr", chr, ".txt"))
  
  alleles <- read_tsv(file = allelesfile, col_names = T, col_types = "icc")
  counts <- read_tsv(file = countsfile, col_names = c("chr", "pos", "count_A", "count_C", "count_G", "count_T", "depth"), col_types = "ciiiiii", comment = "#")
  alleles$count_ref <- ifelse(alleles$a0 == "A", counts$count_A,
                              ifelse(alleles$a0 == "C", counts$count_C, 
                                     ifelse(alleles$a0 == "G", counts$count_G, counts$count_T)))
  alleles$count_alt <- ifelse(alleles$a1 == "A", counts$count_A,
                              ifelse(alleles$a1 == "C", counts$count_C, 
                                     ifelse(alleles$a1 == "G", counts$count_G, counts$count_T)))
  
  output <- alleles[alleles$count_ref >= mindepth & alleles$count_alt >= mindepth, ]
  write_tsv(x = cbind(chr = rep(chr, nrow(output)), output), path = outfile_allelecounts, col_names = F)
  return(NULL)
}


combine_loci_nomatch <- function(countsdir, sample_id, bsgnom) {
  
  allelecounts_files <- list.files(path = countsdir, pattern = "_inform_alleles_nomatch_chr", full.names = T)
  allelecounts <- do.call(rbind, lapply(X = allelecounts_files, FUN = read_tsv, col_names = c("chr", "pos", "ref", "alt", "count_ref", "count_alt"), col_types = "ciccii"))
  
  outfile <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.txt"))
  write_tsv(x = allelecounts, path = outfile, col_names = T)
  
  locivcf <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.vcf"))
  idxfile <- paste0(locivcf, ".idx")
  ## annoyingly, transcriptome is mapped to hg19, while exomes are mapped to 1000G_GRCh37d5 (positions of main chroms match, but names differ of course)
  allelecounts_vr <- VRanges(seqnames = allelecounts$chr, ranges = IRanges(start = allelecounts$pos, end = allelecounts$pos), ref = allelecounts$ref, alt = allelecounts$alt, seqinfo = seqinfo(bsgnom))
  seqlevelsStyle(allelecounts_vr) <- "UCSC"
  
  allelecounts_vr <- sort(allelecounts_vr)
  sampleNames(allelecounts_vr) <- sample_id
  
  if (file.exists(idxfile))
    file.remove(idxfile)
  writeVcf(obj = allelecounts_vr, filename = locivcf, index = F)
  
  return(NULL)
}


# 
# combine_loci_nomatch_newnonfunc <- function(countsdir, sample_id, chrs = chrs) {
#   
#   allelecounts_files <- paste0(countsdir, "/", sample_id, "_inform_alleles_nomatch_chr", chrs, ".txt")
#   # allelecounts_files <- list.files(path = countsdir, pattern = "_inform_alleles_nomatch_chr", full.names = T)
#   allelecounts <- do.call(rbind, lapply(X = allelecounts_files, FUN = read_tsv, col_names = c("chr", "pos", "ref", "alt", "count_ref", "count_alt"), col_types = "ciccii"))
#   
#   outfile <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.txt"))
#   write_tsv(x = allelecounts, path = outfile, col_names = T)
#   
#   locivcf <- file.path(countsdir, paste0(sample_id, "_hetSNPs_nomatch.vcf"))
#   # file.copy(VCFTEMPLATE, locivcf, overwrite = T)
#   allelecounts_vcf <- data.frame(CHR = paste0("chr", allelecounts$chr), POS = allelecounts$pos, ID = ".", allelecounts[, c("ref","alt")], QUAL = ".", FILTER = ".", INFO = ".", FORMAT = "GT", SAMPLE = "1/1")
#   
#   write_tsv(x = allelecounts_vcf, path = locivcf, append = F, col_names = T)
#   
#   return(NULL)
# }
# 



ASEReadCount <- function(hetSNPvcf, bamfile, refgenome, outfile, minBaseQ = 20, minMapQ = 35) {
  cmd <- paste0(JAVA, " -jar /srv/sw/eb/software/GATK/3.8-1-Java-1.8.0_162/GenomeAnalysisTK.jar",
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

## read in loci from the tumour exome, get counts and construct Beta distrib (theta)
## then read in loci from the tumour transcriptome, compute Betabin(X>x|theta)
compute_pvals_nomatch <- function(toutdir, tsample, filtercutoff = 0.01, exclude_bad_snps = F) {
  asecountsfile <- file.path(toutdir, paste0(tsample, "_asereadcounts_nomatch.rtable"))
  genomecountsfile <- file.path(toutdir, paste0(tsample, "_hetSNPs_nomatch.txt"))
  
  asecounts <- read_tsv(file = asecountsfile, col_names = T, col_types = "ciccciiiiiiii")
  genomecounts <- read_tsv(file = genomecountsfile, col_names = T, col_types = "ciccii")
  colnames(genomecounts) <- c("chr", "pos", "ref", "alt", "refCountGenome", "altCountGenome")
  
  asecounts$contig <- sub(pattern = "chr", replacement = "", x = asecounts$contig)
  asedf <- merge(x = asecounts, y = genomecounts, by.x = c("contig", "position"), by.y = c("chr", "pos"))
  
  asedf <- asedf[ , c("contig", "position", "refAllele", "altAllele", "refCountGenome", "altCountGenome",
                      "refCount", "altCount")]
  
  asedf$filter <- apply(X = asedf[, c("refCountGenome", "altCountGenome", "refCount", "altCount")], MARGIN = 1,
                        FUN = function(x) min(betabinom.test.ab(q = 0, size = x["refCount"] + x["altCount"], shape1 = x["refCountGenome"] + 1, shape2 = x["altCountGenome"] + 1, alternative = "two.sided"),
                                              betabinom.test.ab(q = x["refCount"] + x["altCount"], size = x["refCount"] + x["altCount"], shape1 = x["refCountGenome"] + 1, shape2 = x["altCountGenome"] + 1, alternative = "two.sided")))
  
  asedf$pval <- apply(X = asedf[, c("refCountGenome", "altCountGenome", "refCount", "altCount")], MARGIN = 1,
                      FUN = function(x) betabinom.test.ab(q = x["refCount"], size = x["refCount"] + x["altCount"], shape1 = x["refCountGenome"] + 1, shape2 = x["altCountGenome"] + 1, alternative = "two.sided"))
  
  asedf$padj <- 1
  is_testworthy <- asedf$filter <= filtercutoff
  
  # consider lifting over bad SNPs? If a SNP is bad in hg19, not necessarily in hg38.
  # if (exclude_bad_snps) {
  #   probloci <- read_tsv(file = "/srv/shared/vanloo/pipeline-files/human/references/battenberg/battenberg_problem_loci/probloci_270415.txt.gz", col_types = "ci")
  #   isbadsnp <- paste0(asedf$contig, "_", asedf$position) %in% paste0(probloci$Chr, "_", probloci$Pos)
  #   is_testworthy <- is_testworthy & !isbadsnp
  # }
  
  # asedf$is_testworthy <- asedf$filter <= filtercutoff
  asedf[is_testworthy, "padj"] <- p.adjust(asedf[is_testworthy, "pval"], method = "fdr")
  
  return(asedf)
}


## Gene annotations:

ase_annotate <- function(asedf, bsgnom) {
  chromregions <- paste(asedf$contig, asedf$position, asedf$position, sep = ":", collapse = ",")
  
  # identify and get data from ensembl74 biomart
  # listMarts(host="grch37.ensembl.org")
  # ensembl=useMart(host="www.ensembl.org",biomart = "ENSEMBL_MART_ENSEMBL")
  # listDatasets(ensembl)
  ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  # filters <- listFilters(ensembl)
  # attributes <- listAttributes(ensembl)
  annot <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
                 filters = c("chromosomal_region"), mart = ensembl, values = chromregions)
  
  asegr <- GRanges(seqnames = asedf$contig, ranges = IRanges(start = asedf$position, end = asedf$position),
                   mcols = asedf[, -c(1,2)], seqinfo = seqinfo(bsgnom))
  annotgr <- GRanges(seqnames = annot$chromosome_name, ranges = IRanges(start = annot$start_position, end = annot$end_position),
                     mcols = annot[, -c(1:3)], seqinfo = seqinfo(bsgnom))
  annothits <- findOverlaps(query = asegr, subject = annotgr)
  
  asedf$gene <- NA
  asedf[unique(queryHits(annothits)), "gene"] <- do.call(rbind, by(data = annot[subjectHits(annothits), "external_gene_name"], INDICES = queryHits(annothits),
                                                                   FUN = function(x) data.frame(genes = paste0(unique(x), collapse = ","), stringsAsFactors = F)))
  
  return(asedf)
}


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
