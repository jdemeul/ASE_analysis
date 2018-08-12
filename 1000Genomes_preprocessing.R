## preprocessing

library(readr)
library(parallel)
# library(BSgenome.Hsapiens.1000genomes.hs37d5)
# genome <- BSgenome.Hsapiens.1000genomes.hs37d5
# chrominfo <- genome@seqinfo

## filtered 1000G loci vcf via:
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz
# zcat ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz | grep -P -v "^#" | cut -f1,2,4,5 > ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.txt
# for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do grep -P "^${chr}\t" ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.txt > ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.chr${chr}.txt; done


clean_1000gloci <- function(chr, path = "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_20130501_v5b/", prefix = "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.chr") {
  g1000loci <- read_tsv(file = file.path(path, paste0(prefix, chr, ".txt")), col_names = c("chr", "pos", "ref", "alt"), col_types = "cicc", progress = F)
  g1000loci_dups <- g1000loci[duplicated(g1000loci$pos), "pos"]
  g1000loci <- g1000loci[g1000loci$ref %in% c("A", "C", "T", "G") & g1000loci$alt %in% c("A", "C", "T", "G") & !g1000loci$pos %in% g1000loci_dups, ]
  write_tsv(x = g1000loci, path = file.path(path, paste0(prefix, chr, ".loci.txt")))
  return(NULL)
}

allowed_chroms <- c(1:22, "X")
# debug(clean_1000gloci)
mclapply(X = allowed_chroms, FUN = clean_1000gloci, mc.cores = 6)
