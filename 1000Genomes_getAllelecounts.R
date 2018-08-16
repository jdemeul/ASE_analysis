## allelecounter
# library(parallel)

#' Run alleleCount
#' 
#' Count the alleles for specified locations in the loci file. Expects alleleCount binary in $PATH
#' @param locifile A file with at least chromsome and position columns of the locations to be counted
#' @param bam A bam file
#' @param outfile Where to store the output
#' @param min_baq The minimum base quality required for a read to be counted
#' @param min_maq The minimum mapping quality required for a read to be counted
#' @author sd11
#' @export
alleleCount <- function(locifile, bam, outfile, min_baq=20, min_maq=35) {
  cmd = paste(ALLELECOUNTER,
              "-b", bam,
              "-o", outfile,
              "-l", locifile,
              "-m", min_baq,
              "-q", min_maq,
              "--dense-snps", sep=" ")
  system(cmd, wait=T)
}

# alleleCount(locifile = "/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_20130501_v5b/1000genomesloci2013v5b_chr1.txt",
#             bam = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/Exome_Data/Aligned_20170616-mmansour-EXT2/mmansour-ext-WES_1/mmansour-ext-WES_1_sorted_unique_realigned.bam",
#             outfile = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/WES_1_chr1_alleleCounts.txt")

# ALLELECOUNTER <- "/home/jdemeul/bin/alleleCounter"

# chrs <- c(1:22, "X")
# locifiles <- paste0("/srv/shared/vanloo/pipeline-files/human/references/1000genomes/1000genomes_20130501_v5b/1000genomesloci2013v5b_chr", chrs, ".txt")
# outfiles <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/", MWESID, "/", MWESID, "_alleleCounts_chr", chrs, ".txt")
# bamfile <- paste0("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/Exome_Data/Aligned_20170616-mmansour-EXT2/mmansour-ext-", MWESID, "/mmansour-ext-", MSAMPLE, "_sorted_unique_realigned.bam")
# 
# dir.create(path = dirname(outfiles))
# 
# mcmapply(FUN = alleleCount, locifile = locifiles, outfile = outfiles, MoreArgs = list(bam = bamfile), mc.cores = 12)