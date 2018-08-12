## check called FAT1 variants and overlay with samples having FAT allelic imbalance + overexpression

airesults <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20171120_allelic_imbalance_pooledsamples.txt",
                        sep = "\t", as.is = T)
airesults$sample <- factor(airesults$sample, levels = paste0("WES_", 1:18))
fatsamples <- table(airesults[airesults$gene_name == "FAT1", "sample"])

snvcallfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/Exome_Data/", pattern = "*.raw.snps.indels_VARIANT_FILTRATION.vcf.avinput.hg19_multianno.txt", full.names = T, recursive = T)
snvcallfiles


getFATvariants <- function(snvfiles, WESID) {
  snvcallfile <- grep(pattern = paste0("-WES_", WESID, "."), x = snvfiles, fixed = T, value = T)
  snvcalls <- read.delim(file = snvcallfile, as.is = T, skip = 1)
  colnames(snvcalls) <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene",
                          "genomicSuperDups", "1000g2012apr_all", "esp6500si_all", "esp6500si_aa", "esp6500_ea", "snp129", "snp137", "cg69", "ExAC_Freq",
                          "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "LJB_PhyloP", "LJB_PhyloP_Pred", "LJB_SIFT",
                          "LJB_SIFT_Pred", "LJB_PolyPhen2", "LJB_PolyPhen2_Pred", "LJB_LRT", "LJB_LRT_Pred", "LJB_MutationTaster", "LJB_MutationTaster_Pred",
                          "LJB_GERP++", "clinvar_20140211", "caddgt20", "DMN_Exomes_Freq",
                          "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GT")
  
  fat1calls <- snvcalls[snvcalls$Gene.refGene == "FAT1", ]
  # fat1calls <- snvcalls
  if (nrow(fat1calls) > 0)
    fat1calls$sample <- paste0("WES_", WESID)
  return(fat1calls)
}


allfatcalls <- do.call(rbind, lapply(X = 1:18, FUN = getFATvariants, snvfile = snvcallfiles))
allfatsnvs <- allfatcalls[allfatcalls$snp137 == "NaN", ]

write.table(x = allfatcalls, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20171122_allFAT1_SNVs.txt", quote = F, sep = "\t", row.names = F)
write.table(x = allfatsnvs, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20171122_allFAT1_somaticSNVs.txt", quote = F, sep = "\t", row.names = F)
