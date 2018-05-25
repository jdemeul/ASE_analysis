## check allelically imbalanced and downregulated genes for nonsense mutations

airesults <- read.delim(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20180525_allelic_imbalance_pooledsamples_vst.txt",
                        sep = "\t", as.is = T)
airesults$sample <- factor(airesults$sample, levels = c(paste0("WES_", 1:18), sampledf[sampledf$cell_line, "sampleid"]))

airecurrences <- table(airesults[airesults$log2fc > 0, "gene_name"])
airecurrences <- sort(airecurrences[airecurrences > 1], decreasing = T)
airecurrences

aioccurrencespergene <- do.call(rbind, by(data = airesults, INDICES = airesults$gene_name, FUN = function(x) data.frame(samples = paste0(x$sample, collapse = ","),
                                                                                                                n_up = sum(x$log2fc > 0), n_down = sum(x$log2fc < 0), 
                                                                                                                updown = paste0(ifelse(x$log2fc > 0, "+", "-"), collapse = ","),
                                                                                                                gene_name = x$gene_name[1])))
write.table(x = aioccurrencespergene[, c("gene_name", "samples", "updown", "n_up", "n_down")],
            file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20180525_alloccurences_vst.txt", quote = F, sep = "\t", row.names = F)



####

snvcallfiles <- list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/data/Exome_Data/", pattern = "*.raw.snps.indels_VARIANT_FILTRATION.vcf.avinput.hg19_multianno.txt", full.names = T, recursive = T)
snvcallfiles


getXvariants <- function(snvfiles, WESID, gene_name) {
  snvcallfile <- grep(pattern = paste0("-WES_", WESID, "."), x = snvfiles, fixed = T, value = T)
  snvcalls <- read.delim(file = snvcallfile, as.is = T, skip = 1)
  colnames(snvcalls) <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene",
                          "genomicSuperDups", "1000g2012apr_all", "esp6500si_all", "esp6500si_aa", "esp6500_ea", "snp129", "snp137", "cg69", "ExAC_Freq",
                          "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS", "LJB_PhyloP", "LJB_PhyloP_Pred", "LJB_SIFT",
                          "LJB_SIFT_Pred", "LJB_PolyPhen2", "LJB_PolyPhen2_Pred", "LJB_LRT", "LJB_LRT_Pred", "LJB_MutationTaster", "LJB_MutationTaster_Pred",
                          "LJB_GERP++", "clinvar_20140211", "caddgt20", "DMN_Exomes_Freq",
                          "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GT")
  
  fat1calls <- snvcalls[snvcalls$ExAC_Freq == "NaN" & snvcalls$'1000g2012apr_all' == "NaN" & snvcalls$snp137 == "NaN" & snvcalls$FILTER == "PASS", ]
  # fat1calls <- snvcalls
  if (nrow(fat1calls) > 0)
    fat1calls$sample <- paste0("WES_", WESID)
  return(fat1calls)
}

allXcalls <- do.call(rbind, lapply(X = 1:18, FUN = getXvariants, snvfile = snvcallfiles, gene_name = "BCLAF1"))
allfatsnvs <- allfatcalls[allfatcalls$snp137 == "NaN", ]

