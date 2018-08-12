## check recurrence of allelically imbalanced + up/downregulated genes

resultfiles <- c(list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/AML_ase/", pattern = "*imbalance_expression_vst.txt", full.names = T, recursive = T))#,
                 # list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/cell_lines/", pattern = "*imbalance_expression_vst.txt", full.names = T, recursive = T))
## omit sample 22 (mismatch) and 23 (failed exome) for now
resultfiles <- resultfiles[-c(2,3)]

ai_results <- lapply(X = resultfiles, FUN = read.delim, as.is = T)
# View(ai_results[[1]])

aidf <- do.call(rbind, ai_results)
aidf$sample <- rep(sub(pattern = "_imbalance_expression_vst.txt", replacement = "", x = basename(resultfiles)), sapply(X = ai_results, nrow))

aidf <- aidf[aidf$padj <= .05 & (aidf$log2fc >= 1 | aidf$log2fc <= -0.73), ]
aidf <- aidf[aidf$contig != "X" & !grepl(aidf$gene_name, pattern = "HLA-*") &
               !grepl(aidf$gene_name, pattern = "^IG[HLK].*", perl = T) &
               !grepl(aidf$gene_name, pattern = "^TR[ABDG][VCDJ].*", perl = T), ]

recurrent_genes <- sort(table(aidf$gene_name), decreasing = T)
head(recurrent_genes, n = 25)

write.table(x = aidf, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/20180706_allelic_imbalance_pooledsamples_vst_AML.txt", sep = "\t", quote = F, row.names = F, col.names = T)
