## create final plot to show recurrences and outlying log2fc

# # get sample IDs mappings
# sampledf <- read.delim(file = "20180626_aml_samples.txt", as.is = T)
# ## omit sample 22 (mismatch) and 23 (failed exome) for now
# sampledf <- sampledf[!sampledf$Number %in% c(22,23), ]
# sampledf$Number <- paste0("X", sampledf$Number)
# sampledf <-  rbind(sampledf, c(Number = "CD34", 'No.' = "", Initials = "", id = "CD34"))


# get log2fc matrix (VST)
l2fcfile <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180706_RNAlog2fc_vst_AML.txt"
l2fcdf <-  read.delim(file = l2fcfile, as.is = T)

# reformat column names
# colnames(l2fcdf)[sapply(X = sub(pattern = "-", replacement = ".", sampledf$sampleid), FUN = grep, x = colnames(l2fcdf), simplify = T)] <- c(paste0("WES_", sampledf[!sampledf$cell_line, "t_wes_id"]), sub(pattern = "-", replacement = ".", x = sampledf[sampledf$cell_line, "sampleid"]))

# get allelic imbalance pooled samples file
airesultsfile <- "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180817_alloccurences_vst_AML.txt"
airesults <- read.delim(file = airesultsfile, sep = "\t", as.is = T)
airesults_sub <- airesults[airesults$n_up >= 2, ]

# splitsamplestring <- strsplit(airesults_sub$samples, split = ",", fixed = T)
# imbalanced_gene_sample_combinations <- paste0(rep(x = airesults_sub$gene_name, lengths(splitsamplestring)), "_", unlist(splitsamplestring))
# 
l2fcdf_sub <- l2fcdf[l2fcdf$gene_name %in% airesults_sub$gene_name, ]
# l2fcdf_sub$exprrank <- (1:nrow(l2fcdf_sub))[order(l2fcdf_sub$mean_expression, decreasing = F)]
l2fcdf_sub$gene_name <- factor(x = l2fcdf_sub$gene_name, levels = l2fcdf_sub[order(l2fcdf_sub$mean_expression, decreasing = F), "gene_name"])
# subset log2fc mat to the represented genes (maybe those that occur at least twice to get overview on X)



### MOD, include ALL AI cases irrespective of log2fc
resultfiles <- c(list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/AML_ase/", pattern = "*imbalance_expression_vst.txt", full.names = T, recursive = T))#,
                 # list.files(path = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/ASE_analysis/cell_lines/", pattern = "*imbalance_expression_vst.txt", full.names = T, recursive = T))

## omit sample 22 (mismatch) and 23 (failed exome) for now
resultfiles <- resultfiles[-c(2)]

airesults_all <- lapply(X = resultfiles, FUN = read.delim, as.is = T)
# View(ai_results[[1]])

airesults_all_df <- do.call(rbind, airesults_all)
airesults_all_df$sample <- sub(pattern = "-", replacement = ".", rep(sub(pattern = "_imbalance_expression_vst.txt", replacement = "", x = basename(resultfiles)), sapply(X = airesults_all, nrow)))

airesults_all_df <- airesults_all_df[airesults_all_df$padj <= .05, ]
airesults_all_df <- airesults_all_df[airesults_all_df$contig != "X" & !grepl(airesults_all_df$gene_name, pattern = "HLA-*") &
                                       !grepl(airesults_all_df$gene_name, pattern = "^IG[HLK].*", perl = T) &
                                       !grepl(airesults_all_df$gene_name, pattern = "^TR[ABDG][VCDJ].*", perl = T), ]

imbalanced_gene_sample_combinations <- paste0(airesults_all_df$gene_name, "_", airesults_all_df$sample)
###



library(reshape2)
l2fcdf_melt <- melt(data = l2fcdf_sub, id.vars = c("gene_name", "mean_expression"), variable.name = "sample_id", value.name = "l2fc")
l2fcdf_melt$sample_id <- sub(pattern = "^X", replacement = "", x = l2fcdf_melt$sample_id)
l2fcdf_melt$is_ai <- paste0(l2fcdf_melt$gene_name, "_", l2fcdf_melt$sample_id) %in% imbalanced_gene_sample_combinations

# plot mean counts (rank axis) vs log2fc_vst for all samples for those genes
# simple dot for each, larger dot for samples in panel and fill for imbalanced ones
library(ggplot2)

p1 <- ggplot(data = l2fcdf_melt, mapping = aes(x = gene_name, y = l2fc))
# p1 <- p1 + geom_point(alpha = .5, shape = 16, size = 1.5, stroke = 0)
p1 <- p1 + geom_violin(scale = "width", alpha =.5)
p1 <- p1 + geom_point(data = l2fcdf_melt[l2fcdf_melt$sample_id %in% c("CD34"), ], colour = "brown", alpha = .5, size = 2, stroke = 0)
p1 <- p1 + geom_point(data = l2fcdf_melt[l2fcdf_melt$sample_id %in% c("HL60", "MOLM16", "OCI_AML2", "OCI_AML3", "TF1", "THP1_S13", "THP1_S6"), ], colour = "green", alpha = .5, size = 2, stroke = 0)
p1 <- p1 + geom_point(data = l2fcdf_melt[!l2fcdf_melt$sample_id %in% c("CD34","HL60", "MOLM16", "OCI_AML2", "OCI_AML3", "TF1", "THP1_S13", "THP1_S6"), ], alpha = .5, size = 2, stroke = 0)
p1 <- p1 + geom_point(data = l2fcdf_melt[l2fcdf_melt$is_ai, ], colour = "red", alpha = .6, shape = 1, size = 2, stroke = 1)
p1 <- p1 + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
p1
ggsave(filename = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180817_l2fc_vst_AIrecurrenceGreaterThan1_AML.png", plot = p1, dpi = 300, width = 15, height = 6)
write.table(x = l2fcdf_melt, file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180817_l2fc_vst_AIrecurrenceGreaterThan1_table_AML.txt", sep = "\t", col.names = T, row.names = F, quote = F)
