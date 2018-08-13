## process TCGA LAML gdc_smaple_table
library(rjson)
tcgasamplesheet <- read.delim("/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/gdc_sample_sheet.2018-08-12.tsv", as.is = T)
tcgajson <- fromJSON(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/files.2018-08-12.json")
# tcgaRNA <- tcgasamples[sapply(tcgasamples, function(x) x$experimental_strategy == "RNA-Seq")]

tcgajsondf <- as.data.frame(do.call(rbind, lapply(tcgajson, unlist)))

colnames(tcgajsondf)
colnames(tcgasamplesheet)

tcgadata <- merge(x = tcgasamplesheet, y = tcgajsondf, by.x = "File.Name", by.y = "file_name")

# omit normal ones and pick first for each
tcgadata <- tcgadata[tcgadata$Sample.Type == "Primary Blood Derived Cancer - Peripheral Blood", ]
tcgadata <- tcgadata[!duplicated(tcgadata[, c("Sample.ID", "experimental_strategy")]), ]
tcgadata <- tcgadata[tcgadata$Sample.ID %in% names(which(table(tcgadata$Sample.ID) == 2)), ]

tcgaWXS <- tcgadata[tcgadata$experimental_strategy == "WXS", ]
tcgaRNA <- tcgadata[tcgadata$experimental_strategy == "RNA-Seq", ]
tcgadata <- merge(x = tcgaWXS, y = tcgaRNA, by = "Sample.ID")
tcgadata <- tcgadata[, c("Case.ID.x", "Sample.ID", "File.Name.x", "File.Name.y")]
colnames(tcgadata) <- c("case_id", "sample_id", "wxs_bam", "rna_bam")

write.table(file = "/srv/shared/vanloo/home/jdemeul/projects/2016_mansour_ASE_T-ALL/results/20180812_tcgadata_formatted.tsv", x = tcgadata, quote = F, sep = "\t")
