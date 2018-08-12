## QC on mapping bias

library(ggplot2)

confhetsnpsfile <- file.path(MOUTDIR, paste0(MSAMPLE, "_hetSNPs.txt"))
confhetsnps <- read.delim(file = confhetsnpsfile, as.is = T)
confhetsnps$baf <- confhetsnps$count_alt / (confhetsnps$count_alt + confhetsnps$count_ref)
q1 <- ggplot(data = confhetsnps, mapping = aes(x = baf)) + geom_histogram(binwidth = 0.01) + scale_x_continuous(limits = c(0,1))
q1

q1 <- ggplot(data = confhetsnps, mapping = aes(x = baf, y = count_ref + count_alt)) + geom_jitter(alpha = .3, shape = ".") + scale_x_continuous(limits = c(0,1))
q1

summary(confhetsnps$baf)
