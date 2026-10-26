#!/usr/bin/env Rscript
print("v1")
library(xtable)
args = commandArgs(trailingOnly=TRUE)
tab <- read.csv(args[1], sep="\t", header=F)
capture.output(summary(tab$V3), file = args[2])
pdf(args[3])
plot(density(tab$V3), main = args[4], xlab = "Depth of coverage")
dev.off()
