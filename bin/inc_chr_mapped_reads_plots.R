#!/usr/bin/env Rscript

# Script para generar plots de las reads mapeadas en los cromosomas y de las proporcion chrY/chrX
# Con genomas de referencia que usen los nombres de cromosoma definidos abajo dara error

args = commandArgs(trailingOnly=TRUE)


file <- args[1] #"/media/scratch0/TEST_ONC4/tmp_ONC4_100x/mapping/ONC4-X2301155_idxstats"
outdir <- args[2] #"~/pruebasJF/pruebas_sexo_ONC4"

library(ggplot2)
library(dplyr)

sample_name <- strsplit(file, "/")
sample_name <- gsub("_idxstats", "", sample_name[[1]][[length(sample_name[[1]])]])

# Cargamos el fichero idxstats
idx <- read.table(file = file, header = F, sep = "\t")
colnames(idx) <- c("seqnames", "length", "mapped", "unmapped") 
# Faltan ordenar los cromosomas


  
plot <- idx %>% mutate(seqnames = factor(seqnames, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                                                              "chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                                                              "chr17","chr18","chr19","chr20","chr21","chr22","chrX",
                                                              "chrY","chrMT","*" ))) %>% 
  ggplot(aes(x = seqnames, y = mapped)) + 
  geom_bar(stat = "identity") + 
  xlab("Cromosome") + ylab("Mapped reads") +
  ggtitle(paste0("Mapped reads per chromosome in sample ",sample_name)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = paste0(outdir,"/",sample_name,"_mapped_by_chromosome.png"), plot = plot)

plot2 <- idx %>% dplyr::filter(seqnames %in% c("chrX", "chrY")) %>% mutate(Sample = sample_name) %>% 
  mutate(Chromosome = seqnames) %>%
  ggplot(aes(x = Sample, y = mapped, fill = Chromosome)) +
  geom_bar(position="fill", stat = "identity") +
  xlab("Sample") + ylab("Proportion of mapped reads") +
  scale_fill_manual(values=c("#1f77b4", "#ff7f0e")) +
  coord_flip()
ggsave(filename = paste0(outdir,"/",sample_name,"_chrY_chrX.png"), plot = plot2)

