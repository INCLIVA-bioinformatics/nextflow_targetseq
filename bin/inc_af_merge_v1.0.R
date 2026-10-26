#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

##########################################################################

## INC AF

##########################################################################


## FUNCTIONS

split_vcf = function(vcf){
  lista <-  strsplit( as.character( vcf$INFO ) , split=";")
  for( for.fila in 1:nrow (vcf) ){
    for ( for.columna in 1:length( lista[[ for.fila ]] ) ){
      vcf[for.fila, paste0( unlist( strsplit(x = lista[[ for.fila]][ for.columna ] ,
                                             split = "=") )[ 1 ] , ".1" )] <- unlist( strsplit(x = lista[[ for.fila ]][ for.columna ] ,
                                                                                               split = "=") )[ 2 ]
    }
  }
  
  v10<-as.character(vcf[[10]])
  v10lista <-  strsplit(v10, split=":")
  v9<-as.character(vcf$FORMAT)
  v9lista<-strsplit(v9, split=":")
  for( for.fila in 1:nrow ( vcf ) ){
    for ( for.columna in 1:length( v9lista[[ for.fila ]] ) ){
      
      vcf[for.fila, paste0( v9lista [[for.fila]] [for.columna] , ".2" )
          ] <- v10lista [[for.fila]] [for.columna]
    }
  }
  
  # Se eliminan las columnas que hemos separado
  vcf <- vcf[, -c(8, 9, 10, 15)]
  
  return(vcf)
}

den_plot = function(path, title, xlabel, df) {
  
  png(path)
  p = ggplot(df, aes(x=V3, color=Sample)) +
    geom_density(linewidth = 0.8) + ggtitle(title) +
    theme_bw() + ylim(0, 2.5) + xlab(xlabel) + xlim(0,1) +
    theme(plot.title = element_text(size = 15, face = "bold"),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          axis.text=element_text(size=15),
          axis.title=element_text(size=14,face="bold"))
  print (p)
  dev.off()
  
}


## MAIN PROGRAM

library(data.table)
library(ggplot2)

vcf_path = args[1]
outdir = args[2]

#vcf_path = "/media/scratch4/UMP17_run209/tmp_100x/variant_calling"
#outdir = "/media/scratch4/UMP17_run209/analysis_100x/QC/mapping"

stats_path = list.files(outdir, pattern="*stats.csv")

stats = read.table(paste0(outdir, "/", stats_path), sep = ",", header = TRUE, stringsAsFactors = FALSE)
samples = stats$Sample

vcf_plot = as.data.frame(matrix(0, ncol = 2, nrow = 0))

i = 0
j = 0
c = 0
for (sample in samples){
  vcf = fread(paste0("zcat ", paste0(vcf_path,  "/", sample, ".vcf.gz")), skip = "CHROM")
  colnames(vcf)[1] = "CHROM"
  
  if (nrow(vcf) != 0){
    vcf = split_vcf(vcf)
    
    colnames(vcf)[which( colnames(vcf)=="AF.2" )] = "V3"
    vcf$V3 = as.numeric(vcf$V3)
    
    vcf_plot = rbind(vcf_plot, data.frame(rep(sample, each = nrow(vcf)), vcf$V3))
    i = i + 1
  }
  else {
    c = c + 1
  }
  
  if (nrow(vcf) != 0){
    if (i %% 6 == 0 || i == (length(samples)-c)){
      j = j + 1
      
      name = sapply(strsplit(vcf_path, "/"), `[`)
      name = name[4]
      # Nextflow: nombre de la corrida como tercer argumento (en rutina se deduce de la ruta)
      if (length(args) >= 3) name = args[3]
      colnames(vcf_plot) = c("Sample", "V3")
      den_plot(paste0(outdir, "/", name, "_af_plot_", as.character(j), ".png"), paste0(name, " AF Density Plot"), "AF", vcf_plot)
      vcf_plot = as.data.frame(matrix(0, ncol = 2, nrow = 0))
    }
  }
}

