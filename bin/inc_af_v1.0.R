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

den_plot = function(path, title, xlabel, df, min, max, min_lim, max_lim) {
  
  png(path)
  p = ggplot(df, aes(x=V3)) +
    geom_density(size = 0.8) + ggtitle(title) +
    theme_bw() + xlim(as.numeric(min) - min_lim, as.numeric(max) + max_lim) + xlab(xlabel) +
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

vcf = fread(paste0("zcat ", vcf_path), skip = "CHROM")
colnames(vcf)[1] = "CHROM"

if (nrow(vcf) != 0){
  vcf = split_vcf(vcf)
  
  colnames(vcf)[which(colnames(vcf) == "AF.2")] = "V3"
  vcf$V3 = as.numeric(vcf$V3)
  
  min = 0
  max = 1
  
  name = sapply(strsplit(vcf_path, "/"), `[`)
  name = name[length(name)]
  name = sapply(strsplit(name, "[.]"), `[`)[1]
  
  den_plot(paste0(outdir, "/", name, "_af_plot.png"), paste0(name, " AF Density Plot"), "AF", vcf, min, max, 0, 0)
}

