#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

##########################################################################

## CHECK SAMPLES CONTAMINATION
# Version for pipe v1.16

##########################################################################


## FUNCTIONS

check_dir = function(output_dir=NULL){
  ## Check if input directory ends with '/' character
  
  if (output_dir[length(output_dir)] != "/"){
    output_dir = paste0(output_dir, "/")
  }
  return(output_dir)
}

## MAIN PROGRAM

dirData = args[1]
dirData = check_dir(dirData)
outdir = paste0(dirData, "QC/mapping/")

dirData = paste0(dirData, 'variant_calling/')

# Read input CSV
csv_file = list.files(dirData, pattern="*annotated.csv")
csv = read.csv(paste0(dirData, csv_file), header = TRUE)

samples2 <- csv[,grep("_FREQUENCY$", colnames(csv))]
samples2[samples2 == "-"] <- 0
samples2[samples2 != "0"] <- 1

matrix = samples2


lista_samples = list()
for (i in 1:ncol(matrix)){
  lista_samples[[i]] = matrix[,i]
}

perc_matrix = as.data.frame(matrix(0, ncol = ncol(samples2), nrow = ncol(samples2)))
colnames(perc_matrix) = sapply(strsplit(colnames(matrix), "_"), `[`, 1)
rownames(perc_matrix) = sapply(strsplit(colnames(matrix), "_"), `[`, 1)
lista2 = lista_samples
names(lista2) = colnames(matrix)
names(lista_samples) = colnames(matrix)
l = 0
n = 0
for (i in 1:(length(lista_samples))){
  n = l + 1
  l = l + 1
  
  for (j in 1:length(lista2)){
    m = 0
    max1 = length(which(lista_samples[[i]] == 1))
    max2 = length(which(lista2[[j]] == 1))
    name_max = names(lista_samples)[i]
    name_min = names(lista2)[j]
    max = lista_samples[[i]]
    min = lista2[[j]]
    if (max2 > max1){
      max = lista2[[j]]
      min = lista_samples[[i]]
      name_min = names(lista_samples)[i]
      name_max = names(lista2)[j]
    }
    for (k in 1:nrow(matrix)){
      if (max[k] == 1){
        if (lista_samples[[i]][k] == lista2[[j]][k]){
          m = m + 1
        }
      }
    }
    
    perc = m/length(which(max == 1))
    perc2 = m/length(which(min == 1))
    perc = round(perc, digits = 4)
    perc2 = round(perc2, digits = 4)
    
    if (is.na(perc)){
      perc = 0
    }
    if (is.na(perc2)){
      perc2 = 0
    }
    
    perc_matrix[n, l] = perc
    perc_matrix[l, n] = perc2
    
    n = n + 1

  }
  lista2[[1]] = NULL
  j = j - 1
}

name = sapply(strsplit(dirData, "/"), `[`, 4)
# Nextflow: nombre de la corrida como segundo argumento (en rutina se deduce de la ruta)
if (length(args) >= 2) name = args[2]
write.csv(perc_matrix, file = paste0(outdir, name, "_sample_contamination.csv"), row.names=TRUE)

