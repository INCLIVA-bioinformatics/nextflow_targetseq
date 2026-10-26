#!/usr/bin/env Rscript

# inc_split_annotated_csv_v1.0.R

# Este script se crea a partir de inc_csv_format_v1.1.R. Ese script y inc_csv_format_conf_v1.1.R eran iguales, por lo
# que se crea esta versión para unificar esa parte del pipeline. Además, el nombre de este archivo es mucho más
# descriptivo respecto a su función.

# Versión: 1.0 (created for pipe v1.15)
# Autor: UB INCLIVA - Basado en inc_csv_format_v1.1.R
################################################################################
library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("-c", "--csv", dest = "input_csv", type = "character", help = "Ruta al archivo csv de entrada.")
parser$add_argument("-o", "--output_directory", dest = "output_directory", type = "character", help = "Ruta al directorio donde se generarán los archivos de salida.")

# Parsear los argumentos de línea de comandos
args <- parser$parse_args()

## FUNCTIONS

check_dir = function(output_dir=NULL){
    ## Check if input directory ends with '/' character

    if (output_dir[length(output_dir)] != "/"){
        output_dir = paste0(output_dir, "/")
        }
    return(output_dir)
}

select_sample_columns_and_save <- function(csv, sample){
  sample_columns <- colnames(csv)[str_detect(colnames(csv),
                                             pattern = str_replace_all(sample,
                                                                       pattern = "\\.",
                                                                       replacement = "\\\\."),
                                             negate = F)]
  sample_columns <- c(general_columns, sample_columns)
  df <- csv[,sample_columns]
  
  # Filtrar filas donde muestra_DP y muestra_FREQUENCY no son "-"
  dp_column <- paste0(sample, "_DP")
  frequency_column <- paste0(sample, "_FREQUENCY")
  
  df <- df %>% filter(!!sym(dp_column) != "-", !!sym(frequency_column) != "-")
  
  write.csv(df, file = paste0(output_dir, sample, ".csv"), row.names = FALSE)
}

## MAIN PROGRAM

# User ask for csv and outdir path
csv_path = args$input_csv
output_dir = args$output_directory

output_dir = check_dir(output_dir)

# Read input CSV
csv = read.csv(csv_path, header = TRUE, colClasses = c("character"))

sample_names <- colnames(csv)[str_detect(colnames(csv), pattern = "_FREQUENCY")] %>%
  str_remove(pattern = "_FREQUENCY")

csv <- csv %>%
  separate(col = "FILTERMUTECTCALLS", sep = "\\|", into = paste(sample_names,"FILTERMUTECTCALLS",sep = "_")) %>%
  separate(col = "STRAND_BIAS_BY_CALLER", sep = "\\|", into = paste(sample_names,"STRAND_BIAS_BY_CALLER",sep = "_"))

str <- str_replace_all(paste(sample_names, collapse = "|"), pattern = "\\.", replacement = "\\\\.")
general_columns <- colnames(csv)[str_detect(colnames(csv), pattern = str, negate = T)]

lapply(sample_names, function(x){select_sample_columns_and_save(csv,x)})
