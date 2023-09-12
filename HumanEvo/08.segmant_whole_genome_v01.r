#################################################################
##
##
##  input:
##           "methylation+chip"
##
##
##
##  output: methylation+chip
##
##
##  v_01 21.08.2023
##  Author: Daniel Batyrev 777634015
#################################################################
#Clear R working environment
rm(list = ls())
cluster <- TRUE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/"
  picuture_file_extension <- "pdf"
} else{
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "png"
}
setwd(this.dir)
detachAllPackages <- function() {
  basic.packages <-
    c(
      "package:stats",
      "package:graphics",
      "package:grDevices",
      "package:utils",
      "package:datasets",
      "package:methods",
      "package:base"
    )
  package.list <-
    search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]
  package.list <- setdiff(package.list, basic.packages)
  if (length(package.list) > 0)
    for (package in package.list)
      detach(package, character.only = TRUE)
}
detachAllPackages()
# #################################### Libs ########################################
library(foreach)
library(doParallel)
library(ggplot2)
#library(Rtsne)
#library(ggrepel)
##################################### CONSTANTS ########################################
INPUT_FOLDER <- "07.methylation_p_val"
OUTPUT_FOLDER <- "08.genome_segemntation"
N_SAMPLES <- 39
ANNOTATION <- "hg19"

CHR_NAMES <-
  c(
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22"
  )

##################################################

window_width <- 400
min_delta <- 0.25
offset = 0

# paralle processing reigstration
numCores <- detectCores() - 1
registerDoParallel(numCores)

# for each chromosome
for (chr in rev(CHR_NAMES)) {
  #chr <- CHR_NAMES[22]
  
  start_time <- Sys.time()

  print(chr)
  file_name <-
    paste(chr,
          ANNOTATION,
          N_SAMPLES,
          "samples",
          "p_val",
          "delta_meth",
          "rds",
          sep = ".")
  
  df_chr <- readRDS(file = file.path(INPUT_FOLDER, file_name))
  # number of windows
  #n_windows <- floor(df_chr$end[nrwo(df_chr)] - df_chr$start[1]) / window_width)
  start_coordinates <-
    seq(from =  df_chr$start[1],
        to = df_chr$end[nrow(df_chr)],
        by = window_width)
  
  df_results <-
    foreach(st = start_coordinates, .combine = rbind) %dopar% {
      #for (st in start_coordinates[1:10]) {
      # st <- start_coordinates[100]
      indices <-
        (st <= df_chr$start  & df_chr$end <= st + window_width)
      n_CpG <- sum(indices)
      if (n_CpG  > 0) {
        delta_meth <- df_chr$delta_meth[indices]
        p_val <- df_chr$p_val[indices]
        # relevant CpGs
        index_CpG_relevant <- delta_meth >= min_delta & !is.na(p_val)
        
        if (sum(index_CpG_relevant) > 0) {
          start <- df_chr$start[indices]
          end <- df_chr$end[indices]
          name <- df_chr$name[indices]
          score <- df_chr$score[indices]
          pearson_cor <- df_chr$pearson_cor[indices]
          
          min_p <- min(p_val[index_CpG_relevant])
          min_index <- p_val == min_p
          
          return(
            data.frame(
              chr = chr,
              start = st,
              end = st + window_width,
              name = paste(window_width, "bp", sep = ""),
              score = n_CpG,
              strand = ".",
              best_start = start[min_index],
              best_end = end[min_index],
              best_name = name[min_index],
              best_n_sample = score[min_index],
              best_pearson_cor = pearson_cor[min_index],
              best_delta_meth = delta_meth[min_index],
              p_val = min_p
            )
          )
        }
      }
      return(
        data.frame(
          chr = chr,
          start = st,
          end = st + window_width,
          name = paste(window_width, "bp", sep = ""),
          score = n_CpG,
          strand = ".",
          best_start = NA,
          best_end = NA,
          best_name = NA,
          best_n_sample = NA,
          best_pearson_cor = NA,
          best_delta_meth = NA,
          p_val = NA
        )
      )
    }
  saveRDS(object = df_results,
          file = file.path(
            OUTPUT_FOLDER,
            paste(
              chr,
              as.character(window_width),
              "bp",
              as.character(offset),
              "offset",
              min_delta,
              "min_delta",
              "pearson",
              "rds",
              sep = "."
            )
          ))
  end_time <- Sys.time()
  print(end_time - start_time)
  gc()
}
stopImplicitCluster()
gc()