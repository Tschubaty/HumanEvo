#################################################################
##
##
##  input:
##            Chromatin\processed
##            methylation_data_chr_merged
##
##
##  output:   methylation+chip
##
##
##  v_01 19.04.2023
##  Author: Daniel Batyrev 777634015
#################################################################
#Clear R working environment
rm(list = ls())
cluster <- FALSE
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
#library(ggplot2)
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- "methylation+chip"
INPUT_METH <- "methylation_data_chr_merged"
INPUT_CHROMATIN <- file.path("Chromatin","processed")
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

################ CODE ###########################

df_H3K27ac <- readRDS(file = file.path(INPUT_CHROMATIN,
                                       paste("H3K27ac",ANNOTATION,"RDS",sep = ".")))
df_merged <- data.frame()
df_full_merged <- data.frame()
#chr <- CHR_NAMES[22]

for(chr in CHR_NAMES){
  
  print(chr)
  # only one chromosome
  df_chr_chip <- df_H3K27ac[df_H3K27ac$chr == chr,]
  # sort dataframe by coodinates
  df_chr_chip <- df_chr_chip[order(df_chr_chip$start),]
  # load mthylation for chromosome 
  df_chr_meth <- readRDS(file = file.path(INPUT_METH,paste(chr,N_SAMPLES,"samples","rds",sep = ".")))
  df_chr_meth$name <- "NO_CHIP"
  
  for(r in 1:nrow(df_chr_chip)){
    p_start <- df_chr_chip$start[r]
    p_end <- df_chr_chip$end[r]
    df_chr_meth$name[p_start <= df_chr_meth$start & df_chr_meth$end <= p_end] <- paste(chr,p_start,p_end,"H3K27ac",sep = ".")
  }
  # save chromosome file
  saveRDS(object = df_chr_meth,file = file.path(OUTPUT_FOLDER,paste(chr,ANNOTATION,N_SAMPLES,"samples","rds",sep = ".")))
  df_merged <- rbind(df_merged,df_chr_meth[df_chr_meth$name != "NO_CHIP",])
  df_full_merged <- rbind(df_full_merged,df_chr_meth)
}
# save merged only H3K27ac
saveRDS(object = df_merged ,file = file.path(OUTPUT_FOLDER,paste("H3K27ac","only",N_SAMPLES,"samples","merged",ANNOTATION,"rds",sep = ".")))
# save all
saveRDS(object = df_full_merged ,file = file.path(OUTPUT_FOLDER,paste("all_CpG",N_SAMPLES,"samples","merged",ANNOTATION,"rds",sep = ".")))
