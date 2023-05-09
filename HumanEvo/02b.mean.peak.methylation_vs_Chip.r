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
#chr <- CHR_NAMES[22]

df_merge_chip <- data.frame()

for(chr in CHR_NAMES){
  
  print(chr)
  # only one chromosome
  df_chr_chip <- df_H3K27ac[df_H3K27ac$chr == chr,]
  # sort dataframe by coodinates
  df_chr_chip <- df_chr_chip[order(df_chr_chip$start),]
  # set name and score variable to NaN
  df_chr_chip$name <- "H3K27ac"
  # number of CpG in peak bounderies
  df_chr_chip$score <- NaN
  # load mthylation for chromosome
  df_chr_meth <- readRDS(file = file.path(INPUT_METH,paste(chr,N_SAMPLES,"samples","rds",sep = ".")))
  df_chr_chip[,c(colnames(df_chr_meth)[7:ncol(df_chr_meth)])] <- NaN
  
  for(r in 1:nrow(df_chr_chip)){
    p_start <- df_chr_chip$start[r]
    p_end <- df_chr_chip$end[r]
    CpGs <- df_chr_meth[p_start <= df_chr_meth$start & df_chr_meth$end <= p_end,] 
    df_chr_chip$score[r]<- nrow(CpGs)
    df_chr_chip[r,c(colnames(CpGs)[7:ncol(CpGs)])] <- apply(X = CpGs[,c(colnames(CpGs)[7:ncol(CpGs)])],
                                                            MARGIN = 2,
                                                            FUN = function(input) {return(mean(x = input,,na.rm = TRUE))}
                                                              )
  }
# merge
  df_merge_chip <- rbind(df_merge_chip,df_chr_chip)
}
# save merged mean peak mewthylation H3K27ac
saveRDS(object = df_merge_chip ,file = file.path(OUTPUT_FOLDER,paste("H3K27ac","peak","mean","methylation",N_SAMPLES,"samples","merged",ANNOTATION,"rds",sep = ".")))