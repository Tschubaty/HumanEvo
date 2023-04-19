#################################################################
##
##
##  input:
##           "methylation+chip"
##            
##
##
##  output:   
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
# library(foreach)
# library(doParallel)
library(ggplot2)
library(Rtsne)
library(ggrepel)
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- 
INPUT_FOLDER <- "methylation+chip"
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


df <- readRDS(file = file.path(INPUT_FOLDER,paste("H3K27ac","only",N_SAMPLES,"samples","merged",ANNOTATION,"rds",sep = ".")))

df$score <- apply(df[,c(7:ncol(df))], 1 ,FUN = function(r) {sum(!is.na(r))})

ggplot(data = df,mapping =  aes(score) )+geom_histogram(bins = N_SAMPLES+1)+theme_minimal()+ggtitle("# samples with values for a CpG")

# NO DATA for I2978
#[1] TRUE TRUE

df_complete <- df[df$score == N_SAMPLES,]

