#################################################################
##
##
##  input:
##
##   choose_best_P_val_per_peak
##
##  output:
##
##
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
library(foreach)
library(doParallel)
library(ggplot2)
#library(Rtsne)
#library(ggrepel)
library("parallel")
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- "gatherData"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
ANNOTATION <- "hg19"
N_SAMPLES <- 37
CUTVALUE <- 300

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

meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
real_sample_typisation <-
  readRDS(file.path("04.methylation_vs_Type", "real_sample_typisation.rds"))
meta37 <- meta[meta$sample %in% real_sample_typisation$sample,]

peak_CpG <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/ChipSegmentCut/all_unfiltered_peak_CpG_hg19_cut_peak_annotation.rds")

# count number of datapoints
peak_CpG$score <-
  apply(
    X = peak_CpG[,meta$sample],
    MARGIN = 1,
    FUN = function(r) {
      return(sum(!is.na(r)))
    }
  )

saveRDS(object = peak_CpG,file = file.path(OUTPUT_FOLDER,"all_unfiltered_peak_CpG_hg19_cut_peak_annotation.rds"))
print(nrow(peak_CpG))
saveRDS(object = real_sample_typisation,file = file.path(OUTPUT_FOLDER,"37real_sample_typisation.rds"))
saveRDS(object = meta37,file = file.path(OUTPUT_FOLDER,"37meta.rds"))
saveRDS(object = meta,file = file.path(OUTPUT_FOLDER,"meta.rds"))

peak_37CpG <- peak_CpG[, c("chrom","start","end","name","score","strand", real_sample_typisation$sample)]

saveRDS(object = peak_37CpG[,real_sample_typisation$sample],file = file.path(OUTPUT_FOLDER,"37_meth_values_unfiltered_peak_CpG_hg19_cut_peak_annotation.rds"))


# count number of datapoints
peak_37CpG$score <-
  apply(
    X = peak_37CpG[,real_sample_typisation$sample],
    MARGIN = 1,
    FUN = function(r) {
      return(sum(!is.na(r)))
    }
  )

saveRDS(object = peak_37CpG,file = file.path(OUTPUT_FOLDER,"37_unfiltered_peak_CpG_hg19_cut_peak_annotation.rds"))

saveRDS(object = peak_37CpG[N_SAMPLES == peak_37CpG$score, ],file = file.path(OUTPUT_FOLDER,"37_completecase_peak_CpG_hg19_cut_peak_annotation.rds"))

saveRDS(object = peak_37CpG[N_SAMPLES == peak_37CpG$score,real_sample_typisation$sample],file = file.path(OUTPUT_FOLDER,"37_meth_values_completecase_peak_CpG_hg19_cut_peak_annotation.rds"))



