#################################################################
##
##
##  input:
##
##            \HumanEvo\amitai\Matlab_Data\Matlab_data_28.2.19\chr_divided_data
##
##
##  output:   \HumanEvo\HumanEvo\amitai\01.sample_clustering
##
##
##  v_01 14.03.2023
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
INPUT_FOLDER <- "methylation_data"
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


meta_from_Liran <- as.data.frame(readxl::read_xlsx(path = "Meth samples - human.xlsx"))

meta_from_Reich <- as.data.frame(readxl::read_xlsx(path = "AGDP.metadata.xlsx"))

print("1. High coverage published ")
high_coverage_published <- meta_from_Reich[meta_from_Reich$Coverage > 17 & meta_from_Reich$`Year shotgun published` != "Unpublished (this is a prepublication release)" ,]
print(nrow(high_coverage_published))
print("2. High coverage unpublished ")
high_coverage_unpublished <-meta_from_Reich[meta_from_Reich$Coverage > 17 & meta_from_Reich$`Year shotgun published` == "Unpublished (this is a prepublication release)" ,]
print(nrow(high_coverage_unpublished))

printout <- high_coverage_unpublished[,1:3]
colnames(printout) <- c("I-ID", "Coverage", "age")

high_coverage_unpublished_available_data <- high_coverage_unpublished[high_coverage_unpublished$`I-ID` %in% list.files(INPUT_FOLDER),]
print("high_coverage_unpublished_available_data")
print(nrow(high_coverage_unpublished_available_data))

high_coverage_unpublished_missing_data <- high_coverage_unpublished[!high_coverage_unpublished$`I-ID` %in% list.files(INPUT_FOLDER),]
print("high_coverage_unpublished_missing_data")
print(nrow(high_coverage_unpublished_missing_data))

printout <- high_coverage_unpublished_missing_data[,1:3]
colnames(printout) <- c("I-ID", "Coverage", "age")
printout

#nrow(high_coverage_published$`I-ID` %in% meta_from_Liran$Sample , ])

reich_unpublished_not_in_liran <- meta_from_Reich[meta_from_Reich$Coverage > 17 & 
                                  !meta_from_Reich$`I-ID` %in% meta_from_Liran$Sample & 
                                  meta_from_Reich$`Year shotgun published` == "Unpublished (this is a prepublication release)", ]
print("reich_samples not in lirans")
print(paste( nrow(reich_unpublished_not_in_liran)))

meta_from_Reich[meta_from_Reich$Coverage > 17 & 
                  !meta_from_Reich$`I-ID` %in% meta_from_Liran$Sample & 
                  meta_from_Reich$`Year shotgun published` == "Unpublished (this is a prepublication release)", 
                c("I-ID","Coverage","Date mean in BP [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]")]

nrow(meta_from_Reich[meta_from_Reich$Coverage > 17 & meta_from_Reich$`I-ID` %in% meta_from_Liran$Sample, ])
nrow(meta_from_Liran[meta_from_Liran$Sample %in% meta_from_Reich$`I-ID`, ])

nrow(meta_from_Liran[!meta_from_Liran$Sample %in% meta_from_Reich$`I-ID` & meta_from_Liran$Coverage > 17, ])

nrow(meta_from_Liran[!meta_from_Liran$Sample %in% meta_from_Reich$`I-ID` & 
                       meta_from_Liran$Coverage > 17 & 
                       meta_from_Liran$`Lab sequenced`== "Reich",])


meta_from_Liran[!meta_from_Liran$Sample %in% meta_from_Reich$`I-ID` & 
                  meta_from_Liran$Coverage > 17 & 
                  meta_from_Liran$`Lab sequenced`== "Reich",c("Sample","Coverage","age of sample (kyo)","tissue")]
##############################################################


sum(meta_from_Reich$`I-ID` %in% list.files(INPUT_FOLDER))

