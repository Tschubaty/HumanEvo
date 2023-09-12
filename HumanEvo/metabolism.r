#################################################################
##
##
##  input:"05.annotation"
##        
##  output: "metabolsim"    
##
##
##  v_01 21.04.2023
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
#library(Rtsne)
#library(ggrepel)
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- "metabolsim"
INPUT_FOLDER <- "05.annotation"
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

min_delta <- 0.25

meta <- readRDS(file = file.path("03.plots",paste("meta","rds",sep = ".")))
sample_order <- as.character(meta$sample[order(meta$age_mean_BP)])
age <- meta$age_mean_BP[order(meta$age_mean_BP)]

alpha <- 0.001

######################## metabolsim

df <- as.data.frame(readxl::read_xlsx(path = file.path(INPUT_FOLDER,"delta_0.25.annotaed_peak_with_best_CpG_results.p.0.001.xlsx")))

# https://esbl.nhlbi.nih.gov/Databases/KSBP2/Targets/Lists/MetabolicEnzymes/MetabolicEnzymeDatabase.html
df_metabolism <- as.data.frame(readxl::read_xlsx(path = file.path(OUTPUT_FOLDER,"Mammalian_Metabolic_Final.xlsx")))

intersect(toupper(df_metabolism$`Gene Symbol`),toupper(df$annotation_narrow))

write.table(
  x = df_metabolism[toupper(df_metabolism$`Gene Symbol`) %in% toupper(df$annotation_narrow), ],
  file = file.path(OUTPUT_FOLDER,"results_in_MetabolicEnzymeDatabase.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


