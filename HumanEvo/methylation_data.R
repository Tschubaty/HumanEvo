# Documentation -----------------------------------------------------------
##
##
##  input: "methylation_data"
##
##
##  output: "methylation_data"
##
##
##  v_01 20.03.2024
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------

#Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
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
library(data.table)

# Constants defined previously
INPUT_FOLDER <- "methylation_data"
OUTPUT_FOLDER <- "methylation_data"
CHR_NAMES <- c(
  "chr1", "chr2", "chr3", "chr4", "chr5",
  "chr6", "chr7", "chr8", "chr9", "chr10",
  "chr11", "chr12", "chr13", "chr14", "chr15",
  "chr16", "chr17", "chr18", "chr19", "chr20",
  "chr21", "chr22"
)

convert_experiment_bed_to_wig <- function(exp_dir, output_wig_file) {
  # Exclude .wig files when listing .bed files
  bed_files <- list.files(exp_dir, pattern = "\\.bed$", full.names = TRUE)
  
  wig_conn <- file(output_wig_file, "w")
  writeLines("track type=wiggle_0 name=\"Converted from BED\" description=\"Variable Step Format\"", wig_conn)
  
  all_data <- rbindlist(lapply(bed_files, function(bed_file) {
    dt <- fread(bed_file, select = c("chrom", "start", "score"))
    dt[, .(chrom, position = as.numeric(start) + 1, score)]
  }))
  
  for (chr in CHR_NAMES) {
    if (chr %in% all_data$chrom) {
      chr_data <- all_data[chrom == chr]
      writeLines(paste0("variableStep chrom=", chr, " span=1"), wig_conn)
      for (i in 1:nrow(chr_data)) {
        writeLines(paste(chr_data$position[i], chr_data$score[i]), wig_conn)
      }
    }
  }
  
  close(wig_conn)
}

experiment_dirs <- list.files(INPUT_FOLDER, pattern = "^[^\\.].*", full.names = TRUE) # Ignore hidden files/directories

for (exp_dir in experiment_dirs) {
  exp_id <- basename(exp_dir)
  # Skip processing if the directory name ends with '.wig' to avoid processing output from previous runs
  if(!grepl("\\.wig$", exp_id)) {
    output_wig_file <- file.path(OUTPUT_FOLDER, paste0(exp_id, ".wig"))
    convert_experiment_bed_to_wig(exp_dir, output_wig_file)
    cat("Converted BED files to WIG for experiment:", exp_id, "\n")
  }
}


# Read metadata from Excel file
meta <- as.data.frame(read_xlsx("AGDP.39.metadata2.xlsx"))

# Specify source and destination directories
source_dir <- INPUT_FOLDER
dest_dir <- "methylation_wig"

# Loop over each file name and move/rename the file
for (s in meta$sample) {
  print(s)
  # Construct source and destination paths
  source_path <- file.path(source_dir, paste(s,"wig",sep = "."))
  dest_path <- file.path(dest_dir, paste(s,"hg38","wig",sep = "."))
  
  # Move and rename the file
  file.rename(source_path, dest_path)
}
