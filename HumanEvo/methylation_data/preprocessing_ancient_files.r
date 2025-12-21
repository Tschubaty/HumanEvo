# Documentation -----------------------------------------------------------
##
##  input: ancient files hitogramm matching 
##
##  output: Methylation files in hg19 with + and - strands combined
##
##  v_04 09.04.2025
##  Author: Daniel Batyrev (HUJI 777634015)
##

# Set up Work Environment --------------------------------------------------

# Clear R working environment
rm(list = ls())
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)#
library(data.table)  # Faster data handling

#Set working directory ####################
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
} else {
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
}
setwd(this.dir)

# set refewrnce directory 
reference_dir <- file.path(this.dir,"GSE276666_RAW/methylation_data/hg19reference")

input_dir <- file.path(this.dir, "Ancientmethylation2025")

# Set output directory
output_dir <- file.path(input_dir,"merged")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## only big sample groups
TYPES <- c("Farmer", "Steppe", "HG","Modern")

meta <- as.data.frame(readxl::read_xlsx(file.path(dirname(this.dir),"AGDP.45.metadata.2.1.xlsx")))
meta <- meta[order(meta$age_mean_BP), ]
meta$sample <- factor(x = meta$sample, levels = meta$sample)

# 1. Compute mean age per Type
type_order <- meta %>%
  group_by(Type) %>%
  summarise(mean_age = mean(age_mean_BP, na.rm = TRUE)) %>%
  arrange(mean_age)

# 2. Reorder factor levels by ascending mean age
meta$Type <- factor(meta$Type, levels = type_order$Type)


deaminationRates<- as.data.frame(readxl::read_xlsx(file.path(dirname(this.dir),"deaminationRatesChen2025.xlsx")))

meta <- left_join(meta, deaminationRates, by = "sample") 

# Write to new Excel file
writexl::write_xlsx(meta, file.path(dirname(this.dir),"AGDP.45.metadata.2.2.xlsx"))

library(data.table)

# Define input directory
modern_dir <- file.path(this.dir,"GSE276666_RAW")

# List all coverage files
coverage_files <- list.files(modern_dir, pattern = "^coverage_sorted_S\\d+\\.txt$", full.names = TRUE)

# Initialize result list
coverage_list <- list()

# Loop over files
for (file in coverage_files) {
  sample_id <- sub("^coverage_sorted_(S\\d+)\\.txt$", "\\1", basename(file))
  
  dt <- fread(file, skip = 1, col.names = c(
    "rname", "startpos", "endpos", "numreads", "covbases",
    "coverage", "meandepth", "meanbaseq", "meanmapq"
  ))
  
  # Filter to chr1–22
  dt <- dt[rname %in% as.character(1:22)]
  dt[, length := endpos - startpos + 1]
  
  # Weighted mean coverage
  Coverage <- sum(dt$meandepth * dt$length) / sum(dt$length)
  
  coverage_list[[sample_id]] <- Coverage
}

# Combine into data.frame
coverage_df <- data.frame(
  sample = names(coverage_list),
  Coverage = unlist(coverage_list),
  row.names = NULL
)

print(coverage_df)

library(dplyr)

meta <- left_join(meta, coverage_df, by = "sample") %>%
  mutate(Coverage = ifelse(is.na(Coverage.x), Coverage.y, Coverage.x)) %>%
  select(-Coverage.x, -Coverage.y)


# Write to new Excel file
writexl::write_xlsx(meta, file.path(dirname(this.dir),"AGDP.45.metadata.2.3.xlsx"))

# Load or Create CpG Sites Data --------------------------------------------
# if (FALSE) {
#   # Ensure output directory and subfolder exist
#   hg19_ref_dir <- file.path(output_dir, "hg19reference")
#   if (!dir.exists(hg19_ref_dir))
#     dir.create(hg19_ref_dir, recursive = TRUE)
#   
#   # Save combined CpG Data in BED format
#   bed_file_combined <- file.path(hg19_ref_dir, "hg19_CpG_sites.bed")
#   
#   message("Generating CpG sites data...")
#   
#   chromosomes <- paste0("chr", 1:22)
#   
#   # Initialize an empty list for per-chromosome storage
#   cpg_list <- list()
#   
#   for (chr in chromosomes) {
#     message("Processing ", chr, "...")
#     chr_seq <- BSgenome.Hsapiens.UCSC.hg19[[chr]]
#     
#     # Convert sequence to uppercase to handle masked sequences
#     chr_seq <- toupper(chr_seq)
#     
#     # Find CpG sites
#     cpg_sites <- matchPattern("CG", chr_seq)
#     
#     # Create data frame for this chromosome
#     chr_cpg_df <- data.table(
#       chrom = chr,
#       start = start(cpg_sites),
#       end = end(cpg_sites)
#     )
#     
#     # Store in list
#     cpg_list[[chr]] <- chr_cpg_df
#     
#     # Save individual chromosome CpG sites in BED format
#     bed_file_chr <- file.path(hg19_ref_dir, paste0("hg19_CpG_sites_", chr, ".bed"))
#     header_line <- paste0("#", paste(colnames(chr_cpg_df), collapse = "\t"))
#     writeLines(header_line, bed_file_chr)  # Write header
#     fwrite(chr_cpg_df,
#            file = bed_file_chr,
#            sep = "\t",
#            col.names = FALSE)
#   }
#   
#   # Combine all chromosomes into one data frame
#   cpg_df <- rbindlist(cpg_list)
#   
#   
#   # Dynamically generate header from column names
#   header_line <- paste0("#", paste(colnames(cpg_df), collapse = "\t"))
#   writeLines(header_line, bed_file_combined)  # Write header
#   fwrite(cpg_df,
#          file = bed_file_combined,
#          sep = "\t",
#          col.names = FALSE)
#   
#   
#   message("CpG site extraction complete! Each chromosome saved separately in BED format.")
# }
read_bedmethyl_file <- function(file_path) {
  df <- fread(
    file_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  setnames(
    df,
    c("V1", "V2", "V3", "V4"),
    c(
      "chrom",
      "start",
      "end",
      "Methylation_Percentage"
    )
  )
  return(df)
}

# Process Bismark Methylation Files ----------------------------------------
if (FALSE) {
  file_list <- list.files(path = input_dir,
                          pattern = "_meth.bed",
                          full.names = TRUE)
  
  samples <- sub(".*\\/([A-Z0-9]+)_meth\\.bed$", "\\1", file_list)
  
  missing_samples <- setdiff(meta$sample, samples)
  print(missing_samples)
  
  hg19_CpG_sites <- fread(
    file.path(reference_dir, "hg19_CpG_sites.bed"),
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  n_ref_CpG <- nrow(hg19_CpG_sites)
  
  for (file_path in file_list) {
    # Extract sample ID
    sample_id <- gsub(pattern = "_meth.bed", replacement = "", basename(file_path))
    
    message("Processing new sample: ", sample_id)
    
    # Read methylation data
    sample_data <- read_bedmethyl_file(file_path)
    
    # Filter to canonical chromosomes
    sample_data <- sample_data[chrom %in% paste0("chr", 1:22)]
    
    # Ensure matching number of CpGs
    stopifnot(n_ref_CpG == nrow(sample_data))
    
    # Append methylation percentage as new column
    hg19_CpG_sites[[sample_id]] <- sample_data$Methylation_Percentage
  }
  
  setnames(
    x = hg19_CpG_sites,
    old = c("V1", "V2", "V3"),
    new = c(
      "chrom",
      "start",
      "end")
  )
  
  library(data.table)
  

  
  # List chr1–chr22 files
  all_chr_modern_files <- file.path(modern_dir,"methylation_data",paste0("modern_sample.chr",c(1:22),".bed"))

  # Extract header line from the first file
  header_line <- readLines(all_chr_modern_files[1], n = 1)
  col_names <- strsplit(gsub("^#", "", header_line), "\t")[[1]]
  
  # Read and rbind all
  meth_list <- lapply(all_chr_modern_files, function(f) {
    fread(f, skip = 1)  # skips #chrom header line
  })
  
  meth_all <- rbindlist(meth_list)
  
  # Sanity check
  stopifnot(nrow(meth_all) == nrow(hg19_CpG_sites))
  
  stopifnot(
    all(hg19_CpG_sites$chr == meth_all$chrom),
    all(hg19_CpG_sites$start == meth_all$start),
    all(hg19_CpG_sites$end == meth_all$end)
  )
  
  
  # Assign correct column names
  setnames(meth_all, col_names)
  
  # cbind to hg19_CpG_sites
  final_dt <- cbind(hg19_CpG_sites, meth_all[, -c("chrom", "start", "end")])
  
  saveRDS(
    final_dt,
    file = file.path(dirname(this.dir),"methylation+chip", "all_CpG.45.samples.histogram.hg19.rds")
  )
  
  all_CpG.45.samples.merged.hg19 <- readRDS(file.path(dirname(this.dir),"methylation+chip","all_CpG.45.samples.merged.hg19.rds"))
  
  # Sanity check: matching rows
  stopifnot(
    all(final_dt$chr == all_CpG.45.samples.merged.hg19$chrom)
  )
  
  setDT(all_CpG.45.samples.merged.hg19)
  
  # Combine the columns
  final_dt_extended <- cbind(
    final_dt[, .(chrom, start, end)],
    all_CpG.45.samples.merged.hg19[, .(name, score, strand)],
    final_dt[, -c("chrom", "start", "end")]
  )  

  final_dt_extended$strand <- "."
  
  # Ensure sample columns are in final_dt_extended
  sample_cols <- intersect(meta$sample, names(final_dt_extended))
  
  # Compute row-wise non-NA count across those columns
  final_dt_extended[, score := rowSums(!is.na(.SD)), .SDcols = sample_cols]
  
saveRDS(object = final_dt_extended, 
        file = file.path(dirname(this.dir),"methylation+chip","all_CpG.45.samples.histogram.merged2025.hg19.rds"))


saveRDS(object = final_dt_extended[final_dt_extended$name != "NO_CHIP",],
        file = file.path(dirname(this.dir),"methylation+chip","H3K27ac.only.45.samples.histogram.merged2025.hg19.rds"))
}


# soothin#######################################################

# load ata if not loade

# load the full object back into final_dt_extended
rds_path_full <- file.path(dirname(this.dir), "methylation+chip", "all_CpG.45.samples.histogram.merged2025.hg19.rds")
if (!file.exists(rds_path_full)) stop("File not found: ", rds_path_full)
final_dt_extended <- readRDS(rds_path_full)

# load the filtered object back (into a new variable name)
rds_path_filtered <- file.path(dirname(this.dir), "methylation+chip", "H3K27ac.only.45.samples.histogram.merged2025.hg19.rds")
if (!file.exists(rds_path_filtered)) stop("File not found: ", rds_path_filtered)
H3K27ac.only.45.samples.histogram.merged2025.hg19 <- readRDS(rds_path_filtered)






library(data.table)
library(zoo)  # For moving average

#' Smooth methylation values using a moving window
#'
#' @param df A data.table containing columns: chrom, start, end, Count_Methylated, Count_Unmethylated.
#' @param winsize Size of the smoothing window (default = 5).
#' @param new_col_name Name of the new column storing smoothed values.
#' @return A data.table with an additional column for smoothed methylation values.
smooth_methylation <- function(df, winsize = 5, new_col_name = "smoothed_meth") {
  # Check required columns
  required_cols <- c("chrom", "start", "end", "Count_Methylated", "Count_Unmethylated")
  if (!all(required_cols %in% names(df))) {
    stop("Input data.table must contain columns: chrom, start, end, Count_Methylated, Count_Unmethylated")
  }
  
  # Compute total reads
  df[, total_reads := Count_Methylated + Count_Unmethylated]
  
  # Compute raw methylation level (C / (C + T)), avoiding division by zero
  df[, meth_level := fifelse(total_reads > 0, Count_Methylated / total_reads, NA_real_)]
  
  # Apply a moving average for smoothing
  df[, (new_col_name) := rollapply(meth_level, width = winsize, FUN = mean, na.rm = TRUE, fill = NA, align = "center")]
  
  # Drop intermediate columns
  df[, c("total_reads", "meth_level") := NULL]
  
  # Return modified data.table
  return(df)
}


#estimate original smoothing: #####################################################################

if(FALSE){
  

  
library(data.table)
library(ggplot2)

  # Custom autocorrelation function that ignores NA values
  custom_acf <- function(x, max_lag = 100) {
    n <- length(x)
    acf_vals <- numeric(max_lag + 1)
    
    for (lag in 0:max_lag) {
      valid_indices <- which(!is.na(x[1:(n - lag)]) & !is.na(x[(lag + 1):n]))  # Only valid pairs
      if (length(valid_indices) > 0) {
        acf_vals[lag + 1] <- cor(x[valid_indices], x[valid_indices + lag], use = "pairwise.complete.obs")
      } else {
        acf_vals[lag + 1] <- NA
      }
    }
    
    return(acf_vals)
  }  
  

df_raw <- final_dt_extended[, "S1153"]

# Compute ACF with missing values handled correctly
acf_values <- custom_acf(df_raw, max_lag = 100)

# Create ACF DataTable
acf_df <- data.table(lag = 0:100, acf = acf_values)

# Plot ACF
ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.37, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = paste0("Autocorrelation Function (ACF) of raw data ", biosample),
       x = "Lag (Genomic Positions)",
       y = "Autocorrelation")

# Estimate smoothing window where ACF drops below 0.37
estimated_window <- min(acf_df[acf < 0.37]$lag, na.rm = TRUE)



# # smoothing
# # Load example data
# df <- fread("C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation_data/GSE276666_RAW/methylation_data/S1148/chr1_mapped.bed",
#             header = FALSE, col.names = c("chrom", "start", "end", "Count_Methylated", "Count_Unmethylated"))
# 
# # Apply smoothing function with a window size of 5
# smoothed_df <- smooth_methylation(df, winsize = 21, new_col_name = "S1148")
# 
# # View results
# head(smoothed_df)
# 
# # Compute ACF with missing values handled correctly
# acf_values <- custom_acf(smoothed_df$S1148, max_lag = 100)
# 
# # Create ACF DataTable
# acf_df <- data.table(lag = 0:100, acf = acf_values)
# 
# # Plot ACF
# ggplot(acf_df, aes(x = lag, y = acf)) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = 0.37, linetype = "dashed", color = "red") +
#   theme_minimal() +
#   labs(title = paste0("Autocorrelation Function (ACF) of window size 11 smoothed data ","S1148"),
#        x = "Lag (Genomic Positions)",
#        y = "Autocorrelation")


smoothing_estimates <- data.frame(biosamples = meta$sample,
                                  estimation = NA)


for(biosample in smoothing_estimates$biosamples) {
  print(biosample)
  df_meth <- final_dt_extended[final_dt_extended$chrom == "chr1", ..biosample]
  
  
  # Compute ACF with missing values handled correctly
  acf_values <- custom_acf(df_meth[[1]], max_lag = 100)
  
  # Create ACF DataTable
  acf_df <- data.table(lag = 0:100, acf = acf_values)
  
  # Plot ACF
  ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.37,
               linetype = "dashed",
               color = "red") +
    theme_minimal() +
    labs(
      title = paste0("Autocorrelation Function (ACF) of ", biosample),
      x = "Lag (Genomic Positions)",
      y = "Autocorrelation"
    )
  
  # Estimate smoothing window where ACF drops below 0.37
  estimated_window <- min(acf_df[acf < 0.37]$lag, na.rm = TRUE)
  
  smoothing_estimates$estimation[which(smoothing_estimates$biosamples == biosample)] = estimated_window
  
  # Print estimated smoothing window
  cat("Estimated Smoothing Window (Custom ACF):",
      estimated_window,
      "\n")
}  


## estimate secodn derivative: #############################

library(data.table)
library(ggplot2)

# Create plots folder if it doesn't exist
plot_dir <- file.path(dirname(this.dir), "methylation+chip", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

smoothing_estimates <- data.frame(biosamples = meta$sample,
                                  estimation = NA)

for(biosample in smoothing_estimates$biosamples) {
  print(biosample)
  df_meth <- final_dt_extended[final_dt_extended$chrom == "chr1", ..biosample]
  df_meth <- df_meth[[1]]
  
  # Compute ACF with missing values handled correctly
  acf_values <- custom_acf(df_meth, max_lag = 100)
  
  
  # Create ACF DataTable
  acf_df <- data.table(lag = 0:100, acf = acf_values)
  
  # Compute discrete second derivative (slope change)
  acf_df[, diff1 := c(NA, diff(acf))]
  acf_df[, diff2 := c(NA, diff(diff1))]
  
  # Find lag with most negative second derivative (most curvature)
  elbow_point <- acf_df[which.max(diff2), lag]
  
  # Plot ACF with elbow point
  p <- ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = elbow_point -1 , linetype = "dotted", color = "blue") +
    geom_hline(yintercept = 0.37, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = paste0("Autocorrelation Function (ACF) of Data ", biosample),
      x = "Lag (Genomic Positions)",
      y = "Autocorrelation"
    )
  
  # Save with informative filename
  ggsave(
    filename = file.path(plot_dir, paste0("acf_", biosample, "_chr1.png")),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )
  

  # Save estimation
  smoothing_estimates$estimation[which(smoothing_estimates$biosamples == biosample)] <- elbow_point
  
  cat("Estimated Smoothing Window (Elbow Point):", elbow_point, "\n")
}
}


library(data.table)
library(ggplot2)
library(gridExtra)

# Create plots folder
plot_dir <- file.path(dirname(this.dir), "methylation+chip", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load old dataset
old_data <- readRDS(file.path(dirname(this.dir), "methylation+chip", "all_CpG.39.samples.merged.hg19.rds"))


smoothing_estimates <- data.frame(biosamples = meta$sample, estimation = NA)

for (biosample in smoothing_estimates$biosamples) {
  print(biosample)
  
  # New data ACF
  df_meth_new <- final_dt_extended[chrom == "chr1", ..biosample][[1]]
  acf_values_new <- custom_acf(df_meth_new, max_lag = 100)
  
  acf_df_new <- data.table(lag = 0:100, acf = acf_values_new)
  acf_df_new[, diff1 := c(NA, diff(acf))]
  acf_df_new[, diff2 := c(NA, diff(diff1))]
  elbow_point <- acf_df_new[which.max(diff2), lag]
  
  p_new <- ggplot(acf_df_new, aes(x = lag, y = acf)) +
    geom_line(color = "darkgreen") +
    geom_point() +
    geom_vline(xintercept = elbow_point - 1, linetype = "dotted", color = "blue") +
    geom_hline(yintercept = 0.37, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = paste0("New Histogram-Matched ACF — ", biosample), x = "Lag", y = "ACF")
  
  # Save elbow point
  smoothing_estimates$estimation[smoothing_estimates$biosamples == biosample] <- elbow_point
  
  # Old data ACF (if sample exists)
  if (biosample %in% colnames(old_data)) {
    df_meth_old <- old_data[old_data$chrom == "chr1", biosample]
    acf_values_old <- custom_acf(df_meth_old, max_lag = 100)
    acf_df_old <- data.table(lag = 0:100, acf = acf_values_old)
    
    p_old <- ggplot(acf_df_old, aes(x = lag, y = acf)) +
      geom_line(color = "darkred") +
      geom_point() +
      theme_minimal() +
      labs(title = paste0("Old ACF — ", biosample), x = "Lag", y = "ACF")
    
    # Combine and save both plots
    ggsave(
      filename = file.path(plot_dir, paste0("acf_", biosample, "_chr1_combined.png")),
      plot = grid.arrange(p_new, p_old, ncol = 2),
      width = 10,
      height = 4,
      dpi = 300
    )
  } else {
    # Save only new if old not available
    ggsave(
      filename = file.path(plot_dir, paste0("acf_", biosample, "_chr1_newonly.png")),
      plot = p_new,
      width = 6,
      height = 4,
      dpi = 300
    )
  }
  
  cat("Estimated Smoothing Window (Elbow Point):", elbow_point, "\n")
}


########### same code as grid plto output 

# Minimal edits to collect all per-sample plots into one big grid and save once at the end.
library(data.table)
library(ggplot2)
library(gridExtra)

# Create plots folder
plot_dir <- file.path(dirname(this.dir), "methylation+chip", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load old dataset
old_data <- readRDS(file.path(dirname(this.dir), "methylation+chip", "all_CpG.39.samples.merged.hg19.rds"))

smoothing_estimates <- data.frame(biosamples = meta$sample, estimation = NA)

# NEW: collect combined grobs here
plot_list <- list()

for (biosample in smoothing_estimates$biosamples) {
  
  #biosample <- smoothing_estimates$biosamples[1]
  print(biosample)
  # New data ACF
  df_meth_new <- final_dt_extended[chrom == "chr1", ..biosample][[1]]
  acf_values_new <- custom_acf(df_meth_new, max_lag = 100)
  
  acf_df_new <- data.table(lag = 0:100, acf = acf_values_new)
  acf_df_new[, diff1 := c(NA, diff(acf))]
  acf_df_new[, diff2 := c(NA, diff(diff1))]
  elbow_point <- acf_df_new[which.max(diff2), lag]
  
  p_new <- ggplot(acf_df_new, aes(x = lag, y = acf)) +
    geom_line(color = "darkgreen") +
    geom_point() +
    geom_vline(xintercept = elbow_point - 1, linetype = "dotted", color = "blue") +
    geom_hline(yintercept = acf_df_new$acf[elbow_point], linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = paste0(biosample," Hist. Matched ACF elbow_point: ",elbow_point), x = "Lag", y = "ACF")
  
  # Save elbow point
  smoothing_estimates$estimation[smoothing_estimates$biosamples == biosample] <- elbow_point
  
  # Old data ACF (if sample exists)
  if (biosample %in% colnames(old_data)) {
    df_meth_old <- old_data[old_data$chrom == "chr1", biosample]
    acf_values_old <- custom_acf(df_meth_old, max_lag = 100)
    acf_df_old <- data.table(lag = 0:100, acf = acf_values_old)
    
    p_old <- ggplot(acf_df_old, aes(x = lag, y = acf)) +
      geom_line(color = "darkred") +
      geom_point() +
      theme_minimal() +
      labs(title = paste0(biosample, " Old ACF"), x = "Lag", y = "ACF")
    
    # MINIMAL CHANGE: create a grob (do not save per-sample), append to list
    plot_list[[length(plot_list) + 1]] <- arrangeGrob(p_new, p_old, ncol = 2)
    }else {
      
      # ----------------------------
      # Smoothed data (31-CpG ±15)
      # ----------------------------
      df_meth_new_raw <- final_dt_extended[chrom == "chr1", ..biosample][[1]]
      df_meth_new_ma31 <- zoo::rollmean(df_meth_new_raw, k = 31, align = "center", fill = NA)
      
      acf_values_ma31 <- custom_acf(df_meth_new_ma31[!is.na(df_meth_new_ma31)], max_lag = 100)
      acf_df_ma31 <- data.table(lag = 0:100, acf = acf_values_ma31)
      
      # elbow detection (IDENTICAL logic to p_new)
      acf_df_ma31[, diff1 := c(NA, diff(acf))]
      acf_df_ma31[, diff2 := c(NA, diff(diff1))]
      elbow_point_ma31 <- acf_df_ma31[which.max(diff2), lag]
      
      p_ma31 <- ggplot(acf_df_ma31, aes(x = lag, y = acf)) +
        geom_line(color = "darkgreen") +
        geom_point() +
        geom_vline(
          xintercept = elbow_point_ma31 - 1,
          linetype = "dotted",
          color = "blue"
        ) +
        geom_hline(
          yintercept = acf_df_ma31$acf[elbow_point_ma31],
          linetype = "dashed",
          color = "red"
        ) +
        theme_minimal() +
        labs(
          title = paste0(
            biosample,
            " NEW ACF (31-CpG MA) elbow_point: ",
            elbow_point_ma31
          ),
          x = "Lag",
          y = "ACF"
        )
      
      # ----------------------------
      # Combine: raw vs MA31
      # ----------------------------
      plot_list[[length(plot_list) + 1]] <- arrangeGrob(
        p_ma31,
        p_new,
        ncol = 2,
        widths = c(1, 1)
      )
  }
    
  
  cat("Estimated Smoothing Window (Elbow Point):", elbow_point, "\n")
}

# AFTER loop: arrange all collected grobs into one big grid and save
# Arrange all sample grobs in a 5 x 9 grid (5 columns, 9 rows)
# If fewer/more samples, grid will adapt but preferred layout is 5x9 for 45 samples.
final_grob <- do.call(arrangeGrob, c(plot_list, ncol = 5))

# Save big PNG and PDF. Adjust height per sample if needed.
n_samples <- length(plot_list)

out_png <- file.path(plot_dir, "acf_chr1_all_samples_grid.png")
# width ≈ 5 * 2.4, height ≈ 9 * 2.4 (tweak 2.4 if you want larger/smaller panels)
ggsave(out_png, final_grob, width = 5 * 2.4, height = 9 * 2.4, dpi = 300, limitsize = FALSE)

out_pdf <- file.path(plot_dir, "acf_chr1_all_samples_grid.pdf")
# width ≈ 5 * 2.4, height ≈ 9 * 2.4 (tweak 2.4 if you want larger/smaller panels)
ggsave(out_pdf, final_grob, width = 5 * 7, height = 9 * 4, dpi = 300, limitsize = FALSE)

message("Saved big grid plots:\n - ", out_png, "\n - ", out_pdf)
