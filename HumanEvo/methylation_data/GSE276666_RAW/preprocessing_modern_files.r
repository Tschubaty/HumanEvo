# Documentation -----------------------------------------------------------
##
##  input: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276666
##
##  output: Methylation files in hg19 with + and - strands combined
##
##  v_04 16.03.2025
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

# Set output directory
output_dir <- file.path(this.dir, "methylation_data")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Load or Create CpG Sites Data --------------------------------------------
if (FALSE) {
  # Ensure output directory and subfolder exist
  hg19_ref_dir <- file.path(output_dir, "hg19reference")
  if (!dir.exists(hg19_ref_dir))
    dir.create(hg19_ref_dir, recursive = TRUE)
  
  # Save combined CpG Data in BED format
  bed_file_combined <- file.path(hg19_ref_dir, "hg19_CpG_sites.bed")
  
  message("Generating CpG sites data...")
  
  chromosomes <- paste0("chr", 1:22)
  
  # Initialize an empty list for per-chromosome storage
  cpg_list <- list()
  
  for (chr in chromosomes) {
    message("Processing ", chr, "...")
    chr_seq <- BSgenome.Hsapiens.UCSC.hg19[[chr]]
    
    # Convert sequence to uppercase to handle masked sequences
    chr_seq <- toupper(chr_seq)
    
    # Find CpG sites
    cpg_sites <- matchPattern("CG", chr_seq)
    
    # Create data frame for this chromosome
    chr_cpg_df <- data.table(
      chrom = chr,
      start = start(cpg_sites),
      end = end(cpg_sites)
    )
    
    # Store in list
    cpg_list[[chr]] <- chr_cpg_df
    
    # Save individual chromosome CpG sites in BED format
    bed_file_chr <- file.path(hg19_ref_dir, paste0("hg19_CpG_sites_", chr, ".bed"))
    header_line <- paste0("#", paste(colnames(chr_cpg_df), collapse = "\t"))
    writeLines(header_line, bed_file_chr)  # Write header
    fwrite(chr_cpg_df,
           file = bed_file_chr,
           sep = "\t",
           col.names = FALSE)
  }
  
  # Combine all chromosomes into one data frame
  cpg_df <- rbindlist(cpg_list)
  
  
  # Dynamically generate header from column names
  header_line <- paste0("#", paste(colnames(cpg_df), collapse = "\t"))
  writeLines(header_line, bed_file_combined)  # Write header
  fwrite(cpg_df,
         file = bed_file_combined,
         sep = "\t",
         col.names = FALSE)
  
  
  message("CpG site extraction complete! Each chromosome saved separately in BED format.")
}
read_bismark_file <- function(file_path) {
  df <- fread(
    file_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  setnames(
    df,
    c("V1", "V2", "V3", "V4", "V5", "V6"),
    c(
      "chrom",
      "start",
      "end",
      "Methylation_Percentage",
      "Count_Methylated",
      "Count_Unmethylated"
    )
  )
  
  # **Ensure chromosome format matches "chr1", "chr2", etc.**
  df[, chrom := paste0("chr", chrom)]  #
  return(df)
}

# Process Bismark Methylation Files ----------------------------------------
if (FALSE) {
  file_list <- list.files(path = this.dir,
                          pattern = "deduplicated.bismark.cov$",
                          full.names = TRUE)
  
  for (file_path in file_list) {
    # Extract sample ID
    sample_id <- gsub(".*?(S[0-9]+).*", "\\1", basename(file_path))
    sample_dir <- file.path(output_dir, sample_id)  # Create sample-specific folder
    
    message("Processing new sample: ", sample_id)
    
    # Create folder for the sample
    dir.create(sample_dir, recursive = TRUE)
    
    # Read methylation data
    sample_data <- read_bismark_file(file_path)
    
    # Convert to `data.table` for speed
    setDT(sample_data)
    
    # Split by chromosome
    sample_list <- split(sample_data, sample_data$chrom)
    
    for (chr in paste0("chr", 1:22)) {
      message("Processing chromosome: ", chr)
      
      chr_output_file <- file.path(sample_dir, paste0(chr, ".bed"))
      
      # if (file.exists(chr_output_file)) {
      #   message("Chromosome ", chr, " already processed: Skipping")
      #   next
      # }
      
      if (chr %in% names(sample_list)) {
        # Assuming sample_data is a data.table
        # Add placeholder columns
        
        sample_data_chr <- sample_list[[chr]]
        
        setDT(sample_data_chr)
        
        sample_data_chr[, `:=`(name = ".",
                               # Placeholder for the 'name' field
                               score = 0,
                               # Placeholder for the 'score' field
                               strand = "."    # Placeholder for the 'strand' field
                               )
                               ]
                               
                               # Reorder columns to match BED format
                               sample_data_chr <- sample_data_chr[, .(
                                 chrom,
                                 start,
                                 end,
                                 name,
                                 score,
                                 strand,
                                 Methylation_Percentage,
                                 Count_Methylated,
                                 Count_Unmethylated
                               )]
                               
                               
      } else {
        print(paste0(chr, "not canonical chromosem"))
        sample_data_chr <- NA
        
      }
      
      # Save processed chromosome data in bed format
      header_line <- paste0("#", paste(colnames(sample_data_chr), collapse = "\t"))
      writeLines(header_line, chr_output_file)  # Write header
      fwrite(
        sample_data_chr,
        file = chr_output_file,
        sep = "\t",
        col.names = FALSE
      )
    }
    
    # Clean up memory
    # rm(sample_data, merged_chr)
    gc()
  }
  
  # Final Message
  message("Processing complete! Each chromosome saved separately in sample-specific folders.")
}

dirs <- list.dirs(path = file.path(this.dir, "methylation_data"),
                  recursive = FALSE)
sample_dirs <- grep("^S.*", basename(dirs), value = TRUE)

# merge to on edataframe #########################################
if (FALSE) {
  
  library(data.table)
  library(zoo)  # For moving average
  
  #' Smooth methylation values using a moving window
  #'
  #' @param df A data.table containing columns: chrom, start, end, Count_Methylated, Count_Unmethylated.
  #' @param winsize Size of the smoothing window (default = 5).
  #' @param new_col_name Name of the new column storing smoothed values.
  #' @return A data.table with an additional column for smoothed methylation values.
  smooth_methylation <- function(df, winsize = 11, new_col_name = "smoothed_meth") {
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
  
  
  # Set output directory
  output_dir <- "C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation_data_chr_merged"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  

  
  
  for (chr in paste0("chr", 1:22)) {
    # chr <- paste0("chr", 1:22)[22]
    message("Processing chromosome: ", chr)
    
    df_modern <- data.frame()
    
    for (biosample in sample_dirs) {
      # biosample <-  sample_dirs[1]
      
      
      
      chr_file <- file.path(this.dir,
                            "methylation_data",
                            biosample,
                            paste0(chr, "_mapped.bed"))
      chr_output_file <- file.path(this.dir,
                                   "methylation_data",
                                   biosample,
                                   paste0(biosample, ".", chr, ".bed"))
      
      df <- read.delim(
        chr_file,
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE
      )
      print(dim(df))
      colnames(df) <- c("chrom",
                        "start",
                        "end" ,
                        "Count_Methylated",
                        "Count_Unmethylated")
      
      
      df[, biosample] <- df$Count_Methylated / (df$Count_Methylated + df$Count_Unmethylated)
      
      # Save processed chromosome data in bed format
      header_line <- paste0("#", paste(colnames(df[, c("chrom", "start", "end", biosample)]), collapse = "\t"))
      writeLines(header_line, chr_output_file)  # Write header
      fwrite(df[, c("chrom", "start", "end", biosample)],
             file = chr_output_file,
             sep = "\t",
             col.names = FALSE,
             append = TRUE)
      
      if (biosample == sample_dirs[1]) {
        df_modern <- df[, c("chrom", "start", "end", biosample)]
      } else{
        #df_modern <- cbind(df_modern, biosample = df[, biosample])
        df_modern[[biosample]] <- df[[biosample]]
      }
      
    }
    chr_output_file <- file.path(this.dir,
                                 "methylation_data",
                                 paste0("modern_sample", ".", chr, ".bed"))
    header_line <- paste0("#", paste(colnames(df_modern), collapse = "\t"))
    writeLines(header_line, chr_output_file)  # Write header
    fwrite(df_modern,
           file = chr_output_file,
           sep = "\t",
           col.names = FALSE,
           append = TRUE)
    
    
  }
  # Clean up memory
  # rm(sample_data, merged_chr)
  gc()
}

if (FALSE) {
  # Set output directory
  output_dir <- "C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation_data_chr_merged"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  input_dir_modern <- file.path(this.dir,
                                "methylation_data")

  
  for (chr in paste0("chr", 1:22)) {
    # chr <- paste0("chr", 1:22)[22]
    message("Processing chromosome: ", chr)
    
    df_modern <- read.delim(file = file.path(input_dir_modern,paste0("modern_sample.",chr,".bed")))
    df_ancient <- readRDS(file = file.path(output_dir,paste0(chr,".39.samples.rds")))
    df_combined <- cbind(df_ancient,df_modern[, -(1:3)])
    saveRDS(object = df_combined,file = file.path(output_dir,paste0(chr,".45.samples.rds")))  
    
  }
  # Clean up memory
  # rm(sample_data, merged_chr)
  gc()
}

if(FALSE) {
  df_combined_modern <- data.frame()
  for (chr in paste0("chr", 1:22)) {
    # chr <- paste0("chr", 1:22)[22]
    message("Processing chromosome: ", chr)
    
    df_combined <- read.delim(file = file.path(output_dir, paste0("modern_sample.",chr, ".bed")))
    df_combined_modern <- rbind(df_combined_modern, df_combined)
  }
  # Clean up memory
  # rm(sample_data, merged_chr)
  gc()
  
}

all_CpG.39.samples.merged.hg19 <- readRDS("C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/all_CpG.39.samples.merged.hg19.rds")
all_CpG.45.samples.merged.hg19 <- cbind(all_CpG.39.samples.merged.hg19, df_combined_modern[,sample_dirs])
saveRDS(object = all_CpG.45.samples.merged.hg19, 
        file = "C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/all_CpG.45.samples.merged.hg19.rds")


saveRDS(object = all_CpG.45.samples.merged.hg19[all_CpG.45.samples.merged.hg19$name != "NO_CHIP",],
        file = "C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/H3K27ac.only.45.samples.merged.hg19.rds")

# soothin#######################################################
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
  
  
all_CpG.39.samples.merged.hg19 <- readRDS("C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/all_CpG.39.samples.merged.hg19.rds")
raw_data <- read.delim(
  "C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation_data/GSE276666_RAW/methylation_data/S1148/S1148.chr1.bed",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)
df_raw <- raw_data$S1148

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
  labs(title = paste0("Autocorrelation Function (ACF) of raw data ","S1148"),
       x = "Lag (Genomic Positions)",
       y = "Autocorrelation")

# Estimate smoothing window where ACF drops below 0.37
estimated_window <- min(acf_df[acf < 0.37]$lag, na.rm = TRUE)



# smoothing
# Load example data
df <- fread("C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation_data/GSE276666_RAW/methylation_data/S1148/chr1_mapped.bed",
            header = FALSE, col.names = c("chrom", "start", "end", "Count_Methylated", "Count_Unmethylated"))

# Apply smoothing function with a window size of 5
smoothed_df <- smooth_methylation(df, winsize = 21, new_col_name = "S1148")

# View results
head(smoothed_df)

# Compute ACF with missing values handled correctly
acf_values <- custom_acf(smoothed_df$S1148, max_lag = 100)

# Create ACF DataTable
acf_df <- data.table(lag = 0:100, acf = acf_values)

# Plot ACF
ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.37, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = paste0("Autocorrelation Function (ACF) of window size 11 smoothed data ","S1148"),
       x = "Lag (Genomic Positions)",
       y = "Autocorrelation")


smoothing_estimates <- data.frame(biosamples = colnames(all_CpG.39.samples.merged.hg19)[7:ncol(all_CpG.39.samples.merged.hg19)],
                                  estimation = NA)
for(biosample in smoothing_estimates$biosamples) {
  print(biosample)
  df_meth <- all_CpG.39.samples.merged.hg19[all_CpG.39.samples.merged.hg19$chrom == "chr1", biosample]
  
  
  # Compute ACF with missing values handled correctly
  acf_values <- custom_acf(df_meth, max_lag = 100)
  
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
      title = paste0("Autocorrelation Function (ACF) of Smoothed Data ", biosample),
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

smoothing_estimates <- data.frame(
  biosamples = colnames(all_CpG.39.samples.merged.hg19)[7:ncol(all_CpG.39.samples.merged.hg19)],
  estimation = NA
)

for(biosample in smoothing_estimates$biosamples) {
  print(biosample)
  df_meth <- all_CpG.39.samples.merged.hg19[all_CpG.39.samples.merged.hg19$chrom == "chr1", biosample]
  
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
  ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = elbow_point -1 , linetype = "dotted", color = "blue") +
    geom_hline(yintercept = 0.37, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = paste0("Autocorrelation Function (ACF) of Smoothed Data ", biosample),
      x = "Lag (Genomic Positions)",
      y = "Autocorrelation"
    )
  
  # Save estimation
  smoothing_estimates$estimation[which(smoothing_estimates$biosamples == biosample)] <- elbow_point
  
  cat("Estimated Smoothing Window (Elbow Point):", elbow_point, "\n")
}
}
