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

# Set working directory
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
      
      sample_data_chr[, `:=`(
        name = ".",     # Placeholder for the 'name' field
        score = 0,      # Placeholder for the 'score' field
        strand = "."    # Placeholder for the 'strand' field
      )]
      
      # Reorder columns to match BED format
      sample_data_chr <- sample_data_chr[, .(chrom, start, end, name, score, strand, Methylation_Percentage, Count_Methylated, Count_Unmethylated)]
      
      
    } else {
      print(paste0(chr,"not canonical chromosem"))
      sample_data_chr <- NA

    }
    
    # Save processed chromosome data in bed format
    header_line <- paste0("#", paste(colnames(sample_data_chr), collapse = "\t"))
    writeLines(header_line, chr_output_file)  # Write header
    fwrite(sample_data_chr,
           file = chr_output_file,
           sep = "\t",
           col.names = FALSE)
  }
  
  # Clean up memory
  # rm(sample_data, merged_chr)
  gc()
}

# Final Message
message("Processing complete! Each chromosome saved separately in sample-specific folders.")
#####################################################################################################################################

# # **Optimized Function: Merge CpG Data in Linear Time with Debugging**
# sum_methylation_fast <- function(cpg_chr, sample_chr) {
#   # **Initialize Progress Bar**
#   pb <- txtProgressBar(min = 0,
#                        max = nrow(cpg_chr),
#                        style = 3)
#   
#   # debug:
#   # cpg_chr <- cpg_list[[chr]]
#   # sample_chr <- sample_list[[chr]]
#   #
#   # sample_chr[(sample_index - 2):(sample_index + 2),]
#   # test_indx <- which(sample_chr$start[sample_index-5] <= cpg_chr$start   & cpg_chr$end  <= sample_chr$start[sample_index+5])
#   #
#   # cpg_chr[test_indx,]
#   
#   # **Ensure both inputs contain only one chromosome**
#   unique_cpg_chrom <- unique(cpg_chr$chrom)
#   unique_sample_chrom <- unique(sample_chr$chrom)
#   
#   if (length(unique_cpg_chrom) > 1 ||
#       length(unique_sample_chrom) > 1 ||
#       unique_cpg_chrom != unique_sample_chrom) {
#     stop(
#       "Error: cpg_chr and sample_chr must contain data for exactly one and the same chromosome."
#     )
#   }
#   
#   #message("Processing Chromosome: ", unique_cpg_chrom)
#   
#   # **Ensure both data.tables are sorted**
#   setkey(cpg_chr, chrom, start)
#   setkey(sample_chr, chrom, start)
#   
#   # **Initialize pointers**
#   sample_index <- 1
#   n_samples <- nrow(sample_chr)
#   
#   # **Initialize results table**
#   cpg_chr[, `:=`(Count_Methylated = 0,
#                  Count_Unmethylated = 0)]
#   
#   # **Loop through CpG sites (linear time merge)**
#   for (i in seq_len(nrow(cpg_chr))) {
#     # **Update Progress Bar**
#     setTxtProgressBar(pb, i)
#     
#     reference_indx_start <- cpg_chr$start[i]
#     reference_indx_end <- cpg_chr$end[i]
#     
#     # Debug message for current CpG site
#     #message("Processing CpG site: ", reference_indx_start, "-", reference_indx_end)
#     
#     # Initialize an empty data.table to store unmatched rows
#     #unmatched_samples <- data.table()
#     
#     # **Collect unmatched rows instead of throwing an error**
#     while (sample_index <= n_samples &&
#            sample_chr$end[sample_index] < reference_indx_start) {
#       # message(paste("Error: Sample CpG missed:", sample_chr$start[sample_index],
#       #            "exceeds CpG region", reference_indx_start, "-", reference_indx_end,
#       #            "at sample index", sample_index," and cpg_chr index ",i))
#       
#       # Append the unmatched row to `unmatched_samples`
#       #unmatched_samples <- rbind(unmatched_samples, sample_chr[sample_index, ], fill = TRUE)
#       
#       # Move to the next sample row
#       sample_index <- min(sample_index + 1, n_samples)
#     }
#     
#     
#     while (sample_chr$start[sample_index] <= reference_indx_end) {
#       # **Accumulate counts**
#       cpg_chr$Count_Methylated[i] <- cpg_chr$Count_Methylated[i] + sample_chr$Count_Methylated[sample_index]
#       cpg_chr$Count_Unmethylated[i] <- cpg_chr$Count_Unmethylated[i] + sample_chr$Count_Unmethylated[sample_index]
#       
#       # **Move to next sample**
#       sample_index <- min(sample_index + 1, n_samples)  # Prevent out-of-bounds access
#     }
#   }
#   
#   # **Handle sites with no reads (set to NA instead of 0)**
#   cpg_chr[Count_Methylated == 0 & Count_Unmethylated == 0, `:=`(
#     Count_Methylated = NA,
#     Count_Unmethylated = NA,
#     Methylation_Percentage = NA
#   )]
#   
#   # **Compute Methylation Percentage**
#   cpg_chr[, Methylation_Percentage := ifelse(Count_Methylated + Count_Unmethylated > 0,
#                                              (Count_Methylated / (Count_Methylated + Count_Unmethylated)) * 100,
#                                              NA)]
#   
#   # **Close Progress Bar**
#   close(pb)
#   
#   return(cpg_chr)
# }


#' # Process Bismark Methylation Files ----------------------------------------
#' file_list <- list.files(path = this.dir,
#'                         pattern = "deduplicated.bismark.cov$",
#'                         full.names = TRUE)
#' 
#' for (file_path in file_list) {
#'   # Extract sample ID
#'   sample_id <- gsub(".*?(S[0-9]+).*", "\\1", basename(file_path))
#'   sample_dir <- file.path(output_dir, sample_id)  # Create sample-specific folder
#'   
#'   message("Processing new sample: ", sample_id)
#'   
#'   # Create folder for the sample
#'   dir.create(sample_dir, recursive = TRUE)
#'   
#'   # Read methylation data
#'   sample_data <- read_bismark_file(file_path)
#'   
#'   # Convert to `data.table` for speed
#'   setDT(cpg_df)
#'   setDT(sample_data)
#'   
#'   # Split by chromosome
#'   cpg_list <- split(cpg_df, cpg_df$chrom)
#'   sample_list <- split(sample_data, sample_data$chrom)
#'   
#'   for (chr in names(cpg_list)) {
#'     message("Processing chromosome: ", chr)
#'     
#'     chr_output_file <- file.path(sample_dir, paste0(chr, ".rds"))
#'     
#'     # if (file.exists(chr_output_file)) {
#'     #   message("Chromosome ", chr, " already processed: Skipping")
#'     #   next
#'     # }
#'     
#'     if (chr %in% names(sample_list)) {
#'       # Use fast summation method
#'       merged_chr <- sum_methylation_fast(cpg_list[[chr]], sample_list[[chr]])
#'     } else {
#'       # If no sample data for this chromosome, keep NA values
#'       merged_chr <- cpg_list[[chr]]
#'       merged_chr[, `:=`(
#'         Count_Methylated = NA,
#'         Count_Unmethylated = NA,
#'         Methylation_Percentage = NA
#'       )]
#'     }
#'     
#'     # Save processed chromosome data
#'     saveRDS(merged_chr, file = chr_output_file)
#'   }
#'   
#'   # Clean up memory
#'   # rm(sample_data, merged_chr)
#'   gc()
#' }
#' 
#' # Final Message
#' message("Processing complete! Each chromosome saved separately in sample-specific folders.")
#' 
#' 
#' #' Smooth Methylation Data Using Moving Average Filter
#' #'
#' #' This function applies a smoothing window to methylation counts using a moving average filter.
#' #'
#' #' @param chrom_data A data table containing columns `Count_Methylated` and `Count_Methylated + Count_Unmethylated`.
#' #' @param winsize The window size for smoothing (must be an odd integer).
#' #'
#' #' @return A data table with smoothed methylation counts and total reads.
#' #'
#' #' @importFrom stats filter
#' #' @export
#' smooth_methylation <- function(chrom_data, winsize) {
#'   if (winsize %% 2 == 0)
#'     stop("Window size must be an odd integer")
#'   
#'   # Define convolution kernel (equal weights for moving average)
#'   kernel <- rep(1 / winsize, winsize)
#'   
#'   # Apply moving average filter (preserve NA values)
#'   chrom_data[, Smoothed_Count_Methylated := stats::filter(Count_Methylated,
#'                                                           kernel,
#'                                                           sides = 2,
#'                                                           circular = FALSE)]
#'   chrom_data[, Smoothed_Total_Reads := stats::filter(
#'     Count_Methylated + Count_Unmethylated,
#'     kernel,
#'     sides = 2,
#'     circular = FALSE
#'   )]
#'   
#'   return(chrom_data)
#' }
