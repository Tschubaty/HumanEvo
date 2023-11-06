#################################################################
##
##
##  input:
##
##   ChipSegmentCut
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
OUTPUT_FOLDER <- "ChipSegmentCut"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
ANNOTATION <- "hg19"
N_SAMPLES <- 37


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

CUTVALUE <- 300
meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
real_sample_typisation <-
  readRDS(file.path("04.methylation_vs_Type", "real_sample_typisation.rds"))
meta37 <- meta[meta$sample %in% real_sample_typisation$sample,]

#load_datapoints
load_datapoints <- FALSE
if (load_datapoints) {
  # load data
  df_all <-
    readRDS(
      "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/all_CpG.39.samples.merged.hg19.rds"
    )
  data_meth <- df_all[, real_sample_typisation$sample]
  
  # make smaller data frame
  df_CpG_coordinates <-
    df_all[, c("chrom", "start",   "end",    "name", "score")]
  
  # count number of datapoints
  df_CpG_coordinates$score <-
    apply(
      X = data_meth,
      MARGIN = 1,
      FUN = function(r) {
        return(sum(!is.na(r)))
      }
    )

  saveRDS(object = df_CpG_coordinates,
          file = file.path(OUTPUT_FOLDER, "all_CpG_coordinate.rds"))
} else{
  df_CpG_coordinates <-
    readRDS(file = file.path(OUTPUT_FOLDER, "all_CpG_coordinate.rds"))
}

# paralle processing reigstration
# numCores <- detectCores() - 1
# registerDoParallel(numCores)
df_all_peaks <- data.frame()

for (chr in CHR_NAMES) {
  # start timing chr <-  CHR_NAMES[1]
  print(chr)
  start_time <- Sys.time()
  df_chr <- df_CpG_coordinates[df_CpG_coordinates$chrom == chr,]
  
  df_peaks <- data.frame()
  #foreach(peak = unique(df_chr$name)[-1],
  #       .combine = rbind) %dopar% {
  for (bigpeak in unique(df_chr$name)[-1]) {
    # peak <- "chr1.1079526.1080063.H3K27ac"
    #peak <- "chr1.713321.714754.H3K27ac"
    
    split_name <- unlist(strsplit(x = bigpeak, "[.]"))
    
    peak_len <- as.numeric(split_name[3]) - as.numeric(split_name[2])
    n_parts <- peak_len %/% CUTVALUE
    
    peak_array <- list()
    if (n_parts > 1) {
      for (i in 1:n_parts) {
        peak_array[[i]] <-
          paste(
            split_name[1],
            as.numeric(split_name[2]) + (i - 1) * floor(peak_len / n_parts) + (i - 1),
            as.numeric(split_name[2]) + i * floor(peak_len / n_parts) + (i - 1),
            paste("part",i,"H3K27ac",sep = "_"),
            sep = "."
          )
      }
      peak_array <- unlist(peak_array)
    }else{
      peak_array <- c(bigpeak)
    }
    
    for(peak in peak_array){
      
      split_name <- unlist(strsplit(x = peak, "[.]"))
      
      start = as.numeric(split_name[2])
      end = as.numeric(split_name[3])
      peak_len <- end - start
    indices <- start <=  df_chr$start & df_chr$end <= end
    n_CpG <- sum(indices)
    n_complete_data <-
      sum(N_SAMPLES == df_chr$score[indices])
    before_CpG_coordinate <-
      df_chr$end[min(which(indices)) - 1]
    next_CpG_coordinate <-
      df_chr$end[max(which(indices)) + 1]
    
    max_span <- next_CpG_coordinate - before_CpG_coordinate

    df_chr$name[indices] <- peak
    #return(
    df_peaks <- rbind(
      df_peaks,
      
      data.frame(
        chr = split_name[1],
        start = start,
        end = end,
        name = peak,
        score = n_CpG,
        strand = ".",
        n_complete = n_complete_data,
        before_CpG = before_CpG_coordinate,
        next_CpG_ = next_CpG_coordinate,
        length = peak_len
      )
    )
    
    }
  }
  saveRDS(object = df_peaks,
          file = file.path(
            OUTPUT_FOLDER,
            paste(
              chr,
              "cut",
              "peak_coordinates",
              ANNOTATION,
              N_SAMPLES,
              "n_samples",
              "rds",
              sep = "."
            )
          ))
  
  # stopImplicitCluster()
  gc()
  ######################################################################
  df_all_peaks <- rbind(df_all_peaks, df_peaks)
  df_CpG_coordinates$name[df_CpG_coordinates$chrom == chr] <- df_chr$name
  end_time <- Sys.time()
  print(end_time - start_time)
}
saveRDS(object = df_all_peaks,
        file = file.path(
          OUTPUT_FOLDER,
          paste(
            "all",
            "cut",
            "peak_coordinates",
            ANNOTATION,
            N_SAMPLES,
            "n_samples",
            "rds",
            sep = "."
          )
        ))

saveRDS(object = df_CpG_coordinates,file = file.path(OUTPUT_FOLDER,"all_CpG_coordinates_cut_peak_annotation.rds"))

p <- ggplot(data = df_all_peaks[df_all_peaks$length >= 300,],mapping = aes(x = length,y = score))+
  geom_point()+                            
  stat_smooth(method = "loess", 
              formula = y ~ x)+
  ylab("# CpG in segment")+
  xlab("lenght of segemnt in bp")+
  ggtitle("#CpG vs length")+
  theme_minimal()

ggsave(plot = p,filename = file.path(OUTPUT_FOLDER,"plot", "#CpG vs length.png"))

p <- ggplot(data = df_all_peaks[df_all_peaks$length >= 300,],mapping = aes(x = length))+
  geom_histogram(bins = 100)+                            
  xlab("lenght of segemnt in bp")+
  ggtitle("length hist")+
  theme_minimal()


ggsave(plot = p,filename = file.path(OUTPUT_FOLDER, "plot","cut peak length hist.png"))


if (load_datapoints) {
  
  df_all$name <- df_CpG_coordinates$name
  
  saveRDS(object = df_all,
          file = file.path(OUTPUT_FOLDER, "all_CpG_hg19_cut_peak_annotation.rds"))
  
  saveRDS(object = df_all[df_all$name != "NO_CHIP",],
          file = file.path(OUTPUT_FOLDER, "all_unfiltered_peak_CpG_hg19_cut_peak_annotation.rds"))
  
}

# df_CpG <-
#   readRDS(
#     "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/H3K27ac.only.39.samples.merged.hg19.rds"
#   )
# length(unique(df_CpG$name))
