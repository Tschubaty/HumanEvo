#################################################################
##
##
##  input:
##
##  output:
##
##  v_02 type
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
##################################### CONSTANTS ########################################
INPUT_FOLDER <- "08.genome_segemntation"
OUTPUT_FOLDER <- ""
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

##################################################

#window_width <- 300
min_delta <- 0.10
#offset = 0
INPUT_FOLDER <- "08.genome_segemntation"
OUTPUT_FOLDER <- "09.statistic"

df_peak_CpG <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/04.methylation_vs_Type/peak_CpG_with_kruskal_valis_p_val.rds")

# paralle processing reigstration
numCores <- detectCores() - 1
registerDoParallel(numCores)

# tsart timing
start_time <- Sys.time()

df_peaks <-
  foreach(peak = unique(df_peak_CpG$name), .combine = rbind) %dopar% {
    # peak <- "chr1.1079526.1080063.H3K27ac"
    indices <- df_peak_CpG$name == peak
    split_name <- unlist(strsplit(x = peak, "[.]"))
    n_CpG <- sum(indices)
    type_delta <- df_peak_CpG$type_delta[indices]
    p_val <- df_peak_CpG$KW_p_val[indices]
    # relevant CpGs
    index_CpG_relevant <-
      type_delta >= min_delta & !is.na(p_val)
    
    if (sum(index_CpG_relevant) > 0) {
      start <- df_peak_CpG$start[indices]
      end <- df_peak_CpG$end[indices]
      name <- df_peak_CpG$name[indices]
      score <- df_peak_CpG$score[indices]
      stat <- df_peak_CpG$KW_statistic[indices]
      
      min_p <- min(p_val[index_CpG_relevant])
      min_index <- p_val == min_p

      return(
        data.frame(
          chr = split_name[1],
          start = as.numeric(split_name[2]),
          end = as.numeric(split_name[3]),
          name = peak,
          score = n_CpG,
          strand = ".",
          best_start = start[min_index],
          best_end = end[min_index],
          best_name = name[min_index],
          best_n_sample = score[min_index],
          best_stat = stat[min_index],
          best_type_delta = type_delta[min_index],
          p_val = min_p
        )
      )
    } else{
      return(
        data.frame(
          chr = split_name[1],
          start = as.numeric(split_name[2]),
          end = as.numeric(split_name[3]),
          name = peak,
          score = n_CpG,
          strand = ".",
          best_start = NA,
          best_end = NA,
          best_name = NA,
          best_n_sample = NA,
          best_stat = NA,
          best_type_delta = NA,
          p_val = NA
        )
      )
    }
  }
saveRDS(object = df_peaks,
        file = file.path(
          OUTPUT_FOLDER,
          paste(
            "peak_best_p_val",
            ANNOTATION,
            N_SAMPLES,
            "n_samples",
            min_delta,
            "min_delta",
            "KW",
            "rds",
            sep = "."
          )
        ))


stopImplicitCluster()
gc()
######################################################################

#saveRDS(object = df_peaks,file = file.path(OUTPUT_FOLDER,""))
# 
# Bp_tolerance <- 25 
# 
# df_peaks_for_analysis <- df_peaks[!is.na(df_peaks$best_name),]
# df_peaks_for_analysis$length <- df_peaks_for_analysis$end - df_peaks_for_analysis$start
# 
# ggplot(data = df_peaks_for_analysis, mapping = aes(x = end - start)) +
#   geom_histogram(breaks = seq(form = 0, to = 3000, by = 10)) + 
#   geom_vline(xintercept = 300, color = "red") + 
#   theme_minimal()
# 
# df_peaks_for_analysis_with_window_width <- df_peaks_for_analysis[ abs(df_peaks_for_analysis$length - window_width) <= Bp_tolerance,]
# 
# ggplot(data = df_peaks_for_analysis_with_window_width, mapping = aes(x = score)) +
#   geom_histogram(breaks = seq(form = 0, to = 100, by = 1)) + 
#   theme_minimal()+
#   xlab("# CpG")+
#   ggtitle(paste("peaks with window width" , "=" ,window_width))
# 
# 
# # load data
# df_segmentation <- data.frame()
# for (chr in CHR_NAMES) {
#   df_temp <-readRDS(file.path(
#     INPUT_FOLDER,
#     paste(
#       chr,
#       window_width,
#       "bp",
#       offset,
#       "offset",
#       min_delta,
#       "min_delta",
#       "pearson",
#       "rds",
#       sep = "."
#     )
#   ))
#   df_segmentation <-  rbind(df_segmentation, df_temp[!is.na(df_temp$best_n_sample),])
# }
# 
# # ############################ num CpG to in segemnt 
# 
# out_dir <- paste(window_width,
#                  "bp",
#                  offset,
#                  "offset",
#                  min_delta,
#                  "min_delta",
#                  "pearson",sep = ".")
# 
# dir.create(file.path(OUTPUT_FOLDER, out_dir), showWarnings = FALSE)
# 
# 
# meta_stat <- data.frame(matrix(nrow = length(unique(df_peaks_for_analysis_with_window_width$score)),ncol = 7))
# 
# i <- 0
# for(n_CpG in unique(df_peaks_for_analysis_with_window_width$score)){
# #n_CpG <- 20
# i <- i+1
# 
# df_segmentation_n_CpG <-
#   df_segmentation[df_segmentation$score == n_CpG, ]
# 
# p <- ggplot(data = df_segmentation_n_CpG, mapping = aes(x = -log(p_val))) + geom_histogram() +
#   ggtitle(paste("genome_segemntaion",window_width, "=window_width","\n", n_CpG, "=n_CpG", sep = " "))+theme_minimal()
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   out_dir,
#   paste(
#     "genome_segemntaion",
#     window_width,
#     "=window_width",
#     n_CpG,
#     "=n_CpG",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# p <- ggplot(data = df_peaks_for_analysis_with_window_width[df_peaks_for_analysis_with_window_width$score == n_CpG,], mapping = aes(x = -log(p_val))) + geom_histogram() +
#   ggtitle(paste("H3K27ac",window_width, "=window_width","\n", n_CpG, "=n_CpG", sep = " "))+theme_minimal()
# 
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   out_dir,
#   paste(
#     "H3K27ac_peaks",
#     window_width,
#     "=window_width",
#     n_CpG,
#     "=n_CpG",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# meta_stat$n_CpG[i] <- n_CpG 
# meta_stat$total_fragments[i] <- nrow(df_segmentation_n_CpG)
# meta_stat$NO_CHIP[i] <- sum(df_segmentation_n_CpG$best_name == "NO_CHIP")
# meta_stat$H3K27ac[i] <- sum(df_segmentation_n_CpG$best_name != "NO_CHIP")
# sorted_pval <- sort(df_segmentation_n_CpG$p_val) 
# meta_stat$percentil95[i] <- sorted_pval[length(sorted_pval)*0.05]
# meta_stat$H3K27ac_fraction_in_data[i] <- meta_stat$H3K27ac[i]/meta_stat$total_fragments[i]
# meta_stat$golbal_min_pval[i] <- min(df_segmentation_n_CpG$p_val)
# 
# }
# 
# p <- ggplot(data = meta_stat,mapping = aes(x = n_CpG,y = -log(percentil95)))+geom_point()+theme_minimal()
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#    paste(
#      "95_intervall_",
#     window_width,
#     "=window_width",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# p <- ggplot(data = meta_stat,mapping = aes(x = n_CpG,y = H3K27ac_fraction_in_data))+geom_point() +theme_minimal()
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "H3K27ac_fraction_in_data",
#     window_width,
#     "=window_width",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# 
# p <- ggplot(data = meta_stat,mapping = aes(x = n_CpG,y = -log(golbal_min_pval)))+geom_point() +theme_minimal()
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "golbal_min_pval_nCpg",
#     window_width,
#     "=window_width",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# p <- ggplot(data = meta_stat,mapping = aes(x = log(total_fragments),y = -log(golbal_min_pval)))+geom_point() +theme_minimal()
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "golbal_min_pval_ndata",
#     window_width,
#     "=window_width",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# p <- ggplot(data = meta_stat,mapping = aes(color = n_CpG,x = total_fragments*n_CpG,y = -log(meta_stat$percentil95)))+geom_point() +theme_minimal()
# # ggsave(filename = file.path(
# #   OUTPUT_FOLDER,
# #   paste(
# #     "golbal_min_pval_vs_ndata*nCpG",
# #     window_width,
# #     "=window_width",
# #     ".",
# #     picuture_file_extension,
# #     sep = ""
# #   )
# # ),
# # plot = p)
# 
# # total peaks
# print(nrow(df_peaks))
# # peaks with best CoG
# print(nrow(df_peaks_for_analysis))
# # with min delta 
# print(sum(df_peaks_for_analysis$best_delta_meth >= min_delta ))
# # woth low p-val 
# print(sum(df_peaks_for_analysis$p_val <= 0.01 ))
# # 
# df_peaks_best <-
#   df_peaks_for_analysis[df_peaks_for_analysis$best_delta_meth >= min_delta &
#                           df_peaks_for_analysis$p_val <= 0.01 &
#                           df_peaks_for_analysis$best_n_sample == N_SAMPLES, ]
# print(nrow(df_peaks_best))
# #hist(df_peaks_best$best_delta_meth)
# 
# p <-
#   ggplot(data = df_peaks_best , mapping = aes(x = length, y = -log(p_val))) +
#   geom_point() + ggtitle("min delata", min_delta)
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "H3K27AC",
#     "min delata",
#     min_delta,
#     "lenght_vs_-logP",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# p <- ggplot(data = df_peaks_for_analysis ,mapping = aes(x = length,y = -log(p_val)))+geom_point()
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "H3K27AC",
#     "all_peaks",
#     "lenght_vs_-logP",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# p <- ggplot(data = df_peaks_best ,mapping = aes(x = score,y = -log(p_val)))+geom_point()+xlab("# CpG")+ggtitle("min delata", min_delta)
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "H3K27AC",
#     "min delata",
#     min_delta,
#     "n_CpG_vs_-logP",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# p <- ggplot(data = df_peaks_for_analysis ,mapping = aes(x = score,y = -log(p_val)))+geom_point()+xlab("# CpG")
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "H3K27AC",
#     "all_peak",
#     "n_CpG_vs_-logP",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# p <- ggplot(data = df_peaks_best ,mapping = aes(x = score/length,y = -log(p_val)))+geom_point()+xlab("CpG density") +ggtitle("min delata", min_delta)
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "H3K27AC",
#     "min delata",
#     min_delta,
#     "CpG density_vs_-logP",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# p <- ggplot(data = df_peaks_for_analysis ,mapping = aes(x = score/length,y = -log(p_val)))+geom_point()+xlab("CpG density")
# ggsave(filename = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "H3K27AC",
#     "all_peaks",
#     "CpG density_vs_-logP",
#     ".",
#     picuture_file_extension,
#     sep = ""
#   )
# ),
# plot = p)
# 
# saveRDS(object = df_peaks_best,file = file.path(OUTPUT_FOLDER,
#         paste(
#           "H3K27AC",
#           "best_peaks",
#           "alpha_0.01",
#           "rds",
#           sep = "."
#         )))
# #ggstatsplot::ggscatterstats(data = df_peaks_best, x = length, y = -log(p_val))
# 
# 
# # library(ggpubr)
# # ggplot(data = df_peaks_for_analysis ,mapping = aes(x = score/length,y = -log(p_val)))+geom_point()+xlab("CpG density")+stat_cor(method="pearson")