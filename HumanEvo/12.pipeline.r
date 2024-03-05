# Documentation -----------------------------------------------------------
##
##
##  input: methylation+chip
##
##
##  output: pipeline
##
##
##  v_01 29.11.2023
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------

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
# Load Libraries ----------------------------------------------------------
library(foreach)
# library(doParallel)
library(ggplot2)
#library(Rtsne)
#library(ggrepel)
library("snow")
library("parallel")
# Load the cowplot library
library(cowplot)
library(tidyr)

# set constants  ----------------------------------------------------------

#supplemntary material https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02572-z
df_universal_annotation <-
  readxl::read_xlsx("13059_2021_2572_MOESM3_ESM.xlsx",
                    sheet = "characterize_full_stack_state",
                    n_max = 100)
# Specify the columns you want to import
columns_to_import <-
  c(
    "state_order_by_group",
    "States  presented in paper",
    "Long annotations",
    "Group",
    "color"
  )
# Read the specified columns from the Excel file
df_universal_annotation <-
  df_universal_annotation[, columns_to_import]

df_universal_annotation$color[df_universal_annotation$color == "light puple"] <-
  "#8A91D0"
df_universal_annotation$color[df_universal_annotation$color == "silver"] <-
  "#808080"
df_universal_annotation$color[df_state_summery$state_group == "weak transcription"] <-
  "#009600"
df_universal_annotation$color[df_state_summery$state_group == "transcription"] <-
  "#008000"

# matching_indices <- match(df_state_summery$state,paste(df_universal_annotation$state_order_by_group,df_universal_annotation$`States  presented in paper`,sep = "_"))
# df_state_summery$state_color <- df_universal_annotation$color[matching_indices]

#  "white"
# "light puple"
#
# "lemonchiffon" "#FFFFCC"
# "yellow"       "#FFFF00"
# "orange"
# "#ADFF2F"
# "green"
# "green"
# "#3cb371"
# "#7fffd4"
# "#fff44f"
# "purple"
# "#ff4500"
# "red"

meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
meta <- meta[order(meta$age_mean_BP),]
meta$sample <- factor(x = meta$sample, levels = meta$sample)

real_age_typisation <-
  meta[order(meta$age_mean_BP), c("sample", "age_mean_BP")]
real_food_typisation <-
  meta[meta$Type %in% TYPES, c("sample", "Type")]

## only big sample groups
TYPES <- c("Farmer", "Steppe", "HG")

OUTPUT_FOLDER <- "12.pipeline"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
INPUT_FOLDER <- "methylation+chip"
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
pearson_test_colnames <-
  c("pearson.p_val", "pearson.statistic", "pearson.delta")
KW_test_colnames <- c("kw.p_val", "kw.statistic", "kw.delta")

# set running parameters -------------------------------------------------------
H3K27ac_Analysis <- FALSE
create_true_stat <- FALSE
describe_Data <- FALSE
create_permutations_horizonal_columns <- FALSE
create_CpG_permutations_vertical <- FALSE
summerize_CpG_columns_permutations <- FALSE
summerize_horizonal_columns_permutations <- FALSE
summerize_real_values <- FALSE
histogramm_plots <- FALSE
aggregate_CpG_results <- FALSE
aggregate_horizontal_results <- FALSE
plot_kw_landscape <- FALSE
pearson_landscape <- FALSE
OTHER <- FALSE
gene_annotation <- FALSE
plot_PCA <- FALSE
DO_enrichment_analysis <- FALSE

# load Data ---------------------------------------------------------------
# load(".RData")
df_peak_CpG <-
  readRDS(
    "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/H3K27ac.only.39.samples.merged.hg19.rds"
  )

# define functions ---------------------------------------------------------

# Function to pad numerical prefixes with leading zeros
pad_with_zeros <- function(x) {
  parts <- strsplit(x, "_")
  nums <- as.numeric(sapply(parts, "[[", 1))
  max_digits <- max(nchar(nums))
  sprintf(paste0("%0", max_digits, "d_%s"), nums, sapply(parts, "[[", 2))
}

#' computes maximum group means between dataset typisations
#'
#' @param df_row a row with all methlyation values
#' with the appropriate sample column names.
#' @param typisation dataframe storing sample names in $sample
#' and group types in $type.
#' @returns A numeric vector with NA for missing values .
delta_mean_type <- function(df_row, typisation) {
  typisation$meth <- NA
  for (ty in unique(typisation$Type)) {
    meth <-
      as.numeric(df_row[as.character(typisation$sample[typisation$Type == ty])])
    # test prerequisite if at least 1 sample
    if (sum(!is.na(meth)) < 1) {
      return(NA)
    }
    typisation$meth[typisation$Type == ty] <- meth
  }
  type_means <-
    aggregate(
      typisation$meth,
      list(typisation$Type),
      FUN = function(x) {
        mean(x, na.rm = TRUE)
      }
    )
  return(max(type_means$x) - min(type_means$x))
}

#' computes kruskal.test on data.frame row of methylation
#'
#' @param df_row a row with all methlyation values
#' with the appropriate sample column names.
#' @param typisation dataframe storing sample names in $sample
#' and group types in $type one can make permuations
#' by permuting $type or $sample in the typisation parameter
#' @returns A kruskal.test object
meth_vs_type_test <- function(df_row, typisation) {
  typisation$meth <- NA
  for (s in unique(typisation$Type)) {
    meth <- as.numeric(df_row[as.character(typisation$sample[typisation$Type == s])])
    # test prerequisite if at least 5 samples
    if (sum(!is.na(meth)) < 5) {
      return(NA)
    }
    typisation$meth[typisation$Type == s] <- meth
  }
  
  kru_test <- kruskal.test(meth ~ Type, data = typisation)
  return(kru_test)
}


#' computes kruskal.test on data.frame of methylation
#' function is parallized
#' @param df data frame with all methlyation values
#' with the appropriate sample column names.
#' @param food.typisation dataframe storing sample names in $sample
#' and group types in $type one can make permuations
#' by permuting $type or $sample in the food.typisation parameter
#' @returns A data.frame with the follwing columns: kw.p_val kw.statistic kw.delta
parallel_testing_kruskall_valis <- function(df, food.typisation) {
  clust <- makeCluster(parallel::detectCores() - 1)
  # Export it form base to workspace
  clusterExport(clust, "meth_vs_type_test")
  clusterExport(clust, "delta_mean_type")
  clusterExport(clust, deparse(substitute(food.typisation)))
  # cpmpute kw test
  print("Kw test")
  simulated_test <-
    parApply(
      cl = clust,
      X = df,
      MARGIN = 1,
      FUN = function(r) {
        return(meth_vs_type_test(r, food.typisation))
      }
    )
  # compute maximum group delta
  print("Kw delta")
  delta <-
    parApply(
      cl = clust,
      X = df,
      MARGIN = 1,
      FUN = function(r) {
        return(delta_mean_type(df_row = r, typisation = food.typisation))
      }
    )
  
  # extract pval from test
  print("KW pval")
  p_val <-
    as.numeric(
      parSapply(
        cl = clust,
        X = simulated_test,
        FUN = function(t) {
          ifelse(test = is.na(t),
                 yes = return(NA),
                 no = return(t$p.value))
        },
        USE.NAMES = FALSE,
        simplify = TRUE
      )
    )
  # extract statistic from test
  print("KW stat")
  statistic <-
    as.numeric(
      parSapply(
        cl = clust,
        X = simulated_test,
        FUN = function(t) {
          ifelse(test = is.na(t),
                 yes = return(NA),
                 no = return(t$statistic))
        },
        USE.NAMES = FALSE,
        simplify = TRUE
      )
    )
  stopCluster(clust)
  gc()
  df_x <-
    data.frame(kw.p_val = p_val,
               kw.statistic = statistic,
               kw.delta = delta)
  print("KW complete")
  return(df_x)
}

#' computes pearson cor test on data.frame of methylation
#' function is parallized
#' @param df data frame with all methlyation values
#' with the appropriate sample column names.
#' @param age.typisation dataframe storing sample names in $sample
#' and age in $age_mean_BP one can make permuations
#' by permuting $age_mean_BP or $sample in the age.typisation parameter
#' @returns A data.frame with the follwing columns: pearson.p_val pearson.statistic pearson.delta
parallel_testing_pearson_cor <- function(df, age.typisation) {
  #library("snow")
  #Create cluster
  clus <- parallel::makeCluster(parallel::detectCores() - 1)
  # Export to workspace
  clusterExport(clus, deparse(substitute(age.typisation)))
  # get meth values
  df_meth <-
    df[, as.character(age.typisation$sample[order(age.typisation$age_mean_BP)])]
  # get numeric  age vector
  age <-
    age.typisation$age_mean_BP[order(age.typisation$age_mean_BP)]
  
  # Export  to workspace
  # clusterExport(clus, deparse(substitute(age)))
  
  print("pearson delta")
  # delta = difference between max and min recorded methylation score
  pearson.delta <-
    parRapply(
      cl = clus,
      x = df_meth,
      FUN = function(r) {
        return(max(r, na.rm = TRUE) - min(r, na.rm = TRUE))
      }
    )
  
  print("pearson test")
  # compute peasron correlations
  psr_cor_tests <-
    parRapply(
      cl = clus,
      x = df_meth,
      FUN = function(r) {
        if (sum(!is.na(r)) < 3) {
          return(NA)
        } else{
          return(cor.test(
            x = age,
            y =  as.numeric(r),
            method = "pearson"
          ))
        }
      }
    )
  
  print("pearson pval")
  # extract pval
  pearson.p_val <- sapply(
    psr_cor_tests ,
    FUN = function(r) {
      if (is.na(r)) {
        return(NA)
      } else{
        return(r$p.val)
      }
    }
  )
  # extract direction
  print("pearson stat")
  pearson.statistic <-
    sapply(
      psr_cor_tests ,
      FUN = function(r) {
        if (is.na(r)) {
          return(NA)
        } else{
          return(r$statistic)
        }
      }
    )
  
  
  stopCluster(clus)
  gc()
  df_x <-
    data.frame(
      pearson.p_val = pearson.p_val,
      pearson.statistic = pearson.statistic,
      pearson.delta = pearson.delta
    )
  return(df_x)
}


#' score stores number of data points out of all typisation$sample samples
#' function is parallized
#' @param df data frame with all methlyation values
#' with the appropriate sample column names.
#' @param typisation dataframe storing sample names in $sample
#' @returns numeric score vector
parallel_score_sample_number <- function(df, typisation) {
  clus <- parallel::makeCluster(parallel::detectCores() - 1)
  score  <-
    parRapply(
      cl = clus,
      x = df[, as.character(typisation$sample)],
      FUN = function(r) {
        sum(!is.na(r))
      }
    )
  stopCluster(clus)
  gc()
  return(score)
}

#' create_horizontal_permutated_typisation creates a random
#' permutation of given typisation
#' @param typisation with the appropriate sample column name in $sample
#' @returns list with permuted typisation sim_typisation and
#' numeric vector with permutation order
create_horizontal_permutated_typisation <- function(typisation) {
  permuation_order <- sample(nrow(typisation))
  sim_typisation <-  typisation
  sim_typisation$sample <- sim_typisation$sample[permuation_order]
  return(list(sim_typisation = sim_typisation, permuation_order = permuation_order))
}
# list.files(path =  file.path(OUTPUT_FOLDER, "simulation"))



#' counts number of significant values
#' for given p val and delta requirements
#' @param df with the  column name $pval
#' and $delta
#' @param minus_log_alpha define minimum pval = exp(-minus_log_alpha)
#' @param min_delta define minimum delta to be reached to be counted
#' @returns number of signifincat pvals
get_n_significant  <-
  function(df, minus_log_alpha, min_delta) {
    return(sum(
      df$p_val <= exp(-minus_log_alpha) &
        df$delta >= min_delta,
      na.rm = TRUE
    ))
  }

#' counts number of significant values with positive statistic
#' for given p val and delta requirements
#' @param df with the  column name $pval
#' and $delta
#' @param minus_log_alpha define minimum pval = exp(-minus_log_alpha)
#' @param min_delta define minimum delta to be reached to be counted
#' @returns number of signifincat pvals
get_positive_significant  <-
  function(df, minus_log_alpha, min_delta) {
    return(sum(
      df$p_val <= exp(-minus_log_alpha) &
        df$delta >= min_delta &
        df$statistic > 0,
      na.rm = TRUE
    ))
  }


#' counts number of significant values with negative statistic
#' for given p val and delta requirements
#' @param df with the  column name $pval
#' and $delta
#' @param minus_log_alpha define minimum pval = exp(-minus_log_alpha)
#' @param min_delta define minimum delta to be reached to be counted
#' @returns number of signifincat pvals
get_negative_significant  <-
  function(df, minus_log_alpha, min_delta) {
    return(sum(
      df$p_val <= exp(-minus_log_alpha) &
        df$delta >= min_delta &
        df$statistic < 0,
      na.rm = TRUE
    ))
  }

#' counts significant pval with different alpha and delta cut off
#' for given p val = exp(-minus_log_alpha) and delta
#' @param sim_file_names with the  column name $test_type.pval
#' and $test_type.delta
#' @param minus_log_alphas define minimum pval = exp(-minus_log_alpha)
#' @param min_deltas define minimum delta to be reached to be counted
#' @returns dataframe with number of signifincat pvals for each combination
#' in each file
parallel_summerize_permutations <-
  function(sim_file_names,
           min_deltas = seq(0, 0.50, 0.01),
           minus_log_alphas = seq(from = 0, to = 10, by = 0.1)) {
    sim_results_colnames <-
      c(
        "test",
        "permuation_type" ,
        "min_delta",
        "minus_log_alpha",
        "n_signfincant_CpG",
        "negative_signfincant_CpG",
        "positive_signfincant_CpG",
        "name"
      )
    
    #create and register cluster
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(n.cores)
    clusterExport(my.cluster, "get_n_significant")
    clusterExport(my.cluster, "get_negative_significant")
    clusterExport(my.cluster, "get_positive_significant")
    doParallel::registerDoParallel(cl = my.cluster)
    
    sim_results <- foreach(n = 1:length(sim_file_names),
                           .combine = 'rbind') %dopar% {
                             #load data
                             sim_listfile <-
                               readRDS(file = sim_file_names[n])
                             sim <- sim_listfile$data
                             
                             sim_file_name <-
                               basename(sim_file_names[n])
                             
                             
                             # allocated df
                             temp_sim <-
                               data.frame(matrix(
                                 nrow = length(min_deltas) * length(minus_log_alphas),
                                 ncol = length(sim_results_colnames)
                               ))
                             colnames(temp_sim) <-
                               sim_results_colnames
                             
                             temp_sim$name <-
                               gsub(x = sim_file_name,
                                    pattern = ".rds",
                                    replacement = "")
                             temp_sim$test <-
                               gsub(pattern = "\\..*",
                                    replacement = "",
                                    x = sim_file_name)
                             sub_string <-
                               gsub(
                                 pattern = paste(temp_sim$test[1], ".|.rds", sep = ""),
                                 replacement = "",
                                 x = sim_file_name
                               )
                             temp_sim$permuation_type <-
                               gsub(pattern = "\\..*",
                                    replacement = "",
                                    x = sub_string)
                             temp_sim$minus_log_alpha <-
                               rep(
                                 x = minus_log_alphas,
                                 times = length(min_deltas),
                                 length.out = NA,
                                 each = 1
                               )
                             temp_sim$min_delta <-
                               rep(
                                 x = min_deltas,
                                 times = 1,
                                 length.out = NA,
                                 each = length(minus_log_alphas)
                               )
                             
                             colnames(sim) <-
                               c("p_val", "statistic", "delta")
                             
                             # set pval to 1 if no changes in methylation
                             sim$p_val[sim$delta == 0] <- 1
                             
                             temp_sim$n_signfincant_CpG <-
                               mapply(
                                 minus_log_alpha = temp_sim$minus_log_alpha,
                                 min_delta = temp_sim$min_delta,
                                 FUN = get_n_significant,
                                 MoreArgs = list(df = sim)
                               )
                             
                             temp_sim$negative_signfincant_CpG <-
                               mapply(
                                 minus_log_alpha = temp_sim$minus_log_alpha,
                                 min_delta = temp_sim$min_delta,
                                 FUN = get_negative_significant,
                                 MoreArgs = list(df = sim)
                               )
                             
                             temp_sim$positive_signfincant_CpG <-
                               mapply(
                                 minus_log_alpha = temp_sim$minus_log_alpha,
                                 min_delta = temp_sim$min_delta,
                                 FUN = get_positive_significant,
                                 MoreArgs = list(df = sim)
                               )
                             
                             # append results
                             return(temp_sim)
                           }
    parallel::stopCluster(cl = my.cluster)
    
    return(sim_results)
  }


#' counts CpGs passing delta cut off
#' @param sim_file_names with the  column name
#' @param min_deltas define minimum delta to be reached to be counted
#' @returns dataframe with number of CpGs
parallel_count_CpGs_with_delta_cutoff <-
  function(sim_file_names,
           min_deltas = seq(0, 0.50, 0.01)) {
    sim_results_colnames <-
      c("test",
        "permuation_type" ,
        "min_delta",
        "n_CpG",
        "name")
    
    #create and register cluster
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(n.cores)
    #clusterExport(my.cluster, "min_deltas")
    doParallel::registerDoParallel(cl = my.cluster)
    
    sim_results <- foreach(n = 1:length(sim_file_names),
                           .combine = 'rbind') %dopar% {
                             #load data
                             sim_listfile <-
                               readRDS(file = sim_file_names[n])
                             sim <- sim_listfile$data
                             
                             sim_file_name <-
                               basename(sim_file_names[n])
                             
                             # allocated df
                             temp_sim <-
                               data.frame(matrix(
                                 nrow = length(min_deltas),
                                 ncol = length(sim_results_colnames)
                               ))
                             colnames(temp_sim) <-
                               sim_results_colnames
                             
                             temp_sim$name <-
                               gsub(x = sim_file_name,
                                    pattern = ".rds",
                                    replacement = "")
                             temp_sim$test <-
                               gsub(pattern = "\\..*",
                                    replacement = "",
                                    x = sim_file_name)
                             sub_string <-
                               gsub(
                                 pattern = paste(temp_sim$test[1], ".|.rds", sep = ""),
                                 replacement = "",
                                 x = sim_file_name
                               )
                             temp_sim$permuation_type <-
                               gsub(pattern = "\\..*",
                                    replacement = "",
                                    x = sub_string)
                             temp_sim$min_delta <-
                               min_deltas
                             
                             colnames(sim) <-
                               c("p_val", "statistic", "delta")
                             
                             temp_sim$n_CpG <-
                               mapply(
                                 min_delta = temp_sim$min_delta,
                                 FUN = function(df, min_delta) {
                                   return(sum(df$delta >= min_delta))
                                 },
                                 MoreArgs = list(df = sim)
                               )
                             
                             # append results
                             return(temp_sim)
                           }
    parallel::stopCluster(cl = my.cluster)
    return(sim_results)
  }


## Convenience function
#https://stackoverflow.com/questions/7740503/getting-frequency-values-from-histogram-in-r
get_hist <- function(p) {
  d <- ggplot2::ggplot_build(p)$data[[1]]
  data.frame(
    x = d$x,
    xmin = d$xmin,
    xmax = d$xmax,
    y = d$y
  )
}


#' counts CpGs passing delta cut off
#' @param sim_file_names with the  column name
#' @param breaks = seq(from = 0, to = 1, by = 0.01))
#' @returns dataframe with number of CpGs in each bin
parallel_count_CpGs_with_delta_histogramm <-
  function(sim_file_names,
           breaks = seq(from = 0, to = 1, by = 0.01)) {
    #create and register cluster
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(n.cores)
    #clusterExport(my.cluster, "breaks")
    clusterExport(my.cluster, "get_hist")
    doParallel::registerDoParallel(cl = my.cluster)
    
    meth_counts <- foreach(n = 1:length(sim_file_names),
                           .combine = 'rbind') %dopar% {
                             #load data
                             sim_listfile <-
                               readRDS(file = sim_file_names[n])
                             sim <- sim_listfile$data
                             
                             sim_file_name <-
                               basename(sim_file_names[n])
                             
                             colnames(sim) <-
                               c("p_val", "statistic", "delta")
                             
                             p <- ggplot2::ggplot(data = sim,
                                                  mapping = ggplot2::aes(x = delta)) +
                               ggplot2::geom_histogram(breaks = breaks)
                             
                             df <- get_hist(p)
                             df$test <-
                               gsub(pattern = "\\..*",
                                    replacement = "",
                                    x = sim_file_name)
                             df$y_fraction <- df$y / sum(df$y)
                             
                             df$name <- sim_file_name
                             
                             # append results
                             return(df)
                           }
    parallel::stopCluster(cl = my.cluster)
    return(meth_counts)
  }

#' create CpG permuartion for each column in dataframe
#' @param df with the  column name df$sample and
#' methylation values that will be permuted inseide
#' @param typisation define permuted columns
#' @returns list with permuted dataframe $data and
#' list of permuation_order in each column in $permuation_order
create_CpG_permution <- function(df, typisation) {
  permuation_order <- list()
  
  n_CpG <- nrow(df)
  for (i in 1:nrow(typisation)) {
    permuation <- sample(1:n_CpG)
    df[, as.character(typisation$sample[i])] <-
      df[permuation, as.character(typisation$sample[i])]
    permuation_order[[as.character(typisation$sample[i])]] <-
      permuation
  }
  df <- df[,!colnames(df) %in%
             c(pearson_test_colnames, KW_test_colnames)]
  return(list(data = df, permuation_order = permuation_order))
}

#' create a linear correlation plot meth vs age from row  in dataframe row
#' @param df_row is a row in the CpG dataframe with
#' columnnames like the typisation$sample
#' @param typisation age typisation (can be permuted for showing sillulated values)
#' @returns plot with linear correlation line meth cs age
plot_age_correlation <- function(df_row, typisation) {
  df_plot <- typisation
  df_plot$meth <-
    as.numeric(df_row[, as.character(typisation$sample)])
  
  p <- ggplot(data = df_plot,
              mapping = aes(x = age_mean_BP, y = meth)) +
    ggtitle(
      label = paste(df_row$chrom, df_row$start),
      subtitle = paste("p_val: ", format.pval(df_row$pearson.p_val))
    ) +
    geom_point() +
    ylim(0, 1) +
    xlim(0, 11500) +
    xlab("sample age estimation in Years before Present") +
    ylab("CpG methylation") +
    geom_smooth(method = "lm", ) +
    theme_minimal()
  
  
  return(p)
}


#' create a KW test boxplot from row  in dataframe row
#' @param df_row is a row in the CpG dataframe with
#' columnnames like the typisation$sample
#' @param typisation food typisation (can be permuted for simulated values)
#' @returns boxplot with KW pvalue and gemonic position in title
plot_food_correlation <- function(df_row, typisation) {
  df_plot <- typisation
  df_plot$meth <-
    as.numeric(df_row[, as.character(typisation$sample)])
  
  p <- ggplot(data = df_plot,
              mapping = aes(x = Type, y = meth, color = Type)) +
    ggtitle(
      label = paste(df_row$chrom, df_row$start),
      subtitle = paste("p_val: ", format.pval(df_row$kw.p_val))
    ) +
    geom_boxplot() +
    geom_point(position = "jitter") +
    ylim(0, 1) +
    xlab("sample food Type") +
    ylab("CpG methylation") +
    theme_minimal()
  return(p)
}


# H3K27ac Analysis --------------------------------------------------------

if (H3K27ac_Analysis) {
  # load entire methylome
  all_CpG.39.samples.merged.hg19 <-
    readRDS(
      "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/all_CpG.39.samples.merged.hg19.rds"
    )
  
  # create results folder
  folder_name <- file.path(OUTPUT_FOLDER, "results", "H3K27ac")
  if (!file.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }
  
  all_CpG.meth <-
    all_CpG.39.samples.merged.hg19[, as.character(meta$sample)]
  
  percent_data <-
    apply(
      X =  all_CpG.meth ,
      MARGIN = 2,
      FUN = function(meth) {
        sum(!is.na(meth)) / length(meth)
      }
    )
  meta$methylation_coverage <- percent_data
  print(paste(
    "min coverage ",
    min(meta$methylation_coverage),
    "max coverage ",
    max(meta$methylation_coverage)
  ))
  
  for (n_chr in 1:length(CHR_NAMES)) {
    print(n_chr)
    file_name <- file.path(
      OUTPUT_FOLDER,
      "results",
      "H3K27ac",
      paste(CHR_NAMES[n_chr], "CpG_annotation", "rds", sep = ".")
    )
    
    if (file.exists(file_name)) {
      next
    }
    # cut out only coordinates per chromosome
    # n_chr <- 22
    df_chr <-
      all_CpG.39.samples.merged.hg19[all_CpG.39.samples.merged.hg19$chrom == CHR_NAMES[n_chr], c(1:6)]
    # load state annotation
    chr_chromatin_seg <-
      readRDS(
        paste(
          "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/Chromatin/processed/",
          paste(CHR_NAMES[n_chr], "chromatin_seg_annotaion_name", "RDS", sep = "."),
          sep = ""
        )
      )
    
    # pefrom mapping of CpGs to state
    start_time <- Sys.time()
    df_chr$state <- NA
    n_seg <- 1
    seg_end <-  chr_chromatin_seg$end[n_seg]
    pb = txtProgressBar(min = 1,
                        max = nrow(df_chr),
                        initial = 1)
    
    for (r in 1:nrow(df_chr)) {
      cpg_pos <- df_chr$end[r]
      while (cpg_pos > seg_end) {
        n_seg <- n_seg + 1
        seg_end <-  chr_chromatin_seg$end[n_seg]
      }
      df_chr$state[r] <- chr_chromatin_seg$state[n_seg]
      setTxtProgressBar(pb, r)
    }
    close(pb)
    end_time <- Sys.time()
    print(end_time - start_time)
    
    saveRDS(object = df_chr,
            file = file_name)
    
    ggplot(data = df_chr, mapping = aes(x = state, fill = state)) +
      geom_histogram(stat = "count") +
      ggtitle(paste(ANNOTATION, "genome wide Chromatin state destribution", sep = " ")) +
      ggplot2::theme(
        axis.line = ggplot2::element_line(colour = "black"),
        panel.background = ggplot2::element_blank(),
        legend.position = "none",
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      )
  }
  
  
  all_CpG_states_hg19 <- data.frame()
  for (n_chr in 1:length(CHR_NAMES)) {
    print(n_chr)
    file_name <- file.path(
      OUTPUT_FOLDER,
      "results",
      "H3K27ac",
      paste(CHR_NAMES[n_chr], "CpG_annotation", "rds", sep = ".")
    )
    df_chr <- readRDS(file = file_name)
    all_CpG_states_hg19 <- rbind(all_CpG_states_hg19, df_chr)
  }
  
  matching_indices <-
    match(
      all_CpG_states_hg19$state,
      paste(
        df_universal_annotation$state_order_by_group,
        df_universal_annotation$`States  presented in paper`,
        sep = "_"
      )
    )
  all_CpG_states_hg19$state_color <-
    factor(df_universal_annotation$color[matching_indices],
           levels = unique(df_universal_annotation$color))
  all_CpG_states_hg19$state_group <-
    factor(df_state_summery$state_group[matching_indices],
           levels = unique(df_state_summery$state_group))
  all_CpG_states_hg19$state <-
    factor(x = all_CpG_states_hg19$state, levels = state_names)
  
  rm(df_chr)
  rm(matching_indices)
  gc()
  
  saveRDS(
    object = all_CpG_states_hg19,
    file =
      file_name <- file.path(
        OUTPUT_FOLDER,
        "results",
        "H3K27ac",
        "all_CpG_states_hg19.rds"
      )
  )
  
  
  p <-
    ggplot(data = all_CpG_states_hg19,
           mapping = aes(x = state, fill = state_color)) +
    geom_histogram(stat = "count", color = "black") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank(),
      #legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    ) +
    ylab("# of CpGs in chromatin state") +
    xlab("chromatin state labeling") +
    scale_fill_manual(
      values = unique(df_universal_annotation$color),
      labels = unique(df_universal_annotation$Group)
    )
  
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "results",
                         "H3K27ac",
                         "all_CpG_states.png")
    ,
    plot = p,
    width = 10,
    height = 6
  )
  
  
  
  # all_chromatin_seg
  all_chromatin_seg <- data.frame()
  for (n_chr in 1:length(CHR_NAMES)) {
    print(n_chr)
    file_name <- file.path(
      OUTPUT_FOLDER,
      "results",
      "H3K27ac",
      paste(CHR_NAMES[n_chr], "CpG_annotation", "rds", sep = ".")
    )
    
    # load state annotation
    chr_chromatin_seg <-
      readRDS(
        paste(
          "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/Chromatin/processed/",
          paste(CHR_NAMES[n_chr], "chromatin_seg_annotaion_name", "RDS", sep = "."),
          sep = ""
        )
      )
    all_chromatin_seg <- rbind(all_chromatin_seg, chr_chromatin_seg)
  }
  
  
  
  all_chromatin_seg$length <-
    all_chromatin_seg$end - all_chromatin_seg$start
  
  # compute coverage
  coverage <-
    aggregate(
      x = all_chromatin_seg$length,
      by = list(Category = all_chromatin_seg$state),
      FUN = sum
    )
  # compute total CpG
  n_CpG <-
    aggregate(
      all_CpG_states_hg19$strand,
      by = list(all_CpG_states_hg19$state),
      FUN = length
    )
  # compute CpG in H3K27ac
  # H3K27ac_CpG <-
  #   aggregate(
  #     all_CpG_states_hg19$name[all_CpG_states_hg19$name != "NO_CHIP"],
  #     by = list(all_CpG_states_hg19$state[all_CpG_states_hg19$name != "NO_CHIP"]),
  #     FUN = length
  #   )
  
  H3K27ac_CpG <-
    aggregate(
      all_CpG_states_hg19$name != "NO_CHIP" ,
      by = list(
        all_CpG_states_hg19$state,
        all_CpG_states_hg19$name != "NO_CHIP"
      ),
      FUN = length
    )
  
  not_H3K27ac_CpG <-
    H3K27ac_CpG[!H3K27ac_CpG$Group.2, c("Group.1", "x")]
  
  state_segementation_meta <-
    data.frame(
      state = coverage$Category,
      coverage = coverage$x,
      n_CpG = n_CpG$x,
      not_H3K27ac_CpG = not_H3K27ac_CpG$x
    )
  state_segementation_meta$H3K27ac_CpG <-
    state_segementation_meta$n_CpG - state_segementation_meta$not_H3K27ac_CpG
  
  state_segementation_meta$state_fraction_is_H3K27ac <-
    state_segementation_meta$H3K27ac_CpG / state_segementation_meta$n_CpG
  
  state_segementation_meta$state <-
    factor(x = state_segementation_meta$state,
           levels = state_segementation_meta$state[order(as.numeric(
             gsub('_.*', replacement = "", x = state_segementation_meta$state)
           ))])
  
  
  # p <-
  #   ggplot(
  #     data = state_segementation_meta,
  #     mapping = aes(x = state, y = state_fraction_is_H3K27ac * 100, fill = state)
  #   ) +
  #   geom_bar(stat = "identity") +
  #   ggplot2::theme(
  #     axis.line = ggplot2::element_line(colour = "black"),
  #     panel.background = ggplot2::element_blank(),
  #     legend.position = "none",
  #     axis.text.x = element_text(
  #       angle = 90,
  #       vjust = 0.5,
  #       hjust = 1
  #     )
  #   ) +
  #   ylim(0, 100) +
  #   ylab("% of CpGs covered by H3K27ac peaks in bone sample") +
  #   xlab("genome segment chromatin state labeling")
  
  p <- ggplot(data = df_summery_by_state,mapping = aes(x = state,y = n_H3K27AC/n_total*100, fill = state_group))+
           geom_bar(stat = "identity",color = "black")+
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    )+
    ylim(0, 100) +
    ylab("% of CpGs in genome covered by H3K27ac") +
    xlab("state group")+
    scale_fill_manual(values = unique(as.character(df_summery_by_state$state_color)))
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "results",
                         "H3K27ac",
                         "H3K27ac_coverage.png")
    ,
    plot = p,
    width = 10,
    height = 6
  )
  
  print("# of CpG in Genome: ")
  n_CpG_genome <- sum(state_segementation_meta$n_CpG)
  print(n_CpG_genome)
  print("# of CpG in H3K27Ac: ")
  n_CpG_H3K27ac <- sum(state_segementation_meta$H3K27ac_CpG)
  print(n_CpG_H3K27ac)
  print("fraction h3K27AC:")
  print(n_CpG_H3K27ac / n_CpG_genome)
}


# Genome wide Analysis --------------------------------------------------------
if (Genome_wide_Analysis) {
  # load entire methylome
  all_CpG.39.samples.merged.hg19 <-
    readRDS(
      "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/all_CpG.39.samples.merged.hg19.rds"
    )
  
  # create results folder
  folder_name <- file.path(OUTPUT_FOLDER, "results", "WholeGenome")
  if (!file.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }
  
  
  for (n_chr in 1:length(CHR_NAMES)) {
    # n_chr <- 22
    print(n_chr)
    
    file_name <- file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste(CHR_NAMES[n_chr], "_complete_with_test", "rds", sep = ".")
    )
    if (file.exists(file_name)) {
      next
    }
    
    # cut out only coordinates per chromosome
    df_chr <-
      all_CpG.39.samples.merged.hg19[all_CpG.39.samples.merged.hg19$chrom == CHR_NAMES[n_chr],]
    
    # compute sample numbers per CpG position
    start_time <- Sys.time()
    df_chr$score <-
      parallel_score_sample_number(df_chr, real_age_typisation)
    end_time <- Sys.time()
    print(end_time - start_time)
    print(paste("# CpG in chr:", nrow(df_chr)))
    print(paste("# CpG in chr with all 39 samples :",
                sum(
                  df_chr$score == nrow(real_age_typisation)
                )))
    print(paste("discard ", 100 * (
      1 - sum(df_chr$score == nrow(real_age_typisation)) / nrow(df_chr)
    ), "% of CpG positions"))
    
    
    # discard CpG with not full data
    df_chr_complete <-
      df_chr[nrow(real_age_typisation) == df_chr$score,]
    
    # compute real data value kruskall_valis
    start_time <- Sys.time()
    df_kruskall_valis <-
      parallel_testing_kruskall_valis(df_chr_complete, real_food_typisation)
    end_time <- Sys.time()
    print(end_time - start_time)
    
    # compute real data value pearson_cor
    start_time <- Sys.time()
    df_testing_pearson_cor <-
      parallel_testing_pearson_cor(df = df_chr_complete, age.typisation = real_age_typisation)
    end_time <- Sys.time()
    print(end_time - start_time)
    
    
    # save results
    df_chr_complete_with_test <-
      cbind(df_chr_complete,
            df_testing_pearson_cor,
            df_kruskall_valis)
    
    annotation_file_name <- file.path(
      OUTPUT_FOLDER,
      "results",
      "H3K27ac",
      paste(CHR_NAMES[n_chr], "CpG_annotation", "rds", sep = ".")
    )
    chr_annotation <-  readRDS(file = annotation_file_name)
    
    df_chr_complete_with_test$state <-
      chr_annotation$state[nrow(real_age_typisation) == df_chr$score]
    
    saveRDS(object = df_chr_complete_with_test,
            file = file_name)
    
  }
  
}
# all CpGs dataframe -all_CpG_data_hg19 -all_data_complete_with_test -----------------------------------------------------------------------
if (FALSE) {
  # compine all chr data to big df
  all_CpG_data_hg19 <- data.frame()
  for (n_chr in 1:length(CHR_NAMES)) {
    print(n_chr)
    file_name <- file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste(CHR_NAMES[n_chr], "_complete_with_test", "rds", sep = ".")
    )
    df_chr <- readRDS(file = file_name)
    all_CpG_data_hg19 <- rbind(all_CpG_data_hg19, df_chr)
  }
  
  state_names <- unique(as.character(all_CpG_data_hg19$state))
  state_names <- state_names[order(pad_with_zeros(state_names))]
  
  all_CpG_data_hg19$state <- factor(x = all_CpG_data_hg19$state,
                                    levels = state_names)
  
  matching_indices <-
    match(
      all_CpG_data_hg19$state,
      paste(
        df_universal_annotation$state_order_by_group,
        df_universal_annotation$`States  presented in paper`,
        sep = "_"
      )
    )
  all_CpG_data_hg19$state_color <-
    factor(df_universal_annotation$color[matching_indices],
           levels = unique(df_universal_annotation$color))
  all_CpG_data_hg19$state_group <-
    factor(df_state_summery$state_group[matching_indices],
           levels = unique(df_state_summery$state_group))
  
  saveRDS(
    object = all_CpG_data_hg19,
    file =  file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("all_data_complete_with_test", "rds", sep = ".")
    )
  )
  
  # all_CpG_data_hg19 <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/results/WholeGenome/all_data_complete_with_test.rds")
  # state_names <- unique(as.character(all_CpG_data_hg19$state))
  # state_names <- state_names[order(pad_with_zeros(state_names))]
  # df_state_summery <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/results/WholeGenome/df_state_summery.rds")
  # df_summery_by_state <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/results/WholeGenome/df_summery_by_state.rds")
  
  # # Convert your dataframe to a data.table
  #
  #
  # # Melt the dataframe
  # all_CpG_meth <- data.table::melt(data = meth, id.vars = NULL, variable.name = "sample", value.name = "methylation")
  #
  
  
  
  df_plot_no_arf <-
    all_CpG_data_hg19[!grepl(pattern = "GapArtf",
                             x = all_CpG_data_hg19$state,
                             ignore.case = TRUE),]
  print(dim(df_plot_no_arf))
  
  # plot delta
  p <-
    ggplot(
      data = df_plot_no_arf,
      mapping = aes(x = state_group, y = pearson.delta, fill = state_group)
    ) +
    geom_violin() +
    ggplot2::theme(#axis.line = ggplot2::element_line(colour = "black"),
      #panel.background = ggplot2::element_blank(),
      #legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
    ylab(bquote(""  ~ delta[meth])) +
    xlab("chromatin state labeling") +
    scale_fill_manual(values = unique(df_universal_annotation$color)[-1])
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("pearson.delta_WG_complete_with_test", "png", sep = ".")
    )
    ,
    plot = p,
    width = 18,
    height = 6
  )
  
  
  # plot -log p
  p <-
    ggplot(data = df_plot_no_arf,
           mapping = aes(
             x = state_group,
             y = -log(pearson.p_val),
             fill = state_group
           )) +
    geom_violin() +
    ggplot2::theme(
      #axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    ) +
    ylab(bquote(""  ~ -log(p_[value]))) +
    xlab("chromatin state labeling") +
    ylim(0, 5) +
    scale_fill_manual(values = unique(df_universal_annotation$color)[-1])
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("pearson_log_pvalue_WG_complete_with_test", "png", sep = ".")
    )
    ,
    plot =  p,
    width = 18,
    height = 6
  )
  
  p <-
    ggplot(data = df_plot_no_arf,
           mapping = aes(
             x = state_group,
             y = -log(pearson.p_val),
             fill = state_group
           )) +
    geom_violin() +
    ggplot2::theme(
      #axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    ) +
    ylab(bquote(""  ~ -log(p_[value]))) +
    xlab("chromatin state labeling") +
    ylim(0, 5) +
    scale_fill_manual(values = unique(df_universal_annotation$color)[-1])
  
  # ggsave(filename = file.path(
  #   OUTPUT_FOLDER,
  #   "results",
  #   "WholeGenome",
  #   paste("methylation_WG_complete_with_test", "png", sep = "."))
  #   ,plot = p,width = 36,height = 12)
  
  #df_state_summery <- summary_df <- aggregate(. ~ state, data =  all_CpG_data_hg19[,c("state",as.character(meta$sample))], FUN = mean)
  state_summery_column_names <-
    c(
      "state",
      "sample",
      "mean_methylation",
      "sd_methylation",
      "methylation_H3K27AC",
      "methylation_NO_CHIP"
    )
  
  # Create an empty data frame with the specified column names
  df_state_summery <-
    data.frame(matrix(
      NA,
      nrow = length(state_names) * nrow(meta),
      ncol = length(state_summery_column_names)
    ))
  
  # Rename the columns
  colnames(df_state_summery) <- state_summery_column_names
  
  r <- 1
  for (state in  state_names) {
    print(state)
    for (sample in as.character(meta$sample)) {
      df_state_summery$state[r] <- state
      df_state_summery$sample[r] <- sample
      
      index <- all_CpG_data_hg19$state == state
      meth_column <- all_CpG_data_hg19[, sample]
      meth <- meth_column[index]
      
      df_state_summery$mean_methylation[r] <- mean(meth)
      df_state_summery$sd_methylation[r] <- sd(meth)
      
      df_state_summery$methylation_H3K27AC[r] <-
        mean(meth_column[index & all_CpG_data_hg19$name != "NO_CHIP"])
      df_state_summery$methylation_NO_CHIP[r] <-
        mean(meth_column[index & all_CpG_data_hg19$name == "NO_CHIP"])
      
      r <- r + 1
    }
  }
  
  # change state to factor
  df_state_summery$state <-
    factor(x = df_state_summery$state, levels = state_names)
  
  # annotate states
  matching_indices <-
    match(
      df_state_summery$state,
      paste(
        df_universal_annotation$state_order_by_group,
        df_universal_annotation$`States  presented in paper`,
        sep = "_"
      )
    )
  df_state_summery$state_color <-
    df_universal_annotation$color[matching_indices]
  df_state_summery$state_group <-
    factor(df_universal_annotation$Group[matching_indices],
           levels = unique(df_universal_annotation$Group))
  # annotate sample
  matching_indices <- match(df_state_summery$sample, meta$sample)
  df_state_summery$type <- meta$Type[matching_indices]
  
  
  saveRDS(
    object = df_state_summery,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("df_state_summery", "rds", sep = ".")
    )
  )
  
  
  df_plot <-
    df_state_summery[!grepl(pattern = "GapArtf",
                            x = df_state_summery$state,
                            ignore.case = TRUE),]
  
  fill_vec <- df_universal_annotation$color[(100-length(unique(df_plot$state))):length(unique(df_plot$state))]
  
  muh_grob <- grid::rectGrob(
    x=1:length(unique(df_plot$state)), y=0, gp=gpar(
      color='black', fill=rainbow(10) , alpha=0.2))
  
  p <-
    ggplot(data = df_plot ,
           mapping = aes(x = state, y = mean_methylation)) +
    #geom_boxplot(outlier.shape = NA)+
    geom_point(position = "jitter", mapping = aes(color = type)) +
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      panel.background = ggplot2::element_blank()
    ) +
    ylab("mean sample methylation") +
    xlab("chromatin state labeling") #+
    # annotation_custom(
    #   grob=muh_grob, xmin = 0, xmax = 1, ymin = -0.5, ymax=0.1
    # )
  
  
  
  # +
  #   scale_fill_manual(values = unique(df_universal_annotation$color)[-1])
  #
  ggsave(filename = file.path(
    OUTPUT_FOLDER,
    "results",
    "WholeGenome",
    paste("methylation_WG_complete_with_test", "png", sep = "."))
   ,plot = p,width = 24,height = 12)
  
  state_column_names <-
    c(
      "state",
      "mean_methylation",
      "sd_methylation",
      "mean_delta",
      "sd_delta",
      "mean_minus_log_p",
      "sd_minus_log_p",
      "n_total",
      "n_H3K27AC"
    )
  
  # Create an empty data frame with the specified column names
  df_summery_by_state <-
    data.frame(matrix(
      NA,
      nrow = length(state_names),
      ncol = length(state_column_names)
    ))
  
  # Rename the columns
  colnames(df_summery_by_state) <- state_column_names
  
  r <- 1
  for (state in  state_names) {
    print(state)
    
    df_summery_by_state$state[r] <- state
    
    index <- all_CpG_data_hg19$state == state
    meth <- all_CpG_data_hg19[index, as.character(meta$sample)]
    df_summery_by_state$mean_delta[r] <-
      mean(all_CpG_data_hg19$pearson.delta[index])
    df_summery_by_state$sd_delta[r] <-
      sd(all_CpG_data_hg19$pearson.delta[index])
    df_summery_by_state$mean_minus_log_p[r] <-
      mean(-log(all_CpG_data_hg19$pearson.p_val[index]), na.rm = TRUE)
    df_summery_by_state$sd_minus_log_p[r] <-
      sd(-log(all_CpG_data_hg19$pearson.p_val[index]), na.rm = TRUE)
    df_summery_by_state$mean_methylation[r] <-
      mean(as.matrix(meth))
    df_summery_by_state$sd_methylation[r] <-
      sd(as.matrix(meth))
    df_summery_by_state$n_total[r] <- sum(index)
    df_summery_by_state$n_H3K27AC[r] <-
      sum(index & all_CpG_data_hg19$name != "NO_CHIP")
    r <- r + 1
  }
  
  n_all_CpG <- sum(df_summery_by_state$n_total)
  n_H3K27AC_CpG <- sum(df_summery_by_state$n_H3K27AC)
  expected_fraction <- n_H3K27AC_CpG / n_all_CpG
  
  df_summery_by_state$expected <-
    df_summery_by_state$n_total * expected_fraction
  
  df_summery_by_state$state <-
    factor(x = df_summery_by_state$state,
           levels = state_names)
  
  
  matching_indices <-
    match(
      df_summery_by_state$state,
      paste(
        df_universal_annotation$state_order_by_group,
        df_universal_annotation$`States  presented in paper`,
        sep = "_"
      )
    )
  df_summery_by_state$state <-
    factor(df_summery_by_state$state, levels = df_summery_by_state$state)
  df_summery_by_state$state_color <-
    factor(df_universal_annotation$color[matching_indices],
           levels = unique(df_universal_annotation$color))
  df_summery_by_state$state_group <-
    factor(df_state_summery$state_group[matching_indices],
           levels = unique(df_state_summery$state_group))
  
  df_summery_by_state$p_val_test_all_methylation_KW <- NA
  
  for (st in state_names) {
    print(st)
    meth <-
      data.table::setDT(all_CpG_data_hg19[all_CpG_data_hg19$state == st, as.character(meta$sample)])
    meth <-
      data.table::melt(
        data = meth,
        id.vars = NULL,
        variable.name = "sample",
        value.name = "methylation"
      )
    # # Perform Kruskal-Wallis test
    kw_state <- kruskal.test(methylation ~ sample, data = meth)
    df_summery_by_state$p_val_test_all_methylation_KW[df_summery_by_state$state == st] <-
      format.pval(kw_state$p.value, digits = 3)
  }
  
  saveRDS(
    object = df_summery_by_state,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("df_summery_by_state", "rds", sep = ".")
    )
  )
  
  df_plot_summery_by_state <-
    df_summery_by_state[!grepl(pattern = "GapArtf", x = df_summery_by_state$state), ]
  
  p <-
    ggplot(
      data = df_plot_summery_by_state,
      mapping = aes(x = state, y = mean_minus_log_p, fill = state_group)
    ) +
    geom_bar(stat = "identity", color = "black") +
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      panel.background = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ylab("-log_p") +
    xlab("chromatin state labeling") +
    scale_fill_manual(values = unique(df_universal_annotation$color)[-1])
  
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("p_p_val_WG_complete_with_test", "png", sep = ".")
    )
    ,
    plot = p,
    width = 36,
    height = 12
  )
  
  p_enrichment <-
    ggplot(data = df_summery_by_state,
           mapping = aes(
             x = state,
             y = n_H3K27AC / expected,
             fill = state
           )) +
    geom_bar(stat = "identity") +
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      panel.background = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ylab("enrichment: oberserved / expected") +
    xlab("chromatin state labeling")
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("p_enrichment_WG_complete_with_test", "png", sep = ".")
    )
    ,
    plot = p_enrichment,
    width = 36,
    height = 12
  )
  
  df_summery_by_state$p_value_NO_CHIP <- NA
  df_summery_by_state$p_value_CHIP <- NA
  
  for (s in 1:length(state_names)) {
    print(s)
    data <-
      df_state_summery[state_names[s] == df_state_summery$state,]
    data$age <- meta$age_mean_BP
    
    # Perform Pearson correlation test
    cor_test <-
      cor.test(data$age, data$mean_methylation, method = "pearson")
    
    # Extract the p-value from the test results
    p_value <- cor_test$p.value
    
    # Plot scatterplot with trend line
    p <- ggplot(data, aes(x = age, y = mean_methylation)) +
      geom_point() +  # Add scatterplot points
      geom_smooth(method = "lm", se = FALSE) +  # Add trend line (linear regression)
      labs(x = "Age", y = "Mean Methylation") +  # Add axis labels
      ggtitle(paste(state_names[s], " p =", sprintf("%.2f", p_value))) +  # Add title
      theme_minimal()  # Use minimal theme (optional, customize as needed)
    
    ggsave(
      filename = file.path(
        OUTPUT_FOLDER,
        "results",
        "WholeGenome",
        "plots",
        paste("methylation", state_names[s], "png", sep = ".")
      )
      ,
      plot = p,
      width = 6,
      height = 6
    )
    
    
    # Perform Pearson correlation test
    cor_test_NO_CHIP <-
      cor.test(data$age, data$methylation_NO_CHIP, method = "pearson")
    
    
    # Extract the p-value from the test results
    p_value_NO_CHIP <- cor_test_NO_CHIP$p.value
    
    df_summery_by_state$p_value_NO_CHIP[as.character(df_summery_by_state$state) == state_names[s]] <-
      p_value_NO_CHIP
    
    # compare H3K37Ac vs no chip
    n <-
      df_summery_by_state$n_H3K27AC[as.character(df_summery_by_state$state) == state_names[s]]
    
    if (n > 3) {
      # Perform Pearson correlation test
      cor_test_CHIP <-
        cor.test(data$age, data$methylation_H3K27AC, method = "pearson")
      # Extract the p-value from the test results
      p_value_CHIP <- cor_test_CHIP$p.value
      df_summery_by_state$p_value_CHIP[as.character(df_summery_by_state$state) == state_names[s]] <-
        p_value_CHIP
    }
    
    
    
    n1 <-
      df_summery_by_state$n_total[as.character(df_summery_by_state$state) == state_names[s]] - df_summery_by_state$n_H3K27AC[as.character(df_summery_by_state$state) == state_names[s]]
    
    # Plot scatterplot with trend line
    p1 <- ggplot(data, aes(x = age, y = methylation_NO_CHIP)) +
      geom_point() +  # Add scatterplot points
      geom_smooth(method = "lm", se = FALSE) +  # Add trend line (linear regression)
      labs(x = "Age", y = "Mean Methylation NO_CHIP") +  # Add axis labels
      ggtitle(paste(
        "n=",
        n1,
        state_names[s],
        " p =",
        sprintf("%.2f", p_value_NO_CHIP)
      )) +  # Add title
      theme_minimal()  # Use minimal theme (optional, customize as needed)
    
    
    if (n > 3) {
      # Plot scatterplot with trend line
      p2 <- ggplot(data, aes(x = age, y = methylation_H3K27AC)) +
        geom_point() +  # Add scatterplot points
        geom_smooth(method = "lm", se = FALSE) +  # Add trend line (linear regression)
        labs(x = "Age", y = "Mean Methylation H3K27Ac") +  # Add axis labels
        ggtitle(paste(
          "n=",
          n ,
          state_names[s],
          " p =",
          sprintf("%.2f", p_value_CHIP)
        )) +  # Add title
        theme_minimal()  # Use minimal theme (optional, customize as needed)
    } else{
      p2 <- ggplot()
    }
    
    
    
    # Combine the plots vertically
    combined_plot <- plot_grid(p1, p2, ncol = 1)
    
    ggsave(
      filename = file.path(
        OUTPUT_FOLDER,
        "results",
        "WholeGenome",
        "plots",
        paste("comparison.methylation", state_names[s], "png", sep = ".")
      )
      ,
      plot = combined_plot,
      width = 6,
      height = 12
    )
    
    
  }
  
  
  big_labels <- as.character(df_summery_by_state$state)
  big_labels[!as.character(df_summery_by_state$state) %in% head(as.character(df_summery_by_state$state)[order(df_summery_by_state$n_H3K27AC, decreasing = TRUE)], 10)] <-
    " "
  df_summery_by_state$big_labels <- big_labels
  
  # Create the pie chart with borders
  p<-
    ggplot(data = df_summery_by_state, aes(
      x = 1,
      y = n_H3K27AC,
      fill = state,
      label = big_labels
    )) +
    coord_polar(theta = "y", start = 0) +  # Convert the bar plot into a pie chart
    labs(fill = "State", x = NULL) +   # Legend title
    theme_void() +  # Remove unnecessary elements
    geomtextpath::geom_textpath(
      data = df_summery_by_state,
      position = position_stack(vjust = 0.5),
      aes(
        x = 1.75,
      ),
      size = 3.4,
      show.legend = FALSE
    ) +
    geom_bar(width = 1,
             stat = "identity",
             color = "black", # Add borders to segments
             #mapping = aes(fill = state_group)
    ) +  
    #theme(legend.position = "none")+
    scale_fill_manual(values = as.character(df_summery_by_state$state_color))

    ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      "plots",
      paste("pie_chart H3K27Ac", "png", sep = ".")
    )
    ,
    plot = p,
    width = 12,
    height = 12
  )
  
  # plot for each sample mean methylation of chromatin state  
  ggplot(data = df_state_summery,mapping = aes(x = sample,y = mean_methylation,color = type))+
    #geom_boxplot()+
    geom_point(mapping = aes(color = Group))

    
    
}

# state_dependen_analysis --------------------------------------------------------
if (state_dependen_analysis) {
  # load data
  all_CpG_data_hg19 <- readRDS(file = file_name <- file.path(
    OUTPUT_FOLDER,
    "results",
    "WholeGenome",
    paste("all_data_complete_with_test", "rds", sep = ".")
  ))
  
  state_names <- unique(as.character(all_CpG_data_hg19$state))
  state_names <- state_names[order(pad_with_zeros(state_names))]
  
  
  # Specify the folder path
  folder_path <- file.path(OUTPUT_FOLDER,
                           "results",
                           "WholeGenome",
                           "analysis_by_state")
  
  # Check if the folder exists
  if (!file.exists(folder_path)) {
    # If the folder doesn't exist, create it
    dir.create(folder_path)
    print(paste("Folder", folder_path, "created successfully."))
  } else {
    print(paste("Folder", folder_path, "already exists."))
  }
  
  summery_all_state_results <- data.frame()
  
  chromatin_state_names <- levels(all_CpG_data_hg19$state)
  
  for (s in 1:length(chromatin_state_names)) {
    start_time <- Sys.time()
    state <- chromatin_state_names[s]
    print(state)
    df_state <-
      all_CpG_data_hg19[all_CpG_data_hg19$state == state, ]
    
    df_state$pearson.p_val[df_state$pearson.delta == 0] <- 1
    print(paste("# CpG in :", state, nrow(df_state)))
    saveRDS(object = df_state, file = file.path(folder_path, paste(state, "rds", sep = ".")))
    
    pearson.state <-
      list(
        sim_typisation = real_age_typisation,
        permuation_order = state,
        data = df_state[, c("pearson.p_val", "pearson.statistic", "pearson.delta")]
      )
    
    saveRDS(object = pearson.state, file = file.path(folder_path,
                                                     paste("pearson", state, "rds", sep = ".")))
    
    file_names <- list.files(path = folder_path,
                             pattern = paste("^pearson", state, "rds", sep = "."))
    
    state_results <-
      parallel_summerize_permutations(sim_file_names = file.path(folder_path, file_names))
    
    
    saveRDS(object = state_results,
            file = file.path(
              folder_path,
              paste("summery", "pearson", state, "rds", sep = ".")
            ))
    
    max_n <-
      state_results$n_signfincant_CpG[state_results$min_delta == 0 &
                                        state_results$minus_log_alpha == 0]
    state_results$fraction_negative_signfincant_CpG <-
      state_results$negative_signfincant_CpG / max_n
    state_results$fraction_positive_signfincant_CpG <-
      state_results$positive_signfincant_CpG / max_n
    state_results$fraction_signfincant_CpG <-
      state_results$n_signfincant_CpG / max_n
    # Convert to long format
    state_results_long <-
      pivot_longer(
        data = state_results,
        cols = c(
          fraction_negative_signfincant_CpG,
          fraction_positive_signfincant_CpG
        ),
        names_to = "slope",
        values_to = "fraction_direction_signfincant"
      )
    
    summery_all_state_results <-
      rbind(summery_all_state_results, state_results_long)
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  saveRDS(object = summery_all_state_results,
          file = file.path(
            folder_path,
            paste("summery", "pearson", "all", "states", "rds", sep = ".")
          ))
  
  
  
  summery_all_state_results_wide <- pivot_wider(data = summery_all_state_results,
                                                names_from = slope,
                                                values_from = fraction_direction_signfincant)
  
  dummy_value <- 1
  summery_all_state_results_wide$negative_enrichment <-
    (summery_all_state_results_wide$negative_signfincant_CpG) / (dummy_value + summery_all_state_results_wide$positive_signfincant_CpG)
  
  alpha <- 8
  min_delta <- 0.3
  ggplot(
    data =  summery_all_state_results_wide[summery_all_state_results_wide$min_delta == min_delta &
                                             summery_all_state_results_wide$minus_log_alpha >= 6,],
    mapping = aes(x = minus_log_alpha,
                  y = negative_enrichment,
                  color = permuation_type)
  ) +
    geom_point(position = "jitter") +
    theme(legend.position = "none")
  
  
  slice_0 <-
    summery_all_state_results_wide[summery_all_state_results_wide$min_delta == min_delta &
                                     summery_all_state_results_wide$minus_log_alpha == 0,]
  
  slice_alpha <-
    summery_all_state_results_wide[summery_all_state_results_wide$min_delta == min_delta &
                                     summery_all_state_results_wide$minus_log_alpha == alpha,]
  
  slice_alpha$percent_negative_significant <-
    slice_alpha$negative_signfincant_CpG / slice_0$negative_signfincant_CpG
  slice_alpha$percent_positive_significant <-
    slice_alpha$positive_signfincant_CpG / slice_0$positive_signfincant_CpG
  
  slice_alpha_long <-
    pivot_longer(
      slice_alpha,
      cols = c(
        percent_negative_significant,
        percent_positive_significant
      ),
      names_to = "slope",
      values_to = "percent_significant"
    )
  slice_alpha_long$permuation_type <-
    factor(x = slice_alpha_long$permuation_type,
           levels = state_names)
  p <- ggplot(
    data = slice_alpha_long,
    mapping = aes(x = permuation_type, y = percent_significant, fill = slope)
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    ggplot2::theme(
      #axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank(),
      #legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    ) +
    geom_hline(
      yintercept = exp(-alpha),
      linetype = "dashed",
      color = "red"
    ) +  # Add horizontal line exp(-alpha)
    annotate(
      geom = "text",
      x = tail(as.factor(
        unique(slice_alpha_long$permuation_type)
      ))[1],
      y = exp(-alpha),
      label = expression(e ^ (-alpha)),
      hjust = -2,
      vjust = 0.5,
      color = "red"
    ) +
    coord_cartesian(clip = "off") +
    ylab(expression(
      paste("fraction of significant CpG ", p[value] <= e ^ (-8), ", ", delta[meth] >= 0.3)
    )) +
    xlab("chromatin state labeling")
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      "plots",
      paste("delta", min_delta, "alpha", alpha, "png", sep = ".")
    ),
    plot = p,
    width = 12,
    height = 6
  )
  
  
  min_delta <- 0.30
  #all_CpG_data_hg19[all_CpG_data_hg19$pearson.delta >=]
  
  
  ggplot(
    data = summery_all_state_results[summery_all_state_results$min_delta == min_delta &
                                       summery_all_state_results$minus_log_alpha >= 6,],
    mapping = aes(
      x = as.character(minus_log_alpha),
      y = fraction_direction_signfincant,
      color = interaction(permuation_type, slope)
    )
  ) +
    geom_point(position = "jitter") +
    theme(legend.position = "none")
  
  slice <-
    summery_all_state_results[summery_all_state_results$min_delta == min_delta &
                                summery_all_state_results$minus_log_alpha == alpha, ]
  slice[order(slice$fraction_direction_signfincant, decreasing = TRUE), c("permuation_type",
                                                                          "slope",
                                                                          "fraction_direction_signfincant")]
  
}
# genome wide only delta good states --------------------------------------------------
if (FALSE) {
  min_delta <- 0.3
  max_p <- 0.001
  genome_CpG_data_hg19 <-
    all_CpG_data_hg19[!grepl(pattern = "GapArtf",
                             x = all_CpG_data_hg19$state,
                             ignore.case = TRUE) &
                        all_CpG_data_hg19$pearson.delta >= min_delta, ]
  
  counts_min_delta <- genome_CpG_data_hg19 %>%
    count(state)
  
  significant_genome_CpG_data_hg19 <-
    genome_CpG_data_hg19[genome_CpG_data_hg19$pearson.p_val <= max_p, ]
  
  counts_sig_min_delta <- significant_genome_CpG_data_hg19 %>%
    count(state)
  
  #plot_age_correlation(df_row = significant_genome_CpG_data_hg19[1,],typisation = real_age_typisation)
  
  significant_genome_CpG_data_hg19$slope <-
    ifelse(significant_genome_CpG_data_hg19$pearson.statistic <= 0,
           "negative",
           "positive")
  
  library(dplyr)
  # Example calculation
  result <- significant_genome_CpG_data_hg19 %>%
    group_by(state) %>%
    summarize(
      fraction_negative = mean(slope == "negative", na.rm = TRUE),
      fraction_positive = mean(slope == "positive", na.rm = TRUE)
    )
  
  counts_sig_min_delta_slope <- significant_genome_CpG_data_hg19 %>%
    group_by(state) %>%
    count(slope)
  
  counts_sig_min_delta_slope$total <-
    rep(x = counts_min_delta$n, each = 2)
  counts_sig_min_delta_slope$fraction <-
    counts_sig_min_delta_slope$n / counts_sig_min_delta_slope$total
  
  matching_indices <-
    match(
      counts_sig_min_delta_slope$state,
      paste(
        df_universal_annotation$state_order_by_group,
        df_universal_annotation$`States  presented in paper`,
        sep = "_"
      )
    )
  counts_sig_min_delta_slope$state_color <-
    factor(df_universal_annotation$color[matching_indices],
           levels = unique(df_universal_annotation$color))
  counts_sig_min_delta_slope$state_group <-
    factor(df_universal_annotation$Group[matching_indices],
           levels = unique(df_universal_annotation$Group))
  
  # sanity check
  # sum(significant_genome_CpG_data_hg19$slope[significant_genome_CpG_data_hg19$state == "99_TSS2"] == "positive",na.rm = TRUE)
  # sum(significant_genome_CpG_data_hg19$slope[significant_genome_CpG_data_hg19$state == "99_TSS2"] == "negative",na.rm = TRUE)
  
  p <- ggplot(
    data = counts_sig_min_delta_slope,
    mapping = aes(
      x = state,
      y = fraction,
      fill = state_group,
      color = slope
    )
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    ggplot2::theme(
      #axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank(),
      #legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    ) +
    # geom_hline(
    #   yintercept = max_p,
    #   linetype = "dashed",
    #   color = "red"
    # ) +  # Add horizontal line exp(-alpha)
    # annotate(
    #   geom = "text",
    #   x = tail(as.factor(unique(slice_alpha_long$permuation_type)))[1],
    #   y = exp(-alpha),
    #   label = expression(e^(-alpha)),
    #   hjust = -2,
  #   vjust = 0.5,
  #   color = "red"
  # )+
  #coord_cartesian(clip = "off")+
  ylab(expression(
    paste("fraction of significant CpG ", p[value] <= 0.001, ", ", delta[meth] >= 0.3)
  )) +
    xlab("chromatin state labeling") +
    scale_fill_manual(values = unique(df_universal_annotation$color))
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      "plots",
      paste("significant CpG whole genome", "png", sep = ".")
    )
    ,
    plot = p,
    width = 12,
    height = 12
  )
  
  summary_df <- counts_sig_min_delta_slope %>%
    group_by(state_group, slope) %>%
    summarize(n = sum(n),
              total_n = sum(total))
  summary_df_wide <- summary_df %>%
    pivot_wider(
      names_from = slope,
      values_from = c(n),
      names_glue = "{slope}_{.value}"
    )
  summary_df_wide$log_fraction <-
    log(summary_df_wide$negative_n / summary_df_wide$positive_n)
  
  p <- ggplot(
    data = summary_df_wide,
    mapping = aes(
      x = state_group,
      y = log_fraction,
      fill = state_group,
      label = paste("Negative:", negative_n, " ; Positive:", positive_n)
    )
  ) +
    geom_bar(stat = "identity", colour = "black") +
    geom_text(
      position = position_stack(vjust = 0.5),
      size = 3,
      mapping = aes(y  =  3)
    ) +  # Add text annotations
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   # axis.text.x = element_text(
                   #   angle = 90,
                   #   vjust = 0.5,
                   #   hjust = 1
                   # )
                   ) +
                   ylab(expression(paste(
                     log(negative[correlation] / positive[correlation]), " ", p[value] <= 0.001, ", ", delta[meth] >= 0.3
                   ))) +
                     xlab("Chromatin state groups") +
                     scale_fill_manual(values = unique(df_universal_annotation$color)) +
                     coord_flip()
                   
                   
                   ggsave(
                     filename = file.path(
                       OUTPUT_FOLDER,
                       "results",
                       "WholeGenome",
                       paste("enrichment significant CpG per Group", "png", sep = ".")
                     )
                     ,
                     plot = p,
                     width = 12,
                     height = 12
                   )
                   
                   
                   summary_df_state <- counts_sig_min_delta_slope %>%
                     group_by(state, state_color, state_group) %>%
                     summarize(
                       negative_count = sum(ifelse(slope == "negative", n, 0)),
                       positive_count = sum(ifelse(slope == "positive", n, 0)),
                       log_fraction = log(negative_count / positive_count)
                     )
                   
                   
                   p <- ggplot(
                     data = summary_df_state,
                     mapping = aes(
                       x = state,
                       y = log_fraction,
                       fill = state_group,
                       label = paste("Negative:", negative_count, " ; Positive:", positive_count)
                     )
                   ) +
                     geom_bar(stat = "identity", colour = "black") +
                     geom_text(
                       position = position_stack(vjust = 0.5),
                       size = 3,
                       mapping = aes(y  =  max(summary_df_state$log_fraction))
                     ) +  # Add text annotations
                     ggplot2::theme(panel.background = ggplot2::element_blank(),
                                    # axis.text.x = element_text(
                                    #   angle = 90,
                                    #   vjust = 0.5,
                                    #   hjust = 1
                                    # )
                                    ) +
                                    ylab(expression(paste(
                                      log(negative[correlation] / positive[correlation]), " ", p[value] <= 0.001, ", ", delta[meth] >= 0.3
                                    ))) +
                                      xlab("Chromatin state groups") +
                                      scale_fill_manual(values = unique(df_universal_annotation$color)) +
                                      coord_flip()
                                    
                                    ggsave(
                                      filename = file.path(
                                        OUTPUT_FOLDER,
                                        "results",
                                        "WholeGenome",
                                        paste("enrichment significant CpG per state", "png", sep = ".")
                                      )
                                      ,
                                      plot = p,
                                      width = 12,
                                      height = 12
                                    )
}
# Calculate true Data -----------------------------------------------------------
if (create_true_stat) {
  
  print("create_true_stat")
  # compute sample numbers per CpG position
  start_time <- Sys.time()
  df_peak_CpG$score <-
    parallel_score_sample_number(df_peak_CpG, real_age_typisation)
  end_time <- Sys.time()
  print(end_time - start_time)
  print(paste("# CpG in H3K27Ac peaks:", nrow(df_peak_CpG)))
  print(paste(
    "# CpG in H3K27Ac peaks with all 39 samples :",
    sum(df_peak_CpG$score == nrow(real_age_typisation))
  ))
  print(paste("discard ", 100 * (
    1 - sum(df_peak_CpG$score == nrow(real_age_typisation)) / nrow(df_peak_CpG)
  ), "% of CpG positions"))
  # discard CpG with not full data
  df_peak_CpG_complete <-
    df_peak_CpG[nrow(real_age_typisation) == df_peak_CpG$score,]
  
  # compute real data value kruskall_valis
  start_time <- Sys.time()
  df_kruskall_valis <-
    parallel_testing_kruskall_valis(df_peak_CpG_complete, real_food_typisation)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # comupute real data value pearson_cor
  start_time <- Sys.time()
  df_testing_pearson_cor <-
    parallel_testing_pearson_cor(df = df_peak_CpG_complete, age.typisation = real_age_typisation)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # save results
  df_peak_CpG_complete_with_test <-
    cbind(df_peak_CpG_complete,
          df_testing_pearson_cor,
          df_kruskall_valis)
  
  
  
  #   df_peak_CpG_complete_with_test <- all_CpG_data_hg19[all_CpG_data_hg19$name != "NO_CHIP",]
  
  saveRDS(
    object = df_peak_CpG_complete_with_test,
    file = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds")
  )
  
  
  saveRDS(
    object = list(
      permuation_order = NA ,
      sim_typisation = real_food_typisation ,
      data = df_peak_CpG_complete_with_test[, KW_test_colnames]
    ),
    file = file.path(OUTPUT_FOLDER, "KW.real.rds")
  )
  saveRDS(
    object = list(
      permuation_order = NA ,
      sim_typisation = real_age_typisation ,
      data = df_peak_CpG_complete_with_test[, pearson_test_colnames]
    ) ,
    file = file.path(OUTPUT_FOLDER, "pearson.real.rds")
  )
  
  for(s in as.character(meta$sample)){
    print(s)
    meth <- as.numeric(df_peak_CpG_complete_with_test[,s])
    sw_normality <- shapiro.test(sample(x = meth,size = 1000,replace = FALSE)) 
    print(sw_normality)
  }
  
  
} else{
  df_peak_CpG_complete_with_test <-
    readRDS(file = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds"))
}
# describe LOCUS ---------------------------------------------------------------
if(FALSE){
  # all_CpG_data_hg19 <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/results/WholeGenome/all_data_complete_with_test.rds")
  
  GH <- list(ID = "GH05J143394",chr = "chr5", start =142774324, end = 142785910)
  Gene <- "NR3C1"
  H3K27ac <- "chr5.142784176.142785514.H3K27ac"
  
  
  df <- all_CpG_data_hg19[all_CpG_data_hg19$chrom == GH$chr &  
                                   GH$start - 1000 <= all_CpG_data_hg19$start &
                            all_CpG_data_hg19$end <= GH$end + 1000,]
  
  unique(df$name)
  
  loc1 <- data.frame(name = "up Risk of being in a poorly regulated neurobehavioural profile", start =  142783501, end =142783640)
  
  loc2 <- data.frame(name = "up Choline", start = 142783501, end = 142783908)
  
  loc3 <- data.frame(name = "up Adult waist circumference, up adult BMI", start = 142782759, end = 142783164)
  
  loc4 <- data.frame(name = "up Meat/fish and vegetable intake, down bread/potato intake in late pregnancy", start = 142783579, end = 142783714)
  
  loc5 <- data.frame(name = "up Adult blood pressure", start = 142783578, end = 142783714)
  
  xmin <- 142782200
  xmax <- 142785600
  
  p <- ggplot(data = df,mapping = aes(x = start,y = -log(pearson.p_val)))+
    geom_rect(aes(xmin = max(GH$start,xmin), xmax = min(xmax,GH$end), ymin = -1, ymax = -0.5), fill = "yellow", alpha = 0.9)+
    geom_segment(aes(x = xmax, y = -0.7, xend = xmax + 100 , yend = -0.7), arrow = arrow(length = unit(0.3, "cm")), color = "yellow")+
    geom_segment(aes(x = xmin, y = -0.7, xend = xmin + 100 , yend = -0.7), arrow = arrow(length = unit(0.3, "cm")), color = "yellow")+
    geom_rect(aes(xmin = 142782230, xmax = 142783990, ymin = -0.5, ymax = -0.1), fill = "LightSkyBlue", alpha = 0.9)+
    
    geom_rect(aes(xmin = 142784176, xmax = 142785514, ymin = -0.5, ymax = -0.1), fill = "LightSkyBlue", alpha = 0.9)+
    geom_bar(stat = "identity")+
    geom_text(aes(x = (xmin + xmax) / 2, y = -0.75, label = GH$ID), color = "black") +
    geom_text(aes(x = (142782230 + 142783990) / 2, y = -0.3, label = "H3K27Ac"), color = "black") +
    geom_text(aes(x = (142784176 + 142785514) / 2, y = -0.3, label = "H3K27Ac"), color = "black")+
    geom_bar(data = df[which.min(df$pearson.p_val),],stat = "identity",color = "red")+
    ylab(expression(-log(pearson_p[value])))+
    xlab("Chr5: Genomic position hg19")+
    xlim(c(xmin,xmax))+
    theme_minimal()+
    # geom_rect(aes(xmin = loc1$start, xmax = loc1$end, ymin = -1.5, ymax = -1.4), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc2$start, xmax = loc2$end, ymin = -1.4, ymax = -1.3), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc3$start, xmax = loc3$end, ymin = -1.3, ymax = -1.2), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc4$start, xmax = loc4$end, ymin = -1.2, ymax = -1.1), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc5$start, xmax = loc5$end, ymin = -1.1, ymax = -1.0), fill = "green", alpha = 0.9)+
    coord_cartesian(expand = c(0, 0))
  
  ggsave(plot = p,filename = file.path(OUTPUT_FOLDER,"results","pearson",paste(Gene,GH$ID,"png",sep = ".")))
  
  # https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-019-0680-7/tables/2
  DMR <- list(ID = "ENO2",chr = "chr12", start =7023752, end = 7024121)
  DMR <- list(ID = "ZNF226",chr = "chr19", start =44669146, end = 44669354)
  DMR <- list(ID = "CCDC51/TMA7",chr = "chr3", start =48481268, end = 48481793)
  
  df <- all_CpG_data_hg19[all_CpG_data_hg19$chrom == DMR$chr &  
                            DMR$start - 1000 <= all_CpG_data_hg19$start &
                            all_CpG_data_hg19$end <= DMR$end + 1000,]
  
  xmin <- DMR$start - 1000
  xmax <- DMR$end + 1000
  
  unique(df$name)
  
  p <- ggplot(data = df,mapping = aes(x = start,y = -log(pearson.p_val)))+
    geom_rect(aes(xmin = max(GH$start,xmin), xmax = min(xmax,GH$end), ymin = -1, ymax = -0.5), fill = "yellow", alpha = 0.9)+
    # geom_segment(aes(x = xmax, y = -0.7, xend = xmax + 100 , yend = -0.7), arrow = arrow(length = unit(0.3, "cm")), color = "yellow")+
    # geom_segment(aes(x = xmin, y = -0.7, xend = xmin + 100 , yend = -0.7), arrow = arrow(length = unit(0.3, "cm")), color = "yellow")+
    geom_rect(aes(xmin = 48481218, xmax = 48482176, ymin = -0.5, ymax = -0.1), fill = "LightSkyBlue", alpha = 0.9)+
    #geom_rect(aes(xmin = 142784176, xmax = 142785514, ymin = -0.5, ymax = -0.1), fill = "LightSkyBlue", alpha = 0.9)+
    geom_bar(stat = "identity")+
    #geom_text(aes(x = (xmin + xmax) / 2, y = -0.75, label = GH$ID), color = "black") +
    #geom_text(aes(x = (142782230 + 142783990) / 2, y = -0.3, label = "H3K27Ac"), color = "black") +
    #geom_text(aes(x = (142784176 + 142785514) / 2, y = -0.3, label = "H3K27Ac"), color = "black")+
    #geom_bar(data = df[which.min(df$pearson.p_val),],stat = "identity",color = "red")+
    ylab(expression(-log(pearson_p[value])))+
    xlab(paste(DMR$chr,": Genomic position hg19"))+
    xlim(c(xmin,xmax))+
    theme_minimal()+
    # geom_rect(aes(xmin = loc1$start, xmax = loc1$end, ymin = -1.5, ymax = -1.4), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc2$start, xmax = loc2$end, ymin = -1.4, ymax = -1.3), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc3$start, xmax = loc3$end, ymin = -1.3, ymax = -1.2), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc4$start, xmax = loc4$end, ymin = -1.2, ymax = -1.1), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc5$start, xmax = loc5$end, ymin = -1.1, ymax = -1.0), fill = "green", alpha = 0.9)+
    coord_cartesian(expand = c(0, 0))
  
}


# describe Data ---------------------------------------------------------------
if (describe_Data) {
  p <- ggplot(data = df_peak_CpG_complete_with_test,
              mapping = aes(x = pearson.delta)) +
    geom_histogram(breaks = seq(from = 0, to = 1, by = 0.01)) +
    theme_minimal() +
    xlab("maximum methylation difference between the samples")
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "maximum methylation difference between the samples.png"
    )
    ,
    plot = p,
    width = 10,
    height = 6
  )
  
  p <- ggplot(data = df_peak_CpG_complete_with_test,
              mapping = aes(x = pearson.delta)) +
    geom_histogram(breaks = seq(from = 0, to = 1, by = 0.01)) +
    theme_minimal() +
    xlab("maximum methylation difference between the samples")
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "maximum methylation difference between the samples.png"
    )
    ,
    plot = p,
    width = 10,
    height = 6
  )
  
  library(data.table)
  long <-
    melt(
      setDT(df_peak_CpG_complete_with_test[as.character(meta$sample)]),
      id.vars = c(),
      variable.name = "sample"
    )
  meth_long <- as.data.frame(long)
  
  p <- ggplot(data = meth_long,
              mapping = aes(x = sample, y = value, color = sample)) +
    geom_violin(trim = FALSE) +
    theme_minimal() +
    xlab("samples") +
    ylab("methylation") +
    coord_flip()
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "plots",
                         "methylation values in samples.png"),
    plot = p,
    width = 8,
    height = 16
  )
  
  
  print(paste("number of CpG in peaks:", nrow(df_peak_CpG)))
  print(paste(
    "number of CpG in peaks with data:",
    nrow(df_peak_CpG_complete_with_test)
  ))
  min_delta = 0.39
  paste(
    "number of CpG in peaks with data with min_delta:",
    sum(df_peak_CpG_complete_with_test$pearson.delta >= min_delta)
  )
  
  
  p <-
    ggplot(data = meta,
           mapping = aes(x = sample, y = age_mean_BP, color = Type)) +
    geom_point() +
    geom_errorbar(aes(ymin = age_mean_BP - age_std_BP, ymax = age_mean_BP +
                        age_std_BP),
                  width = .1) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 1
    )) +
    scale_y_continuous(
      breaks = seq(0, 11000, by = 1000),
      # labels = rep("",13),
      limits = c(0, 11500),
      expand = c(0, 0)
    ) +
    ylab("years before present (B.P.)") +
    xlab("sample ID") +
    scale_color_discrete(name = "lifestyle/diet") +
    theme(
      strip.text.x = element_blank(),
      strip.background = element_rect(colour = "white", fill = "white"),
      legend.position = c(.9, .4)
    ) +
    coord_flip()
  
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "plots",
                         "Age of samples.png"),
    plot = p,
    width = 8,
    height = 8
  )
  #position_nudge(x = -2)
  
  meta$age_max <- meta$age_mean_BP + meta$age_std_BP
  meta$age_min <- meta$age_mean_BP - meta$age_std_BP
  meta$nudge_overlap <- 0
  
  for (r in 1:nrow(meta)) {
    meta$nudge_overlap[r] <-
      sum(meta$age_min[r] < meta$age_min[r:nrow(meta)] &
            meta$age_min[r:nrow(meta)] < meta$age_max[r]) +
      sum(meta$age_min[r] < meta$age_max[r:nrow(meta)] &
            meta$age_max[r:nrow(meta)] < meta$age_max[r])
  }
  
  meta$nudge_overlap <- meta$nudge_overlap / 10
  
  meta$single_factor <- factor("data")
  p <-
    ggplot(
      data = meta,
      mapping = aes(
        x = single_factor,
        y = age_mean_BP,
        color = Type,
        group = single_factor
      )
    ) +
    geom_errorbar(
      position = position_nudge(x = meta$nudge_overlap),
      aes(ymin = age_mean_BP - age_std_BP, ymax = age_mean_BP +
            age_std_BP),
      width = 0.1
    ) +
    geom_point(position = position_nudge(x = meta$nudge_overlap)) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 0,
      vjust = 0.5,
      hjust = 1
    )) +
    scale_y_continuous(
      breaks = seq(0, 11000, by = 1000),
      # labels = rep("",13),
      limits = c(0, 11500),
      expand = c(0, 0)
    ) +
    ylab("years before present (B.P.)") +
    xlab("sample ID") +
    scale_color_discrete(name = "lifestyle/diet") +
    theme(
      strip.text.x = element_blank(),
      strip.background = element_rect(colour = "white", fill = "white"),
      legend.position = c(.9, .4)
    ) + coord_flip()
  
  
  # ggsave(
  #   filename = file.path(OUTPUT_FOLDER,
  #                        "plots",
  #                        "Age of samples.png"),
  #   plot = p,
  #   width = 8,
  #   height = 8
  # )
  
  x_max <- 12000
  
  ggplot(data.frame(meta), aes(x = age_mean_BP, y = 0, color = Type)) +
    #geom_point(size = 5,shape = 25, mapping = aes(y = 0.05))  +
    annotate(
      "segment",
      x = 1,
      xend = x_max,
      y = 0,
      yend = 0,
      size = 2
    ) +
    annotate(
      "segment",
      x = 1,
      xend = 1,
      y = -0.1,
      yend = 0.1,
      size = 2
    ) +
    annotate(
      "segment",
      x = x_max,
      xend = x_max,
      y = -0.1,
      yend = 0.1,
      size = 2
    ) +
    scale_shape_identity() +
    geom_point(shape = 108, size = 4) +
    ggrepel::geom_label_repel(aes(label = sample), col = "black") +
    scale_x_continuous(breaks = seq(0, x_max, by = 1000),
                       limits = c(0, x_max),
    ) +
    scale_y_continuous(limits = c(-3, 3)) +
    theme(
      panel.background = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) + xlab("years before present (B.P.)") +
    scale_color_discrete(name = "lifestyle/diet")
  
  
  p <-
    ggplot(data.frame(meta), aes(x = age_mean_BP, y = 0, fill = Type)) +
    geom_point(size = 5,
               shape = 25,
               mapping = aes(y = 0.05))  +
    annotate(
      "segment",
      x = 1,
      xend = x_max,
      y = 0,
      yend = 0,
      size = 2
    ) +
    annotate(
      "segment",
      x = 1,
      xend = 1,
      y = -0.1,
      yend = 0.1,
      size = 2
    ) +
    annotate(
      "segment",
      x = x_max,
      xend = x_max,
      y = -0.1,
      yend = 0.1,
      size = 2
    ) +
    #scale_shape_identity() +
    #geom_point(shape = 108, size = 4)+
    #ggrepel::geom_label_repel(aes(label = sample), col = "black") +
    scale_x_continuous(breaks = seq(0, x_max, by = 1000),
                       limits = c(0, x_max),
    ) +
    scale_y_continuous(limits = c(-3, 3)) +
    theme(
      panel.background = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) + xlab("years before present (B.P.)") +
    scale_fill_discrete(name = "lifestyle/diet")
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "plots",
                         "Age of samples 1D.png"),
    plot = p,
    width = 12,
    height = 4
  )
  
}
# Calculate column permutation simulations --------------------------------------------------

if (create_permutations_horizonal_columns) {
  print("create_permutations_horizonal_columns")
  dir.create(path = file.path(OUTPUT_FOLDER, "simulation"),
             showWarnings = FALSE)
  n_repetitions <- 9
  for (rep in 1:n_repetitions) {
    # rep <- 1
    start_time <- Sys.time()
    print(rep)
    
    # permutation_age  <-
    #   create_horizontal_permutated_typisation(real_age_typisation)
    permutation_food  <-
      create_horizontal_permutated_typisation(real_food_typisation)
    
    #sim_age_typisation <- permutation_age$sim_typisation
    sim_food_typisation <- permutation_food$sim_typisation
    
    # permutation_age$data <-
    #   parallel_testing_pearson_cor(df = df_peak_CpG_complete_with_test, age.typisation = sim_age_typisation)
    
    permutation_food$data <-
      parallel_testing_kruskall_valis(df = df_peak_CpG_complete_with_test, food.typisation = sim_food_typisation)
    
    # # save age permutation
    # sim_age_file_name <-
    #   paste(
    #     "pearson",
    #     "horizontal",
    #     paste(permutation_age$permuation_order, collapse = "_"),
    #     "rds",
    #     sep = "."
    #   )
    # saveRDS(
    #   object = permutation_age,
    #   file = file.path(OUTPUT_FOLDER, "simulation", sim_age_file_name)
    # )
    
    # save food permutation
    sim_food_file_name <-
      paste(
        "KW",
        "horizontal",
        paste(permutation_food$permuation_order, collapse = "_"),
        "rds",
        sep = "."
      )
    saveRDS(
      object = permutation_food,
      file = file.path(OUTPUT_FOLDER, "simulation", sim_food_file_name)
    )
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
}

# Calculate CpG permuation simulations --------------------------------------------------

if (create_CpG_permutations_vertical) {
  print("create_CpG_permutations_vertical")
  dir.create(path = file.path(OUTPUT_FOLDER, "simulation"),
             showWarnings = FALSE)
  n_repetitions <- 197
  for (rep in 1:n_repetitions) {
    # rep <- 1
    start_time <- Sys.time()
    print(rep)
    
    time_string <- format(start_time, "%Y_%m_%d_%H_%M_%S")
    
    permutation_CpG <-
      create_CpG_permution(df_peak_CpG_complete_with_test, real_age_typisation)
    
    df_permuted <- permutation_CpG$data
    
    permutation_age <- list()
    permutation_age$typisation <- real_age_typisation
    permutation_age$permuation_order <-
      permutation_CpG$permuation_order
    
    permutation_food <- list()
    permutation_food$typisation <- real_food_typisation
    permutation_food$permuation_order <-
      permutation_CpG$permuation_order
    
    permutation_age$data <-
      parallel_testing_pearson_cor(df = df_permuted, age.typisation = real_age_typisation)
    permutation_food$data <-
      parallel_testing_kruskall_valis(df = df_permuted, food.typisation = real_food_typisation)
    
    # save age permutation
    sim_age_file_name <-
      paste("pearson",
            "CpG_permutation",
            time_string,
            "rds",
            sep = ".")
    saveRDS(
      object = permutation_age,
      file = file.path(OUTPUT_FOLDER, "simulation", sim_age_file_name)
    )
    
    # save food permutation
    sim_food_file_name <-
      paste("KW",
            "CpG_permutation",
            time_string,
            "rds",
            sep = ".")
    saveRDS(
      object = permutation_food,
      file = file.path(OUTPUT_FOLDER, "simulation", sim_food_file_name)
    )
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
}


# Evaluate CpG simulations ----------------------------------------------------
if (summerize_CpG_columns_permutations) {
  print("summerize_CpG_columns_permutations")
  # list.files(path =  file.path(OUTPUT_FOLDER, "simulation"))
  
  file_names <-
    list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
               pattern = "*.CpG_permutation.*")
  start_time <- Sys.time()
  sim_results <-
    parallel_summerize_permutations(sim_file_names = file.path(OUTPUT_FOLDER, "simulation", file_names))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  saveRDS(
    object = sim_results,
    file = file.path(OUTPUT_FOLDER, "CpG_permutation.sim_results.rds")
  )
  
  start_time <- Sys.time()
  CpGs_with_delta_cutoff <-
    parallel_count_CpGs_with_delta_cutoff(
      sim_file_names = file.path(OUTPUT_FOLDER, "simulation", file_names),
      min_deltas = seq(0, 0.50, 0.01)
    )
  end_time <- Sys.time()
  print(end_time - start_time)
  
  saveRDS(
    object = CpGs_with_delta_cutoff,
    file = file.path(
      OUTPUT_FOLDER,
      "CpG_permutation.CpGs_with_delta_cutoff.rds"
    )
  )
  
}

# histogramm_plots ##############
if (histogramm_plots) {
  delta <- 0.39
  cut_off <- 6
  alpha_p <- 9.3
  
  #real_results <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/real_results.rds")
  sim_results <-
    readRDS(file = file.path(OUTPUT_FOLDER, "CpG_permutation.sim_results2023.rds"))
  sim_results_pearson_horizontal <-
    readRDS(file = file.path(OUTPUT_FOLDER, "horizontal.sim_results2023.rds"))
  
  sim_results_pearson_CpG <-
    sim_results[sim_results$test == "pearson",]
  sim_results_pearson_horizontal <-
    sim_results_pearson_horizontal[sim_results_pearson_horizontal$test == "pearson",]
  
  
  df_plot <-
    rbind(sim_results_pearson_CpG[, c("min_delta" , "n_signfincant_CpG", "name")],
          sim_results_pearson_CpG[, c("min_delta" , "n_signfincant_CpG", "name")])
  df_plot_horizontal <-
    rbind(sim_results_pearson_horizontal[, c("min_delta" , "n_signfincant_CpG", "name")],
          sim_results_pearson_horizontal[, c("min_delta" , "n_signfincant_CpG", "name")])
  
  df_plot$significant <-
    c(
      sim_results_pearson_CpG$positive_signfincant_CpG,
      sim_results_pearson_CpG$negative_signfincant_CpG
    )
  df_plot_horizontal$significant <-
    c(
      sim_results_pearson_horizontal$positive_signfincant_CpG,
      sim_results_pearson_horizontal$negative_signfincant_CpG
    )
  
  df_plot$alpha <- c(
    sim_results_pearson_CpG$minus_log_alpha,
    -sim_results_pearson_CpG$minus_log_alpha
  )
  df_plot_horizontal$alpha <-
    c(
      sim_results_pearson_horizontal$minus_log_alpha,
      -sim_results_pearson_horizontal$minus_log_alpha
    )
  
  
  
  # prepare real resluts
  real_results_pearson <-
    real_results[real_results$test == "pearson",]
  
  
  df_plot_real <-
    rbind(real_results_pearson[, c("min_delta" , "n_signfincant_CpG", "name")],
          real_results_pearson[, c("min_delta" , "n_signfincant_CpG", "name")])
  
  df_plot_real$significant <-
    c(
      real_results_pearson$positive_signfincant_CpG,
      real_results_pearson$negative_signfincant_CpG
    )
  
  df_plot_real$alpha <-
    c(real_results_pearson$minus_log_alpha,
      -real_results_pearson$minus_log_alpha)
  
  df_plot_real_1 <- df_plot_real[df_plot_real$min_delta == delta,]
  df_plot_1 <- df_plot[df_plot$min_delta == delta,]
  df_plot_horizontal_1 <-
    df_plot_horizontal[df_plot_horizontal$min_delta == delta,]
  
  
  p <-
    ggplot(data = df_plot_1[abs(as.numeric(df_plot$alpha)) > cut_off,] ,
           mapping = aes(x = alpha, y = significant, group = alpha)) +
    geom_boxplot(color = "green", outlier.shape = NA) +
    geom_boxplot(
      data = df_plot_horizontal_1[abs(as.numeric(df_plot_horizontal_1$alpha)) > cut_off,],
      mapping = aes(x = alpha, y = significant, group = alpha),
      color = "blue",
      outlier.shape = NA
    ) +
    geom_point(data = df_plot_real_1[abs(as.numeric(df_plot_real_1$alpha)) > cut_off,], color = "red") +
    theme_minimal() +
    scale_color_discrete(#name = "Dose",
      labels = c("CpG permuatation", "age permuation", "real data"))
  
  folder_name <- paste("delta", delta, "minLogP", alpha_p, sep = "_")
  file_name <- "CDF_p_value destribution.png"
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER, "results", "pearson", folder_name, file_name),
    plot = p,
    width = 8,
    height = 8
  )
  
  
  # df_plot_single_delta<- df_plot[abs(as.numeric(df_plot$alpha)) == alpha_p,]
  # df_plot_real <- df_plot_real[abs(as.numeric(df_plot_real$alpha)) == alpha_p,]
  # df_plot_horizontal <- df_plot_horizontal[abs(as.numeric(df_plot_horizontal$alpha)) == alpha_p,]
  
  
  # df_plot_CpG<- df_plot[abs(as.numeric(df_plot$alpha)) == alpha_p,]
  # df_plot_real <- df_plot_real[abs(as.numeric(df_plot_real$alpha)) == alpha_p,]
  # df_plot_horizontal <- df_plot_horizontal[abs(as.numeric(df_plot_horizontal$alpha)) == alpha_p,]
  
  df_plot_CpG <- df_plot
  df_plot_CpG$slope <-
    ifelse(df_plot_CpG$alpha < 0, yes =  "negative", no = "positive")
  df_plot_CpG$name <- "pearson.CpG_permutation"
  #df_plot_CpG$full_name <- df_plot$name
  
  df_plot_real$slope <-
    ifelse(df_plot_real$alpha < 0, yes =  "negative", no = "positive")
  
  df_plot_horizontal$slope <-
    ifelse(df_plot_horizontal$alpha < 0, yes =  "negative", no = "positive")
  df_plot_horizontal$name <- "pearson.horizontal"
  
  df_plot_all <-
    rbind(df_plot_CpG , df_plot_real, df_plot_horizontal)
  
  df_plot_all$minus_log_p <- abs(df_plot_all$alpha)
  df_plot_all$alpha_string <- as.character(df_plot_all$minus_log_p)
  
  
  #df_plot_single_delta$alpha <- factor(df_plot_single_delta$alpha,levels = sort(unique(df_plot_single_delta$alpha)))
  #df_plot_single_delta$alpha <- as.character(df_plot_single_delta$alpha)
  
  ggplot(data = df_plot_all[df_plot_all$min_delta == delta &
                              as.numeric(df_plot_all$alpha) > cut_off,],
         mapping = aes(x = alpha, y = significant, fill = name)) +
    geom_boxplot(outlier.shape = NA)
  
  
  df_plot_single <- df_plot_all[df_plot_all$min_delta == delta &
                                  df_plot_all$minus_log_p == alpha_p,
                                c("name", "slope" , "significant")]
  
  p <- ggplot(data = df_plot_single,
              mapping = aes(x = slope, y = significant, colour = name)) +
    geom_boxplot(outlier.shape = NA) + #
    #geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
    theme_minimal() +
    ylab(expression(paste(
      "# of significant CpGs:  ", p[value] <= e ^ (-9.3) , ", ", delta[meth] >= 0.39
    ))) +
    xlab("pearson correlation slope") +
    scale_color_manual(
      values  = RColorBrewer::brewer.pal(3, "Set1"),
      name = "CpGs deemed as singificant",
      labels = c(
        "CpG position permuation",
        "sampel age permuation",
        "unpermuted  data"
      )
    )
  
  
  file_name <- "count destribution age and CpG permutation.png"
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER, "results", "pearson", folder_name, file_name),
    plot = p,
    width = 8,
    height = 8
  )
  
  df_plot_simmulations <-
    df_plot_all[df_plot_all$min_delta == delta &
                  df_plot_all$minus_log_p >=  cut_off &
                  df_plot_all$minus_log_p < 10,
                c("name",
                  "slope" ,
                  "significant",
                  "minus_log_p",
                  "alpha_string")]
  
  #interaction_colors <- c( "brown1", "deepskyblue", "chartreuse", "brown3", "deepskyblue3", "chartreuse3")
  interaction_colors <-
    c("brown1", "deepskyblue",  "brown3", "deepskyblue3")
  interaction_colors1 <- c("brown1",  "brown3")
  interaction_colors2 <- c("deepskyblue", "deepskyblue3")
  
  p <-
    ggplot(
      data = df_plot_simmulations[df_plot_simmulations$name != "pearson.CpG_permutation",],
      mapping = aes(
        x = alpha_string,
        y = significant,
        color = interaction(name, slope)
      )
    ) +
    geom_boxplot(outlier.shape = NA) +
    xlab(label = expression("significance cut off":alpha)) +
    ylab(label = "number of CpG passing the significance threshold") +
    scale_color_manual(
      values = interaction_colors,
      name = expression("Number of p-values" < e ^ -alpha),
      #"number of p values < exp(-alpha)",
      labels = c(
        #'CpG permutation with negative pearson correlation',
        'Age permutation with negative pearson correlation',
        'data with negative pearson correlation',
        #'CpG permutation with positive pearson correlation',
        'Age permutation with positive pearson correlation',
        'data with positive pearson correlation'
      )
    ) +
    theme_minimal() +
    theme(legend.position = c(.8, .5))
  
  file_name <- "all alpha count destribution age permutation.png"
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER, "results", "pearson", folder_name, file_name),
    plot = p,
    width = 12,
    height = 8
  )
  
  real_file_names <- list.files(path = file.path(OUTPUT_FOLDER),
                                pattern = "*.real.*")
  start_time <- Sys.time()
  real_CpGs_with_delta_cutoff <-
    parallel_count_CpGs_with_delta_cutoff(
      sim_file_names = file.path(OUTPUT_FOLDER, real_file_names),
      min_deltas = seq(0, 0.50, 0.01)
    )
  end_time <- Sys.time()
  print(end_time - start_time)
  
  saveRDS(
    object =  real_CpGs_with_delta_cutoff,
    file = file.path(OUTPUT_FOLDER, "real.CpGs_with_delta_cutoff.rds")
  )
  
  CpGs_with_delta <-
    rbind(real_CpGs_with_delta_cutoff[real_CpGs_with_delta_cutoff$test == "pearson", ],
          CpGs_with_delta_cutoff[CpGs_with_delta_cutoff$test == "pearson", ])
  
  # ggplot(
  #   data = CpGs_with_delta,
  #   mapping = aes(
  #     x = min_delta,
  #     y = n_CpG,
  #     group = min_delta,
  #     fill = permuation_type,
  #     color = permuation_type
  #   )
  # ) +
  #   geom_boxplot() +
  #   geom_point() +
  #   #ylab(bquote("Number of CpG with methylation variance " ~ delta[meth] >=)) +
  #   xlab(bquote("methylation variance cut off value: " ~ delta[meth]))
  
  
  # ggplot(
  #   data = CpGs_with_delta,
  #   mapping = aes(
  #     x = min_delta,
  #     y = n_CpG,
  #     group = min_delta,
  #     fill = permuation_type,
  #     color = permuation_type
  #   )
  # ) +
  #   geom_boxplot() +
  #   geom_point() +
  #   theme_minimal()
  
  legend_title <-
    bquote("CDF of CpGs \n with minimum methylation variance" ~ delta[meth])
  
  p <-  ggplot(
    data = CpGs_with_delta,
    mapping = aes(
      x = min_delta,
      y = n_CpG,
      # /max(CpGs_with_delta$n_CpG)
      group = min_delta,
      fill = permuation_type,
      color = permuation_type
    )
  ) +
    geom_boxplot() +
    geom_point() +
    theme_minimal() +
    ylab(bquote("Number of CpG with methylation variance of at least "  ~ delta[meth])) +
    #  ylab(expression(paste("Number of CpG with methylation variance ", delta[meth] >=)))+
    xlab(bquote("minimum required methylation variability value: " ~ delta[meth])) +
    scale_fill_discrete(name = legend_title,
                        labels = c("CpG position permutation", "orignal data")) +
    scale_color_discrete(name = legend_title,
                         labels = c("CpG position permutation", "orignal data")) +
    theme(legend.position = c(.8, .5))
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "pearson",
      folder_name,
      "CFD methylation variablity.png"
    ),
    plot = p,
    width = 12,
    height = 8
  )
  
  
  file_names <-
    list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
               pattern = "*.CpG_permutation.*")
  
  
  start_time <- Sys.time()
  count_CpGs_with_delta_histogramm  <-
    parallel_count_CpGs_with_delta_histogramm(sim_file_names = file.path(OUTPUT_FOLDER, "simulation", file_names))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  
  p <- ggplot(data = df_peak_CpG_complete_with_test,
              mapping = aes(x = pearson.delta)) +
    geom_histogram(breaks = seq(from = 0, to = 1, by = 0.01)) +
    theme_minimal() +
    xlab("maximum methylation difference between the samples")
  
  df <- get_hist(p)
  df$name <- "real"
  df$y_fraction <- df$y / sum(df$y)
  
  ggplot(data = df,
         mapping = aes(x = x, y = y_fraction)) +
    geom_point() +
    theme_minimal() +
    xlab(bquote("methylation variance: "  ~ delta[meth]))
  
  
  
  legend_title <-
    bquote("destribution of methylation variability" ~ delta[meth])
  p <-
    ggplot(data = count_CpGs_with_delta_histogramm[count_CpGs_with_delta_histogramm$test == "pearson",],
           mapping = aes(
             x = x,
             y = y,
             group = x,
             color = test
           )) +
    geom_boxplot() +
    geom_point(data = df,
               mapping = aes(
                 x = x,
                 y = y,
                 group = x,
                 color = name
               )) +
    theme_minimal() +
    xlab(bquote("methylation variance: "  ~ delta[meth])) +
    ylab("number of CpG positions") +
    scale_color_discrete(name = legend_title,
                         labels = c("CpG position permutation", "orignal data")) +
    theme(legend.position = c(.8, .5))
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "pearson",
      folder_name,
      "destribution methylation variablity.png"
    ),
    plot = p,
    width = 12,
    height = 8
  )
  
  
  # ggsave(
  #   filename = file.path(
  #     OUTPUT_FOLDER,
  #     "plots",
  #     "maximum methylation difference between the samples.png"
  #   )
  #   ,
  #   plot = p,
  #   width = 10,
  #   height = 6
  # )
  
  
  
  ggplot(
    data = CpGs_with_delta,
    mapping = aes(
      x = min_delta,
      y = n_CpG,
      group = min_delta ,
      color = permuation_type
    )
  ) +
    geom_point()
  
  
  mean_number_of_CpG_permutation_with_delta <-
    mean(CpGs_with_delta[CpGs_with_delta$min_delta == min_delta &
                           CpGs_with_delta$permuation_type == "CpG_permutation", "n_CpG"])
  number_of_CpG_with_delta <-
    mean(CpGs_with_delta[CpGs_with_delta$min_delta == min_delta &
                           CpGs_with_delta$permuation_type == "real", "n_CpG"])
  
  df_plot_CpG_permutations <-
    df_plot_all[("pearson.CpG_permutation" == df_plot_all$name |
                   "pearson.real" == df_plot_all$name) &
                  df_plot_all$min_delta == min_delta &
                  df_plot_all$minus_log_p >=  cut_off &
                  df_plot_all$minus_log_p < 10,
                c("min_delta",
                  "name",
                  "slope" ,
                  "significant",
                  "minus_log_p",
                  "alpha_string")]
  
  
  df_plot_CpG_permutations$percent_significant <- NA
  df_plot_CpG_permutations$percent_significant[df_plot_CpG_permutations$name == "pearson.CpG_permutation"] <-
    df_plot_CpG_permutations$significant[df_plot_CpG_permutations$name == "pearson.CpG_permutation"] /
    mean_number_of_CpG_permutation_with_delta * 100
  df_plot_CpG_permutations$percent_significant[df_plot_CpG_permutations$name == "pearson.real"] <-
    df_plot_CpG_permutations$significant[df_plot_CpG_permutations$name == "pearson.real"] /
    number_of_CpG_with_delta * 100
  
  p <- ggplot(
    data = df_plot_CpG_permutations,
    mapping = aes(
      x = alpha_string,
      y = percent_significant,
      color = interaction(name, slope)
    )
  ) +
    geom_boxplot(outlier.shape = NA) +
    xlab(label = expression("significance cut off":alpha)) +
    ylab(expression(
      paste(
        "% of CpG with at least",
        delta[meth],
        " = 39% passing significance threshold"
      )
    )) +
    scale_color_manual(
      values = interaction_colors,
      name = expression("percent of p-values" < e ^ -alpha),
      #"number of p values < exp(-alpha)",
      labels = c(
        'CpG permutation with negative pearson correlation',
        #'Age permutation with negative pearson correlation',
        'data with negative pearson correlation',
        'CpG permutation with positive pearson correlation',
        #'Age permutation with positive pearson correlation',
        'data with positive pearson correlation'
      )
    ) +
    theme_minimal() +
    theme(legend.position = c(.8, .5))
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "pearson",
      folder_name,
      "CpG permutation vs data.png"
    ),
    plot = p,
    width = 12,
    height = 8
  )
}
# Evaluate horizontal simulations ----------------------------------------------------
if (summerize_horizonal_columns_permutations) {
  test <- "pearson"
  
  print("summerize_horizonal_columns_permutations")
  # list.files(path =  file.path(OUTPUT_FOLDER, "simulation"))
  
  # file_names <-
  #   list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
  #              pattern = "*.horizontal.*")
  
  file_names <-
    list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
               pattern = paste(test,"horizontal.*",sep = "."))
  
  start_time <- Sys.time()
  sim_results <-
    parallel_summerize_permutations(sim_file_names = file.path(OUTPUT_FOLDER, "simulation", file_names))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  # saveRDS(
  #   object = sim_results,
  #   file = file.path(OUTPUT_FOLDER, "horizontal.sim_results.rds")
  # )
  
  saveRDS(
    object = sim_results,
    file = file.path(OUTPUT_FOLDER, paste(test,"horizontal.sim_results.rds",sep = "."))
  )
}

# Evaluate real values ----------------------------------------------------
if (summerize_real_values) {
  print("summerize_real_values")
  file_names <- list.files(path = file.path(OUTPUT_FOLDER),
                           pattern = "*.real.*")
  start_time <- Sys.time()
  real_results <-
    parallel_summerize_permutations(sim_file_names = file.path(OUTPUT_FOLDER, file_names))
  end_time <- Sys.time()
  print(end_time - start_time)
  saveRDS(object = real_results,
          file = file.path(OUTPUT_FOLDER, "real_results.rds"))
  
  
  # df_plot <- df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$pearson.p_val < 0.01,]
  # ggplot(data = df_plot ,mapping = aes(x = pearson.statistic))+
  #          geom_histogram(binwidth = 0.1)
  
  
  
  #pearson.real <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/real_results.rds")
  
}

# aggregate CpG results -------------------------------------------------------
if (aggregate_CpG_results) {
  real_results <-
    readRDS(file = file.path(OUTPUT_FOLDER, "real_results.rds"))
  #sim_results <-
  #  readRDS(file = file.path(OUTPUT_FOLDER, "horizontal.sim_results.rds"))
  sim_results <-
    readRDS(file = file.path(OUTPUT_FOLDER, "CpG_permutation.sim_results.rds"))
  
  
  df_empirical_means <-
    aggregate(
      n_signfincant_CpG ~ min_delta + minus_log_alpha + permuation_type + test,
      sim_results,
      FUN =
        mean
    )
  
  df_empirical_means <-
    df_empirical_means[order(
      df_empirical_means$min_delta,
      df_empirical_means$minus_log_alpha,
      df_empirical_means$test
    ),]
  
  df_empirical_means$real_count <-
    real_results$n_signfincant_CpG[order(real_results$min_delta,
                                         real_results$minus_log_alpha,
                                         real_results$test)]
  
  df_empirical_means$FDR <-
    df_empirical_means$n_signfincant_CpG / df_empirical_means$real_count
  
  saveRDS(
    object = df_empirical_means,
    file = file.path(OUTPUT_FOLDER, "df_CpG_vertical_empirical_means.rds")
  )
  #saveRDS(object = df_empirical_means,file = file.path(OUTPUT_FOLDER,"df_horizontal_column_empirical_means.rds"))
  
  ggplot(
    data = sim_results[sim_results$min_delta == 0.39,],
    mapping = aes(x = minus_log_alpha,
                  y = n_signfincant_CpG,
                  color = name)
  ) +
    geom_point() +
    theme(legend.position = "none")
  
  
  for (delta in unique(sim_results$min_delta)) {
    print(delta)
    sims_delta <- sim_results[sim_results$min_delta == delta, ]
    
    for (minus_log_alpha in seq(from = 5, to = 10, by = 1)) {
      print(head(real_results[real_results$min_delta == delta &
                                real_results$minus_log_alpha == minus_log_alpha, ]))
      for (t in unique(sims_delta$test)) {
        print(t)
        sims <-
          sims_delta[sims_delta$minus_log_alpha == minus_log_alpha &
                       sims_delta$test == t, ]
        sims <-
          sims[order(sims$n_signfincant_CpG, decreasing = TRUE), ]
        print(head(sims, 3))
        print(nrow(sims))
      }
    }
  }
  
  
}
# aggregate horizontal results -------------------------------------------------------
if (aggregate_horizontal_results) {
  test <- "pearson"
  real_results <-
    readRDS(file = file.path(OUTPUT_FOLDER, "real_results.rds"))
  # sim_results <-
  #   readRDS(file = file.path(OUTPUT_FOLDER, "horizontal.sim_results.rds"))
  
  #sim_results <- readRDS(file = file.path(OUTPUT_FOLDER,"CpG_permutation.sim_results.rds"))
  
  
  df_empirical_means <-
    aggregate(
      n_signfincant_CpG ~ min_delta + minus_log_alpha + permuation_type + test,
      sim_results,
      FUN =
        mean
    )
  
  df_empirical_means <-
    df_empirical_means[order(
      df_empirical_means$min_delta,
      df_empirical_means$minus_log_alpha,
      df_empirical_means$test
    ), ]
  
  
  df_empirical_means <-
    df_empirical_means[df_empirical_means$test == test, ]
  
  real_results_test <- real_results[test ==  real_results$test, ]
  #ggplot(data = real_results_test[real_results_test$min_delta == 0.3,],mapping = aes(x = minus_log_alpha,y = n_signfincant_CpG))+geom_point()
  
  
  df_empirical_means$real_count <-
    real_results_test$n_signfincant_CpG[order(real_results_test$min_delta,
                                              real_results_test$minus_log_alpha)]
  

  # df_empirical_means$real_count <-
  #   real_results$n_signfincant_CpG[order(real_results$min_delta,
  #                                        real_results$minus_log_alpha,
  #                                        real_results$test)]

  df_empirical_means$FDR <-
    df_empirical_means$n_signfincant_CpG / df_empirical_means$real_count
  
  #saveRDS(object = df_empirical_means,file = file.path(OUTPUT_FOLDER,"df_CpG_vertical_empirical_means.rds"))
  saveRDS(
    object = df_empirical_means,
    file = file.path(OUTPUT_FOLDER, paste(test,"df_horizontal_column_empirical_means.rds",sep = "."))
  )
}
# plot landscape -------------------------------------------------------
if (plot_landscape) {
  min_alpha <- 6.5
  max_alpha <- 9.9
  min_delta <- 0.2
  test <- "KW"
  min_real_count <- 5
  
  
  # min(best_CpGs$start)
  # [1] 29633535
  # > max(best_CpGs$start)
  # [1] 29633767
  
  search_df <-
    df_empirical_means[df_empirical_means$minus_log_alpha >= min_alpha &
                         df_empirical_means$minus_log_alpha <= max_alpha &
                         df_empirical_means$min_delta > min_delta &
                         df_empirical_means$test == test 
                      & min_real_count <= df_empirical_means$real_count 
                      # &  df_empirical_means$n_signfincant_CpG > 0
                      , ]
  
  plotly::plot_ly(
    search_df ,
    x = ~ min_delta ,
    y = ~ minus_log_alpha,
    z = ~ FDR ,
    marker = list(
      color = ~ FDR,
      colorscale = c('#683531', '#FFE1A1'),
      showscale = TRUE
    )
  )
  
  best_values <- search_df[which.min(search_df$FDR), ]
  best_CpGs <-
    df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$kw.delta >= best_values$min_delta &
                                     df_peak_CpG_complete_with_test$kw.p_val <= exp(-best_values$minus_log_alpha),]
  
  folder_name <- paste("delta",
                       best_values$min_delta,
                       "minLogP",
                       best_values$minus_log_alpha,
                       sep = "_")
  
  dir.create(
    path = file.path(OUTPUT_FOLDER, "results", test, folder_name),
    showWarnings = FALSE
  )
  
  
  saveRDS(
    object = best_CpGs,
    file = file.path(OUTPUT_FOLDER,
                     "results",
                     test,
                     folder_name,
                     "best_CpGs.rds")
  )
  for (r in 1:nrow(best_CpGs)) {
    p <-
      plot_food_correlation(df_row = best_CpGs[r, ], typisation = real_food_typisation)
    ggsave(
      plot = p,
      filename = file.path(
        OUTPUT_FOLDER,
        "results",
        test,
        folder_name,
        paste(best_CpGs$chrom[r],
              best_CpGs$start[r],
              ".png",
              sep = ".")
      ),
      width = 5,
      height = 5
    )
    
    p <-
      plot_age_correlation(df_row = best_CpGs[r, ], typisation = real_age_typisation)
    ggsave(
      plot = p,
      filename = file.path(
        OUTPUT_FOLDER,
        "results",
        test,
        folder_name,
        paste("linear",
              best_CpGs$chrom[r],
              best_CpGs$start[r],
              ".png",
              sep = ".")
      )
    )
    
  }
}
# plot pearson landscape --------------------------------------------------
if (pearson_landscape) {
  min_alpha <- 6.5
  max_alpha <- 9.9
  min_delta <- 0.2
  test <- "pearson"
  min_real_count <- 1
  
  search_df <-
    df_empirical_means[df_empirical_means$minus_log_alpha >= min_alpha &
                         df_empirical_means$minus_log_alpha <= max_alpha &
                         df_empirical_means$min_delta > min_delta &
                         df_empirical_means$test == test &
                         min_real_count <= df_empirical_means$real_count, ]
  
  # debug
  #
  # debug_df <-
  #   df_empirical_means[  df_empirical_means$min_delta == 0.3 &
  #                        df_empirical_means$test == test,]
  #
  #ggplot(data = debug_df,mapping = aes(x = minus_log_alpha,y = real_count))+geom_point()
  
  plotly::plot_ly(
    search_df ,
    x = ~ min_delta ,
    y = ~ minus_log_alpha,
    z = ~ FDR ,
    marker = list(
      color = ~ FDR,
      colorscale = c('#683531', '#FFE1A1'),
      showscale = TRUE
    )
  )
  
  
  # index <- all_CpG_data_hg19$name != "NO_CHIP" &
  #   all_CpG_data_hg19$pearson.delta >= 0.39 &
  #   all_CpG_data_hg19$pearson.p_val <= exp(-9.3)
  # 
  # best_CpGs <-
  #   all_CpG_data_hg19[index,]
  
  # best_CpGs2 <-
  #   df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$pearson.delta >= 0.3 &
  #                                    df_peak_CpG_complete_with_test$pearson.p_val <= exp(-9.3), ]
  
  
  best_values <- search_df[which.min(search_df$FDR),]
  best_CpGs <-
    df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$pearson.delta >= best_values$min_delta &
                                     df_peak_CpG_complete_with_test$pearson.p_val <= exp(-best_values$minus_log_alpha), ]
  
  

  folder_name <- paste("delta",
                       best_values$min_delta,
                       "minLogP",
                       best_values$minus_log_alpha,
                       sep = "_")
  
  dir.create(
    path = file.path(OUTPUT_FOLDER, "results", "pearson", folder_name),
    showWarnings = FALSE
  )
  
  
  saveRDS(
    object = best_CpGs,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      "pearson",
      folder_name,
      "best_CpGs.rds"
    )
  )
  for (r in 1:nrow(best_CpGs)) {
    p <-
      plot_age_correlation(df_row = best_CpGs[r, ], typisation = real_age_typisation)
    ggsave(
      plot = p,
      filename = file.path(
        OUTPUT_FOLDER,
        "results",
        "pearson",
        folder_name,
        paste(best_CpGs$chrom[r],
              best_CpGs$start[r],
              ".png",
              sep = ".")
      )
    )
  }
}
# manual selection  --------------------------------------------------
if (FALSE) {
  test <- "pearson"
  all_data_complete_with_test <-
    readRDS(
      "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/results/WholeGenome/all_data_complete_with_test.rds"
    )
  
  p_significance_cut_off <- exp(-9.3)
  print(paste("p_significance_cut_off", p_significance_cut_off))
  print(paste(
    "-log(p_significance_cut_off)",-log(p_significance_cut_off)
  ))
  min_delta <- 0.39
  print(paste("min_delta ", min_delta))
  
  folder_name <- paste("delta",
                       min_delta,
                       "p_significance_cut_off",
                       p_significance_cut_off,
                       sep = "_")
  
  dir.create(
    path = file.path(OUTPUT_FOLDER, "results", test, folder_name),
    showWarnings = FALSE
  )
  
  
  significant_CpGs <-
    all_data_complete_with_test[all_data_complete_with_test$name != "NO_CHIP" &
                                  all_data_complete_with_test$pearson.delta >= min_delta  &
                                  all_data_complete_with_test$pearson.p_val <= p_significance_cut_off,]
  print(nrow(significant_CpGs))
  
  
  # Create a data frame with counts of each state
  state_counts <- table(significant_CpGs$state)
  state_counts_df <- data.frame(state = names(state_counts), count = as.numeric(state_counts))
  state_counts_df$state_color <- factor(x = df_universal_annotation$color,levels = unique(df_universal_annotation$color))
  state_counts_df$state_group <- factor(x =  df_universal_annotation$Group,levels = unique(df_universal_annotation$Group))
  state_counts_df$state <- factor(x = state_counts_df$state,levels = unique(state_counts_df$state))
  big_labels <- as.character(state_counts_df$state)
  big_labels[state_counts_df$count == 0] <-  " "
  state_counts_df$big_labels <- big_labels
  count_label <- as.character(state_counts_df$count)
  count_label[state_counts_df$count == 0] <-  " "
  # Create the pie chart with borders
  p<-
    ggplot(data = state_counts_df, aes(
      x = 1,
      y = count,
      fill = state,
      label = big_labels
    )) +
    coord_polar(theta = "y", start = 0) +  # Convert the bar plot into a pie chart
    labs(fill = "State", x = NULL) +   # Legend title
    theme_void() +  # Remove unnecessary elements
    #geom_text(mapping = aes(label = count))+
    geom_bar(width = 1,
             stat = "identity",
             color = "black", # Add borders to segments
             #mapping = aes(fill = state_group)
    ) +  
    geomtextpath::geom_textpath(
      data = state_counts_df,
      position = position_stack(vjust = 0.5),
      aes(
        x = 1.75,
      ),
      size = 3.4,
      #show.legend = FALSE
    ) +
    geomtextpath::geom_textpath(
      data = state_counts_df,
      position = position_stack(vjust = 0.5),
      aes(
        x = 1.2,
      ),
      size = 3.4,
      label = count_label
      #show.legend = FALSE
    )+
    #theme(legend.position = "none")+
    scale_fill_manual(values = as.character(df_universal_annotation$color))
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,folder_name,
      paste("pie_chart results", "png", sep = ".")
    )
    ,
    plot = p,
    width = 12,
    height = 12
  )
  
  bed4 <- significant_CpGs[, c("chrom", "start", "end", "name")]
  
  write.table(
    x = bed4,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      test,
      folder_name,
      "significant_CpGs.bed"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  bed_for_meme <- bed4 %>%
    dplyr::distinct(name, .keep_all = TRUE)
  
  add_nucleotides <- 25
  
  for (r in 1:length(unique(bed_for_meme$name))) {
    peak <- bed_for_meme$name[r]
    temp <- significant_CpGs[significant_CpGs$name == peak,]
    
    bed_for_meme$start[r] <- min(temp$start) - add_nucleotides
    bed_for_meme$end[r] <- max(temp$start) + add_nucleotides
  }
  
  bed_for_meme$score <- 1000.
  bed_for_meme$strand <- "."
  
  write.table(
    x = bed_for_meme[,],
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      test,
      folder_name,
      paste(
        "bed_for_meme",
        "add_nucleotides",
        add_nucleotides,
        "bed",
        sep = "."
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  ### make Genhancer query with files https://genome.ucsc.edu/cgi-bin/hgTables
  # C:\Users\Daniel Batyrev\Documents\GitHub\HumanEvo\HumanEvo\12.pipeline\results\pearson\delta_0.39_minLogP_9.3
  # GH Interactio double elite hg19
  
  # Database: hg19    Primary Table: geneHancerInteractionsDoubleElite Data last updated: 2019-01-15
  # Big Bed File Download: /gbdb/hg19/geneHancer/geneHancerInteractionsDoubleElite.v2.hg19.bb
  
  df_GH <- read.delim(
    file =  file.path(OUTPUT_FOLDER,
                      "results",
                      test,
                      folder_name,
                      "hgTables.txt"),
    header = FALSE,
    skip = 1,
    col.names = colnames(significant_CpGs)[1:5]
  )
  
  significant_CpGs$GeneAnnotation <- NA
  significant_CpGs$EnhancerAnnotation <- NA
  
  for (r in 1:nrow(significant_CpGs)) {
    # r <- 1
    annotated_GH <-
      df_GH[significant_CpGs$chrom[r] == df_GH$chrom &
              df_GH$start <= significant_CpGs$start[r] &
              significant_CpGs$end[r] <= df_GH$end, ]
    # if any annotation found
    if (nrow(annotated_GH) > 0) {
      # choose annotation with highest score
      annotation <-
        unlist(strsplit(annotated_GH$name[which.max(annotated_GH$score)], split = "/")) #annotated_GH$name[which.max(annotated_GH$score)]
      # get highest score
      significant_CpGs$score <-
        annotated_GH$score[which.max(annotated_GH$score)]
      significant_CpGs$GeneAnnotation[r] <- annotation[1]
      significant_CpGs$EnhancerAnnotation[r] <- annotation[2]
    }
  }
  
  print(paste("#of positive correlation CpGs : ",   nrow(significant_CpGs[significant_CpGs$pearson.statistic > 0 , ])))
  
  
  print(paste("#of H3K27AC peaks : ", length(unique(
    significant_CpGs$name
  ))))
  
  print(paste("#of annotated enhancers : ", length(
    unique(significant_CpGs$EnhancerAnnotation)
  )))
  print(paste("#of annotated genes : ", length(
    unique(significant_CpGs$GeneAnnotation)
  )))
  
  openxlsx::write.xlsx(
    x = significant_CpGs,
    file.path(
      OUTPUT_FOLDER,
      "results",
      test,
      folder_name,
      "significant_CpGs.xlsx"
    ),
    sheetName = "significant_CpGS",
    colNames = TRUE
  )
  
  
  write.csv2(
    x = unique(significant_CpGs$GeneAnnotation[!is.na(significant_CpGs$GeneAnnotation)]),
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      test,
      folder_name,
      paste(
        "unique.genes",
        min_delta,
        "delta.GH_Interactions",
        "p",
        p_significance_cut_off,
        "csv",
        sep = "."
      )
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  saveRDS(
    object = significant_CpGs,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      "pearson",
      folder_name,
      "significant_CpGs.rds"
    )
  )
  
  dir.create(
    path = file.path(OUTPUT_FOLDER,
                     "results",
                     test,
                     folder_name,
                     "plots"),
    showWarnings = FALSE
  )
  
  
  for (r in 1:nrow(significant_CpGs)) {
    p <-
      plot_age_correlation(df_row = significant_CpGs[r, ], typisation = real_age_typisation)
    ggsave(
      plot = p,
      filename = file.path(
        OUTPUT_FOLDER,
        "results",
        test,
        folder_name,
        "plots",
        paste(
          significant_CpGs$GeneAnnotation[r],
          significant_CpGs$state[r],
          significant_CpGs$chrom[r],
          significant_CpGs$start[r],
          ".png",
          sep = "."
        )
      )
    )
  }
}
# https://maayanlab.cloud/Enrichr/enrich?dataset=47ed0059768b38095dca92d0748c5055

# plot sim vs real --------------------------------------------------------
if (OTHER) {
  library(plotly)
  print("OTHER")
  plot_age_correlation(df_row = df_peak_CpG_complete_with_test[5, ], typisation = real_age_typisation)
  
  plot_food_correlation(df_row = df_peak_CpG_complete_with_test[5, ], typisation = real_food_typisation)
  
  search_df[which.min(search_df$FDR),]
  
  best_ex <-
    df_CpG[!is.na(df_CpG$KW_p_val) &
             df_CpG$type_delta > 0.25 &
             df_CpG$KW_p_val < exp(-9.3),]
  
  ex1 <- best_ex[1,]
  ex2 <- best_ex[2,]
  
  
  df_plot <- data.frame(matrix(nrow = N_SAMPLES, ncol = 0))
  df_plot$sample <- as.character(real_sample_typisation$sample)
  df_plot$meth <- as.numeric(ex1[real_sample_typisation$sample])
  df_plot$type <- real_sample_typisation$type
  
  ggplot(data = df_plot,
         mapping = aes(x = type, y = meth, color = type)) + ggtitle(paste(ex1$chrom, ex1$start)) +
    geom_boxplot() + geom_point(position = "jitter") +
    theme_minimal()
  
  
  df_plot <- data.frame(matrix(nrow = N_SAMPLES, ncol = 0))
  df_plot$sample <- as.character(real_sample_typisation$sample)
  df_plot$meth <- as.numeric(ex2[real_sample_typisation$sample])
  df_plot$type <- real_sample_typisation$type
  
  ggplot(data = df_plot,
         mapping = aes(x = type, y = meth, color = type)) + ggtitle(paste(ex2$chrom, ex2$start)) +
    geom_boxplot() + geom_point(position = "jitter") +
    theme_minimal()
  
  
  
  df_CpG[!is.na(df_CpG$KW_p_val) &
           df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3),]
  
  # min_delta minus_log_alpha n_signfincant_CpG real_count   FDR
  #    0.08             6.4              93.6        300 0.312
  
  length(unique(df_CpG$name[!is.na(df_CpG$KW_p_val) &
                              df_CpG$type_delta > 0.08 &
                              df_CpG$KW_p_val < exp(-6.4)]))
  
  # min_delta minus_log_alpha n_signfincant_CpG real_count   FDR
  #    0.25             9.3              0.17          2    0.085
  df_CpG[!is.na(df_CpG$KW_p_val) &
           df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3),]
  
  saveRDS(
    object = df_empirical_means,
    file = file.path(OUTPUT_FOLDER, "simulation", "empirical_means.rds")
  )
  
  # p <- ggplot(data = df_empirical_means,
  #             mapping = aes(x = min_log_p, y = real_count / count)) + geom_point() + theme_minimal() +
  #   ggtitle(min_delta, "min_delta")
}
# gene annotation --------------------------------------------------------
if (gene_annotation) {
  min_delta <- best_values$min_delta
  minus_log_p <- best_values$minus_log_alpha
  test <- best_values$test
  
  folder_name <- paste("delta",
                       min_delta,
                       "minLogP",
                       minus_log_p,
                       sep = "_")
  
  best_CpGs <- readRDS(
    object = ,
    file = file.path(OUTPUT_FOLDER,
                     "results",
                     test,
                     folder_name,
                     "best_CpGs.rds")
  )
  
  nrow(best_CpGs)
  
  
  bed4 <- best_CpGs[, c("chrom", "start", "end", "name")]
  
  write.table(
    x = bed4,
    file = file.path(OUTPUT_FOLDER,
                     "results",
                     test,
                     folder_name,
                     "best_CpGs.bed"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  bed_for_meme <- bed4 %>%
    dplyr::distinct(name, .keep_all = TRUE)
  
  add_nucleotides <- 25
  
  for (r in 1:length(unique(bed_for_meme$name))) {
    peak <- bed_for_meme$name[r]
    temp <- best_CpGs[best_CpGs$name == peak,]
    
    bed_for_meme$start[r] <- min(temp$start) - add_nucleotides
    bed_for_meme$end[r] <- max(temp$start) + add_nucleotides
  }
  
  bed_for_meme$score <- 1000.
  bed_for_meme$strand <- "."
  
  write.table(
    x = bed_for_meme[,],
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      test,
      folder_name,
      paste(
        "bed_for_meme",
        "add_nucleotides",
        add_nucleotides,
        "bed",
        sep = "."
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  ### make Genhancer query with files https://genome.ucsc.edu/cgi-bin/hgTables
  # C:\Users\Daniel Batyrev\Documents\GitHub\HumanEvo\HumanEvo\12.pipeline\results\pearson\delta_0.39_minLogP_9.3
  # GH Interactio double elite hg19
  
  # Database: hg19    Primary Table: geneHancerInteractionsDoubleElite Data last updated: 2019-01-15
  # Big Bed File Download: /gbdb/hg19/geneHancer/geneHancerInteractionsDoubleElite.v2.hg19.bb
  
  df_GH <- read.delim(
    file =  file.path(OUTPUT_FOLDER,
                      "results",
                      test,
                      folder_name,
                      "hgTables.txt"),
    header = FALSE,
    skip = 1,
    col.names = colnames(best_CpGs)[1:5]
  )
  
  best_CpGs$GeneAnnotation <- NA
  best_CpGs$EnhancerAnnotation <- NA
  
  for (r in 1:nrow(best_CpGs)) {
    # r <- 1
    annotated_GH <-
      df_GH[best_CpGs$chrom[r] == df_GH$chrom &
              df_GH$start <= best_CpGs$start[r] &
              best_CpGs$end[r] <= df_GH$end, ]
    # if any annotation found
    if (nrow(annotated_GH) > 0) {
      # choose annotation with highest score
      annotation <-
        unlist(strsplit(annotated_GH$name[which.max(annotated_GH$score)], split = "/")) #annotated_GH$name[which.max(annotated_GH$score)]
      # get highest score
      best_CpGs$score <-
        annotated_GH$score[which.max(annotated_GH$score)]
      best_CpGs$GeneAnnotation[r] <- annotation[1]
      best_CpGs$EnhancerAnnotation[r] <- annotation[2]
    }
  }
  
  print(paste("#of annotated enhancers : ", length(unique(
    best_CpGs$EnhancerAnnotation
  ))))
  print(paste("#of annotated genes : ", length(unique(
    best_CpGs$GeneAnnotation
  ))))
  
  openxlsx::write.xlsx(
    x = best_CpGs,
    file.path(
      OUTPUT_FOLDER,
      "results",
      test,
      folder_name,
      "best_CpGs.xlsx"
    ),
    sheetName = "significant_CpGS",
    colNames = TRUE
  )
  
  
  write.csv2(
    x = unique(best_CpGs$GeneAnnotation),
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      test,
      folder_name,
      paste(
        "unique.genes",
        min_delta,
        "delta.GH_Interactions",
        "p",
        minus_log_p,
        "csv",
        sep = "."
      )
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# PCA --------------------------------------------------------
if (plot_PCA) {
  normalize_PCA <- TRUE
  if (normalize_PCA) {
    normalized_string <- "normalized"
  } else{
    normalized_string <- ""
  }
  
  
  
  meta37 <-
    meta[as.character(meta$sample) %in% as.character(real_food_typisation$sample), ]
  # truncate to only three groups
  # trauncate min delta pearson.delta
  data_meth <-
    df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$pearson.delta >= 0.39, as.character(real_food_typisation$sample)]
  #print("erase non complete rows")
  data_meth <- data_meth[complete.cases(data_meth), ]
  data_meth <-
    t(data_meth)
  
  #print(sum(complete.cases(df_peak_CpG_complete_with_test)/nrow(data_meth)))
  #data_meth <- as.numeric(data_meth)
  
  data_meth_with_variance <-
    data_meth[, which(apply(data_meth, 2, var) != 0)]
  pca_prcomp_res <-
    prcomp(x = data_meth_with_variance, scale = normalize_PCA)
  
  pca_summery <- summary(pca_prcomp_res)$importance[2, ]
  
  saveRDS(object = pca_prcomp_res, file =  file.path(
    OUTPUT_FOLDER,
    paste(
      is_NORMALIZED,
      "pca_prcomp_res",
      "N_SAMPLES",
      "37",
      "rds",
      sep = "."
    )
  ))
  
  df_pca <- cbind(meta37, pca_prcomp_res$x)
  
  # sex
  p_pca <-
    ggplot(data = df_pca , aes(
      x = PC1,
      y = PC2,
      label = sample,
      color = sex
    )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_1v2 & sex") +
    theme_minimal() +
    xlab(paste("PC1 : ", round(pca_summery[1] * 100, 1), "%")) +
    ylab(paste("PC2 : ", round(pca_summery[2] * 100, 1), "%"))
  
  ggsave(
    plot = p_pca,
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "PCA",
      paste(normalized_string, "PCA_1v2_vs_sex_N_SAMPLES_37.png", sep = "")
    ),
    width = 8,
    height = 6
  )
  
  # sex
  p_pca <-
    ggplot(data = df_pca , aes(
      x = PC3,
      y = PC4,
      label = sample,
      color = sex
    )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_3v4 & sex") +
    theme_minimal() +
    xlab(paste("PC3 : ", round(pca_summery[3] * 100, 1), "%")) +
    ylab(paste("PC4 : ", round(pca_summery[4] * 100, 1), "%"))
  
  ggsave(
    plot = p_pca,
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "PCA",
      paste(normalized_string, "PCA_3v4_vs_sex_N_SAMPLES_37.png", sep = "")
    ),
    width = 8,
    height = 6
  )
  
  
  # locality
  # p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = Locality))+
  #   geom_point()+
  #   ggrepel::geom_text_repel()+
  #   labs(subtitle="PCA_1v2 & locality")
  #
  # ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_1v2_vs_locality_N_SAMPLES_37.png"),width=16, height=6)
  
  df_pca$old <-  df_pca$age_mean_BP >= 5700
  
  
  ##
  # "OLD"
  p_pca <-
    ggplot(data = df_pca , aes(
      x = PC1,
      y = PC2,
      label = sample,
      color = old
    )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_1v2 & OLD >= 5700Y") +
    theme_minimal() +
    xlab(paste("PC1 : ", round(pca_summery[1] * 100, 1), "%")) +
    ylab(paste("PC2 : ", round(pca_summery[2] * 100, 1), "%"))
  
  ##
  # "Country"
  p_pca <-
    ggplot(data = df_pca , aes(
      x = PC1,
      y = PC2,
      label = sample,
      color = Country
    )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_1v2 & Country") +
    theme_minimal() +
    xlab(paste("PC1 : ", round(pca_summery[1] * 100, 1), "%")) +
    ylab(paste("PC2 : ", round(pca_summery[2] * 100, 1), "%"))
  
  ggsave(
    plot = p_pca,
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "PCA",
      paste(normalized_string, "PCA_1v2_vs_Country.png", sep = "")
    ),
    width = 10,
    height = 6
  )
  
  p_pca <-
    ggplot(data = df_pca , aes(
      x = PC3,
      y = PC4,
      label = sample,
      color = Country
    )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_3v4 & Country") +
    theme_minimal() +
    xlab(paste("PC3 : ", round(pca_summery[3] * 100, 1), "%")) +
    ylab(paste("PC4 : ", round(pca_summery[4] * 100, 1), "%"))
  
  ggsave(
    plot = p_pca,
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "PCA",
      paste(normalized_string, "PCA_3v4_vs_Country.png", sep = "")
    ),
    width = 10,
    height = 6
  )
  
  
  ##
  # "age"
  p_pca <-
    ggplot(data = df_pca ,
           aes(
             x = PC1,
             y = PC2,
             label = sample,
             color = age_mean_BP > 5000
           )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_1v2 & age") +
    theme_minimal() +
    xlab(paste("PC1 : ", round(pca_summery[1] * 100, 1), "%")) +
    ylab(paste("PC2 : ", round(pca_summery[2] * 100, 1), "%")) +
    guides(color = guide_legend(title = "sample older then 5000 YBP"))
  
  ggsave(
    plot = p_pca,
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "PCA",
      paste(normalized_string, "PCA_1v2_vs_age.png", sep = "")
    ),
    width = 10,
    height = 6
  )
  
  p_pca <-
    ggplot(data = df_pca ,
           aes(
             x = PC3,
             y = PC4,
             label = sample,
             color = age_mean_BP > 5000
           )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_3v4 & age") +
    theme_minimal() +
    xlab(paste("PC3 : ", round(pca_summery[3] * 100, 1), "%")) +
    ylab(paste("PC4 : ", round(pca_summery[4] * 100, 1), "%")) +
    guides(color = guide_legend(title = "sample older then 5000 YBP"))
  
  ggsave(
    plot = p_pca,
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "PCA",
      paste(normalized_string, "PCA_3v4_vs_age.png", sep = "")
    ),
    width = 10,
    height = 6
  )
  
  ## Type
  
  p_pca <-
    ggplot(data = df_pca , aes(
      x = PC1,
      y = PC2,
      label = sample,
      color = Type
    )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_1v2 & Type") +
    theme_minimal() +
    xlab(paste("PC1 : ", round(pca_summery[1] * 100, 1), "%")) +
    ylab(paste("PC2 : ", round(pca_summery[2] * 100, 1), "%"))
  
  ggsave(
    plot = p_pca,
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "PCA",
      paste(normalized_string, "PCA_1v2_vs_Type.png", sep = "")
    ),
    width = 10,
    height = 6
  )
  
  p_pca <-
    ggplot(data = df_pca , aes(
      x = PC3,
      y = PC4,
      label = sample,
      color = Type
    )) +
    geom_point() +
    ggrepel::geom_text_repel() +
    #labs(subtitle = "PCA_3v4 & Type") +
    theme_minimal() +
    xlab(paste("PC3 : ", round(pca_summery[3] * 100, 1), "%")) +
    ylab(paste("PC4 : ", round(pca_summery[4] * 100, 1), "%"))
  
  ggsave(
    plot = p_pca,
    filename = file.path(
      OUTPUT_FOLDER,
      "plots",
      "PCA",
      paste(normalized_string, "PCA_3v4_vs_Type.png", sep = "")
    ),
    width = 10,
    height = 6
  )
}


# Enrichment analyis results ---------------------------------------------------
if (DO_enrichment_analysis) {
  enrichment_folder_name <-
    file.path(OUTPUT_FOLDER,
              "results",
              "pearson",
              folder_name,
              "enrichment analysis")
  
  library(openxlsx)
  
  # Specify the file path
  file_path <-
    file.path(enrichment_folder_name, "Enrichr_Submissions.xlsx")
  
  # Create a workbook
  wb <- createWorkbook()
  
  # read data
  TF_Gene_Coocurrence <-
    read.delim(
      file = file.path(
        enrichment_folder_name,
        "Enrichr_Submissions_TF-Gene_Coocurrence_table.txt"
      )
    )
  #TF_regulation_text_mining <- read.delim(file = file.path(enrichment_folder_name,"TRRUST_Transcription_Factors_2019_table.txt"))
  ARCHS4_TFs_Coexp <-
    read.delim(file = file.path(enrichment_folder_name, "ARCHS4_TFs_Coexp_table.txt"))
  
  #MGI_Mammalian_Phenotype_Level_4_2021_table.txt
  MGI_Mammalian_Phenotype_Level_4_2021 <-
    read.delim(
      file = file.path(
        enrichment_folder_name,
        "MGI_Mammalian_Phenotype_Level_4_2021_table.txt"
      )
    )
  
  # Add a worksheet
  addWorksheet(wb, sheetName = "TF-Gene_Coocurrence")
  addWorksheet(wb, sheetName = "ARCHS4_TFs_Coexp")
  addWorksheet(wb, sheetName = "Mammalian_Phenotype")
  
  # Write data to the worksheet
  writeData(wb, sheet = 1, x = TF_Gene_Coocurrence[, c(
    "Term",
    "Overlap",
    "P.value",
    "Adjusted.P.value",
    "Odds.Ratio",
    "Combined.Score",
    "Genes"
  )])
  writeData(wb, sheet = 2, x = ARCHS4_TFs_Coexp[, c(
    "Term",
    "Overlap",
    "P.value",
    "Adjusted.P.value",
    "Odds.Ratio",
    "Combined.Score",
    "Genes"
  )])
  
  writeData(wb, sheet = 3, x = MGI_Mammalian_Phenotype_Level_4_2021[, c(
    "Term",
    "Overlap",
    "P.value",
    "Adjusted.P.value",
    "Odds.Ratio",
    "Combined.Score",
    "Genes"
  )])
  
  # Save the workbook to an Excel file
  saveWorkbook(wb, file_path, overwrite = TRUE)
}


# test
#
#
#
# df <- data.frame(SIRT1 = as.numeric(best_CpGs[!is.na(best_CpGs$GeneAnnotation) & best_CpGs$GeneAnnotation == "SIRT1",as.character(meta$sample)]),
#                  FMO5 = as.numeric(best_CpGs[!is.na(best_CpGs$GeneAnnotation) & best_CpGs$GeneAnnotation == "FMO5",as.character(meta$sample)]))
#
#
# # Calculate Pearson correlation
# corr_result <- cor.test(df$SIRT1, df$FMO5, method = "pearson")
#
# # Extract correlation coefficient and p-value
# corr_coef <- corr_result$estimate
# p_value <- corr_result$p.value
#
# # Create the plot
# ggplot(data = df, aes(x = SIRT1, y = FMO5)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) +
#   annotate("text", x = Inf, y = Inf, label = sprintf("r = %.2f (p = %.6f)", corr_coef, p_value),
#            hjust = 1.1, vjust = 2, size = 5, colour = "red")
