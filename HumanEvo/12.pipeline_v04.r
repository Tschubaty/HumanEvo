# Documentation -----------------------------------------------------------
##
##
##  input: methylation+chip
##
##
##  output: pipeline
##  # histogramm matched updat: all_CpG.45.samples.histogram.merged2025.hg19.rds
##  
##  v_04 10.04.2025
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
library(dplyr)
library(data.table)

# set constants  ----------------------------------------------------------

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

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
df_universal_annotation$color[df_universal_annotation$Group == "weak transcription"] <-
  "#009600"
df_universal_annotation$color[df_universal_annotation$Group == "transcription"] <-
  "#008000"

df_universal_annotation$Group[df_universal_annotation$Group == "quescient"] <-
  "quiescent"
df_universal_annotation$Group[df_universal_annotation$Group == "HET"] <-
  "heterochromatin"

state_names <-
  paste(
    df_universal_annotation$state_order_by_group,
    df_universal_annotation$`States  presented in paper`,
    sep = "_"
  )

# matching_indices <- match(df_state_sample$state,paste(df_universal_annotation$state_order_by_group,df_universal_annotation$`States  presented in paper`,sep = "_"))
# df_state_sample$state_color <- df_universal_annotation$color[matching_indices]

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

## only big sample groups
TYPES <- c("Farmer", "Steppe", "HG","Modern")

meta <- as.data.frame(readxl::read_xlsx("AGDP.45.metadata3.xlsx"))
meta <- meta[order(meta$age_mean_BP), ]
meta$sample <- factor(x = meta$sample, levels = meta$sample)

# 1. Compute mean age per Type
type_order <- meta %>%
  group_by(Type) %>%
  summarise(mean_age = mean(age_mean_BP, na.rm = TRUE)) %>%
  arrange(mean_age)

# 2. Reorder factor levels by ascending mean age
meta$Type <- factor(meta$Type, levels = type_order$Type)

real_age_typisation <-
  meta[order(meta$age_mean_BP), c("sample", "age_mean_BP","Type")]


real_food_typisation <-
  meta[meta$Type %in% TYPES, c("sample", "Type")]

# Compute the mean of age_mean_BP by Type
mean_age_by_type <-
  aggregate(meta$age_mean_BP,
            by = list(Type = meta$Type),
            FUN = mean)
# Rename the columns for clarity
colnames(mean_age_by_type) <- c("Type", "Mean_age_mean_BP")
mean_age_by_type <-
  mean_age_by_type[order(mean_age_by_type$Mean_age_mean_BP), ]
mean_age_by_type$color <- gg_color_hue(nrow(mean_age_by_type))

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
compute_df_state_sample <- FALSE
compute_df_state_group_sample <- FALSE
compute_df_state <- FALSE
compute_df_state_group <- FALSE
compute_make_plots <- FALSE
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
Genome_wide_Analysis <- FALSE
state_dependen_analysis <- FALSE
plot_landscape <- FALSE

# load Data ---------------------------------------------------------------
# load(".RData")
df_peak_CpG <-
  readRDS(
    "methylation+chip/H3K27ac.only.45.samples.merged.hg19.rds"
  )

# define functions ---------------------------------------------------------


correct_state_name <- function(x) ifelse(x == "quescient", "quiescent", ifelse(x == "HET", "heterochromatin", x))


# Create a nice printout
interpretation <- function(test_name, test_result, subset_name, fullset_name) {
  cat("\n")
  cat(paste0("=====================================\n"))
  cat(paste0("          ", test_name, " Results\n"))
  cat(paste0("=====================================\n"))
  cat(paste0("Test Statistic: ", formatC(test_result$statistic, format = "e", digits = 2), "\n"))
  cat(paste0("P-value: ", formatC(test_result$p.value, format = "e", digits = 2), "\n"))
  cat(paste0("\nInterpretation:\n"))
  
  if (test_result$p.value < 0.05) {
    cat(paste0("There is a statistically significant difference in the distribution of ", subset_name, " compared to ", fullset_name, ".\n"))
    cat(paste0("The p-values in ", subset_name, " are significantly smaller than those in ", fullset_name, ".\n"))
  } else {
    cat(paste0("There is no statistically significant difference in the distribution of ", subset_name, " compared to ", fullset_name, ".\n"))
    cat(paste0("The p-values in ", subset_name, " are not significantly smaller than those in ", fullset_name, ".\n"))
  }
  cat("\n")
}

# Define the main function
test_genomic_subsets <- function(df_ranges, df_full) {
  # Ensure unique rows
  df_ranges <- df_ranges %>% distinct()
  
  # Convert data frames to GRanges objects
  gr_ranges <- makeGRangesFromDataFrame(df_ranges, seqnames.field = "chr", start.field = "start", end.field = "end")
  gr_full <- makeGRangesFromDataFrame(df_full, seqnames.field = "chrom", start.field = "start", end.field = "end")
  
  # Find overlaps
  overlaps <- findOverlaps(gr_full, gr_ranges)
  
  # Select rows with overlaps
  df_subset <- df_full[queryHits(overlaps), ]
  
  cat(paste0("# CpGs in subset query: ", as.character(nrow(df_subset)), ".\n"))
  
  # Ensure the pearson.p_val and kw.p_val columns exist in both datasets
  if (!("pearson.p_val" %in% names(df_subset) && "pearson.p_val" %in% names(df_full))) {
    stop("pearson.p_val column not found in one of the datasets.")
  }
  if (!("kw.p_val" %in% names(df_subset) && "kw.p_val" %in% names(df_full))) {
    stop("kw.p_val column not found in one of the datasets.")
  }
  

  print(nrow(df_subset))
  print(nrow(df_full))
  
  # Perform Wilcoxon rank-sum test for pearson.p_val
  pearson_test_result <- wilcox.test(df_subset$pearson.p_val, df_full$pearson.p_val, alternative = "less")
  
  # Print the result for pearson.p_val
  interpretation("Wilcoxon Rank-Sum Test for p-values from Pearson correlation", pearson_test_result, "subset", "full dataset")
  
  # Perform Wilcoxon rank-sum test for kw.p_val
  kw_test_result <- wilcox.test(df_subset$kw.p_val, df_full$kw.p_val, alternative = "less")
  
  # Print the result for kw.p_val
  interpretation("Wilcoxon Rank-Sum Test for p-values from KW test", kw_test_result, "subset", "full dataset")
  return(list(pearson_test_result = pearson_test_result,kw_test_result = kw_test_result))
}

#' Load a Variable if Not Already Present in Workspace
#'
#' This function checks if a specified variable is already loaded in the global environment. 
#' If the variable does not exist, it loads the variable from a specified RDS file. This function
#' is particularly useful in managing workspace memory and ensuring data is loaded only when needed.
#'
#' @param variable_name The name of the variable as a string. This is the name of the variable
#'   that will be checked and possibly loaded into the global environment.
#' @param file_path The full file path to the RDS file where the variable data is stored. This path
#'   must be correct and the file must exist; otherwise, an error will occur during the loading process.
#'
#' @return No explicit return value. The function works by side effect, loading a variable into the 
#'   global environment if it does not already exist. Messages will be printed to indicate whether 
#'   the variable was loaded or if it was already present.
#'
#' @examples
#' # Define the path to your RDS file
#' file_path <- "path/to/your/file.rds"
#' # Define the variable name you want to check/load
#' variable_name <- "desired_variable"
#' # Call the function
#' load_variable_if_not_exists(variable_name, file_path)
#'
#' @importFrom utils readRDS
#' @export
load_variable_if_not_exists <- function(variable_name, file_path) {
  if (!exists(variable_name, where = .GlobalEnv)) {
    assign(variable_name, readRDS(file_path), envir = .GlobalEnv)
    message("Variable '", variable_name, "' has been loaded into the workspace.")
  } else {
    message("Variable '", variable_name, "' is already loaded in the workspace.")
  }
}


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
    meth <-
      as.numeric(df_row[as.character(typisation$sample[typisation$Type == s])])
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
  clus <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
  # get meth values
  df_meth <-
    df[, as.character(age.typisation$sample[order(age.typisation$age_mean_BP)])]
  # get numeric  age vector
  age <-
    age.typisation$age_mean_BP[order(age.typisation$age_mean_BP)]
  
  # Export  to workspace
  clusterExport(
    clus, 
    varlist = c("age.typisation", "age"),
    envir   = environment()  # or parent.env(environment()), depending on setup
  )
  
  
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
  pearson.p_val <- sapply(psr_cor_tests, function(x) {
    if (inherits(x, "htest")) {
      # x is a cor.test object, so extract its p-value
      x$p.value
    } else {
      # x was NA or something else, so return NA
      NA
    }
  })
  
  
  # extract direction
  print("pearson stat")
  pearson.statistic <-
    sapply(
      psr_cor_tests ,
      FUN = function(r) {
        if (inherits(r, "htest")) {
          return(r$statistic)
        } else{
          return(NA)
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


#' Computes Pearson correlation test on data.frame of methylation
#' Function is parallelized
#' @param df data frame with all methylation values
#' with the appropriate sample column names.
#' @param age.typisation dataframe storing sample names in $sample
#' and age in $age_mean_BP; one can make permutations
#' by permuting $age_mean_BP or $sample in the age.typisation parameter
#' @returns A data.frame with the following columns: pearson.p_val pearson.statistic pearson.delta
parallel_testing_pearson_cor_new <- function(df, age.typisation) {
  clus <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
  
  # Order samples by age
  sample_order <- order(age.typisation$age_mean_BP)
  df_meth <- df[, as.character(age.typisation$sample[sample_order])]
  age <- age.typisation$age_mean_BP[sample_order]
  
  # Export variables
  parallel::clusterExport(clus, varlist = c("df_meth", "age"), envir = environment())
  
  print("Running parallel Pearson correlation (robust)")
  
  results_list <- parallel::parLapply(clus, X = seq_len(nrow(df_meth)), fun = function(i) {
    r <- df_meth[i, ]
    valid <- !is.na(r) & !is.na(age)
    n <- sum(valid)
    
    if (n < 3) {
      return(list(pval = NA_real_, stat = NA_real_, delta = NA_real_))
    }
    
    r_val <- suppressWarnings(cor(r, age, use = "complete.obs"))
    
    if (is.na(r_val) || abs(r_val) == 1) {
      return(list(pval = NA_real_, stat = NA_real_, delta = NA_real_))
    }
    
    t_stat <- r_val * sqrt((n - 2) / (1 - r_val^2))
    p_val <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
    delta <- max(r, na.rm = TRUE) - min(r, na.rm = TRUE)
    
    return(list(pval = p_val, stat = t_stat, delta = delta))
  })
  
  stopCluster(clus)
  gc()
  
  # Convert list to data.frame
  df_x <- do.call(rbind, lapply(results_list, as.data.frame))
  colnames(df_x) <- c("pearson.p_val", "pearson.statistic", "pearson.delta")
  
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
        df$delta > min_delta,
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
        df$delta > min_delta &
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
        df$delta > min_delta &
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
                                   return(sum(df$delta > min_delta))
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
  df <- df[, !colnames(df) %in%
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
  df_plot$meth <- as.numeric(df_row[, as.character(typisation$sample)])
  
  # Fit a linear model to calculate R-squared
  lm_model <- lm(meth ~ age_mean_BP, data = df_plot)
  r_squared <- summary(lm_model)$r.squared  # Extract R-squared value
  
  p <- ggplot(data = df_plot, mapping = aes(x = age_mean_BP, y = meth)) +
    ggtitle(
      label = paste(df_row$chrom, ":",df_row$start),  # Improved readability
      subtitle = bquote(
        "P-value:" ~ .(format.pval(df_row$pearson.p_val, digits = 3)) ~ 
          ", " ~ italic(R)^2 ~ ":" ~ .(format(r_squared, digits = 3))  # Adding R-squared to subtitle
      )
    ) +
    geom_point(color = "#1f77b4", size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", color = "#ff7f0e", se = FALSE, size = 1) +
    xlim(0, 11000) +  
    ylim(0, 1) +  
    xlab("Sample age estimation (years BP)") +  # More concise label
    ylab("CpG methylation level") +  # More specific label
    scale_x_continuous(breaks = seq(0, 10000, by = 2000)) +  # Set x-axis ticks to match data range
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14, color = "gray40"),  # Subtle color for subtitle
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      panel.grid.minor = element_blank(),  # Reduce clutter from minor grid lines
      panel.grid.major.x = element_line(color = "gray80"),  # Subtle major grid lines
      panel.grid.major.y = element_line(color = "gray80")
    )
  
  return(p)
}


#' Create a linear correlation plot (methylation vs. age) for one CpG row
#'
#' @param df_row A row in the CpG dataframe (with columns matching \code{typisation$sample})
#' @param typisation A data frame storing sample metadata:
#'   - \code{sample}
#'   - \code{age_mean_BP}
#'   - \code{Type} (categorical variable to color points)
#' @returns A ggplot2 object with a scatter plot + linear regression line
plot_age_correlation_type <- function(df_row, typisation) {
  # Prepare plotting data: attach methylation for the relevant row
  df_plot <- typisation
  df_plot$meth <- as.numeric(df_row[, as.character(typisation$sample)])
  
  # Fit a linear model to calculate R-squared
  lm_model <- lm(meth ~ age_mean_BP, data = df_plot)
  r_squared <- summary(lm_model)$r.squared
  
  # Build the ggplot
  p <- ggplot(data = df_plot, 
              aes(x = age_mean_BP, y = meth, color = Type)) +  # <-- COLOR BY Type
    ggtitle(
      label = paste(df_row$chrom, ":", df_row$start),  # Title with chromosomal location
      subtitle = bquote(
        "P-value:" ~ .(format.pval(df_row$pearson.p_val, digits = 3)) * 
          ", " * italic(R)^2 * ":" * .(format(r_squared, digits = 3))
      )
    ) +
    geom_point(size = 3, alpha = 0.6) +  # Removed color from geom_point, replaced with color = Type in aes()
    geom_smooth(method = "lm", color = "#ff7f0e", se = FALSE, size = 1) +
    xlim(0, 11000) +  
    ylim(0, 1) +  
    xlab("Sample age estimation (years BP)") +
    ylab("CpG methylation level") +
    scale_x_continuous(breaks = seq(0, 10000, by = 2000)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14, color = "gray40"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray80"),
      panel.grid.major.y = element_line(color = "gray80")
    )
  
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
  
  # p <- ggplot(data = df_plot,
  #             mapping = aes(x = Type, y = meth, color = Type)) +
  #   ggtitle(
  #     label = paste(df_row$chrom, df_row$start),
  #     subtitle = paste("p_val: ", format.pval(df_row$kw.p_val))
  #   ) +
  #   geom_boxplot() +
  #   geom_point(position = "jitter") +
  #   ylim(0, 1) +
  #   xlab("sample food Type") +
  #   ylab("CpG methylation") +
  #   scale_color_manual(values = setNames(mean_age_by_type$color, mean_age_by_type$Type))+
  #   theme_minimal()

  p <- ggplot(data = df_plot, mapping = aes(x = Type, y = meth, color = Type)) +
    ggtitle(
      label = paste(df_row$chrom, df_row$start),
      subtitle = paste("p_val:", signif(df_row$kw.p_val, digits=3))
    ) +
    geom_boxplot(outlier.shape = NA)  +
    geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.6) +
    ylim(0, 1) +
    ylab("CpG Methylation") +
    xlab("") +
    scale_color_manual(values = setNames(mean_age_by_type$color, mean_age_by_type$Type)) +
    theme_minimal(base_size = 14) +  # Increase base font size for better readability
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_blank(),  # Remove legend title
      legend.text = element_text(size = 12),
      legend.position = "right",
      plot.subtitle = element_text(face = "italic")
    )
  return(p)
}


# done H3K27ac Analysis --------------------------------------------------------
if (H3K27ac_Analysis) {
  
  load_variable_if_not_exists(
    variable_name = "all_CpG.45.samples.merged.hg19",
    file_path = "methylation+chip/all_CpG.45.samples.merged.hg19.rds"
  )
  
  load_variable_if_not_exists(
    variable_name = "df_state_group",
    file_path = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("df_state_group", "rds", sep = ".")
    )
  )
  

  
  # create results folder
  folder_name <- file.path(OUTPUT_FOLDER, "results", "H3K27ac")
  if (!file.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }
  
  all_CpG.meth <-
    all_CpG.45.samples.merged.hg19[, as.character(meta$sample)]
  
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
      all_CpG.45.samples.merged.hg19[all_CpG.45.samples.merged.hg19$chrom == CHR_NAMES[n_chr], c(1:6)]
    # load state annotation
    chr_chromatin_seg <-
      readRDS(
        paste(
          "Chromatin/processed/",
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
    factor(df_universal_annotation$Group[matching_indices],
           levels = unique(df_universal_annotation$Group))
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
  
  # Figure S2
  
  load_variable_if_not_exists(
    variable_name = "all_CpG_states_hg19",
    file_path =       file_name <- file.path(
      OUTPUT_FOLDER,
      "results",
      "H3K27ac",
      "all_CpG_states_hg19.rds"
    )
  )
  
  color_values <- unique(df_universal_annotation$color)
  color_labels <- unique(df_universal_annotation$Group)
  
  # Set the first color to gray
  color_values[1] <- "gray"
  
  change_indices <- which(c(TRUE, df_universal_annotation$color[-1] != df_universal_annotation$color[-length(df_universal_annotation$color)])) -1 
  # Prepare a vector of labels, initially setting all to empty strings
  labels <- rep("", length(df_universal_annotation$color))
  # Set labels at change indices to the corresponding state or number, you can customize it further if needed
  labels[change_indices+1] <- as.character(change_indices)
  
  
  # Prepare the plot
  p <- ggplot(data = all_CpG_states_hg19, aes(x = state, fill = state_color)) +
    geom_bar(stat = "count", color = "black") +
    theme(
      axis.line = element_line(colour = "black"),
      panel.background = element_blank(),
      axis.text.x = element_text(
        angle = 45, vjust = 1, hjust = 1, size = 10, face = "bold"
      ),
      text = element_text(size = 14),  # Adjust size as needed
      axis.text = element_text(face = "bold"),  # Make all axis text bold
      axis.title = element_text(face = "bold"),  # Make axis titles bold
      legend.text = element_text(face = "bold"),  # Make legend text bold
      legend.title = element_blank()  # Remove the legend title
    ) +
    ylab("# of CpGs in chromatin state") +
    xlab("") +  # Empty x-label as per your preference
    scale_fill_manual(
      values = color_values,
      labels = color_labels
    ) +
    scale_x_discrete(
      labels = labels)  # Labels from 0 to 99

  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "results",
                         "H3K27ac",
                         "all_CpG_states.png")
    ,
    plot = p,
    width = 18,
    height = 10
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
          "Chromatin/processed/",
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
  
  state_segementation_meta$state_group <- 
  
  p <-
    ggplot(
      data = state_segementation_meta,
      mapping = aes(x = state, y = state_fraction_is_H3K27ac * 100, fill = state)
    ) +
    geom_bar(stat = "identity") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    ) +
    ylim(0, 100) +
    ylab("% of CpGs covered by H3K27ac peaks in bone sample") +
    xlab("genome segment chromatin state labeling")
  
  # Figure S6
  
  p <- ggplot(data = df_state, aes(x = state, y = n_H3K27AC / n_total * 100, fill = state_group)) +
    geom_bar(stat = "identity", color = "black") +
    theme(
      axis.line = element_line(colour = "black"),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, face = "bold"),  # X-axis labels bold
      axis.text.y = element_text(face = "bold"),  # Y-axis labels bold
      axis.title.x = element_text(face = "bold"),  # X-axis title bold
      axis.title.y = element_text(face = "bold", size = 14),  # Y-axis title bold and potentially adjust size
      legend.title = element_text(face = "bold"),  # Legend title bold
      legend.text = element_text(face = "bold")  # Legend items text bold
    ) +
    ylim(0, 100) +
    ylab("% of CpGs in genome covered by H3K27ac") +
    xlab("") +
    scale_fill_manual(
      values = color_values,
      labels = color_labels
    ) +
    scale_x_discrete(
      labels = labels)  # Labels from 0 to 99

  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "results",
                         "H3K27ac",
                         "H3K27ac_coverage.png")
    ,
    plot = p,
    width = 16,
    height = 8
  )
  
  # Figure 2A
  p <-
    ggplot(
      data = df_state_group,
      mapping = aes(
        x = state_group,
        y = n_H3K27AC / n_total * 100,
        fill = state_group
      )
    ) +
    geom_bar(stat = "identity", color = "black") +
    theme(
      axis.line = element_line(colour = "black"),
      panel.background = element_blank(),
      #legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 12),  # Make x-axis text bold
      axis.text.y = element_text(face = "bold", size = 12),  # Make y-axis text bold
      axis.title.x = element_text(face = "bold", size = 14),  # Make x-axis title bold
      axis.title.y = element_text(face = "bold", size = 14),  # Make y-axis title bold
      plot.title = element_text(face = "bold", size = 16),  # If you have a title, make it bold
      plot.subtitle = element_text(face = "bold"),  # If you have a subtitle, make it bold
      text = element_text(face = "bold")  # Make all other text bold by default
    ) +
    ylim(0, 40) +
    ylab("% of CpGs covered by H3K27ac peaks") +
    xlab("") +
    scale_fill_manual(values = unique(as.character(df_universal_annotation$color)))+
    #scale_fill_manual(values = "black")+
    scale_x_discrete(labels = function(x) ifelse(x == "quescient", "quiescent", ifelse(x == "HET", "heterochromatin", x)))+
    geom_hline(yintercept = 2.4, linetype = "dashed", color = "blue")
  
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "results",
                         "H3K27ac",
                         "15_state_group_H3K27ac_coverage.png")
    ,
    plot = p,
    width = 6,
    height = 10
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
# done Genome wide Analysis --------------------------------------------------------
if (Genome_wide_Analysis) {
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
      paste(CHR_NAMES[n_chr], "complete_with_test","45", "rds", sep = ".")
    )
    # if (file.exists(file_name)) {
    #   next
    # }
    
    # cut out only coordinates per chromosome
    df_chr <-
      all_CpG.45.samples.merged.hg19[all_CpG.45.samples.merged.hg19$chrom == CHR_NAMES[n_chr], ]
    
    # compute sample numbers per CpG position
    start_time <- Sys.time()
    df_chr$score <-
      parallel_score_sample_number(df_chr, real_age_typisation)
    end_time <- Sys.time()
    print(end_time - start_time)
    print(paste("# CpG in chr:", nrow(df_chr)))

    print(paste("# CpG in chr with all samples :",
                sum(
                  df_chr$score == nrow(real_age_typisation)
                )))

    print(paste("discard ", 100 * (
      1 - sum(df_chr$score == nrow(real_age_typisation)) / nrow(df_chr)
    ), "% of CpG positions"))

    # discard CpG with not full data
    df_chr_complete <-
      df_chr[nrow(real_age_typisation) == df_chr$score, ]
    
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
# done compute_all_CpG_complete_with_test.45 -all_CpG_complete_with_test.45 -----------------------------------------------------------------------
if (FALSE) {
  # compine all chr data to big df
  all_CpG_complete_with_test.45 <- data.frame()
  for (n_chr in 1:length(CHR_NAMES)) {
    print(n_chr)
    file_name <- file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste(CHR_NAMES[n_chr], "complete_with_test","45", "rds", sep = ".")
    )
    df_chr <- readRDS(file = file_name)
    all_CpG_complete_with_test.45 <-
      rbind(all_CpG_complete_with_test.45, df_chr)
  }
  
  all_CpG_complete_with_test.45$state <-
    factor(x = all_CpG_complete_with_test.45$state,
           levels = state_names)
  
  matching_indices <-
    match(all_CpG_complete_with_test.45$state,
          state_names)
  all_CpG_complete_with_test.45$state_color <-
    factor(df_universal_annotation$color[matching_indices],
           levels = unique(df_universal_annotation$color))
  all_CpG_complete_with_test.45$state_group <-
    factor(df_universal_annotation$Group[matching_indices],
           levels = correct_state_name(x = unique(df_universal_annotation$Group)))

  
  all_CpG_complete_with_test.45$pearson.p_val_BH <- p.adjust(all_CpG_complete_with_test.45$pearson.p_val, method = "BH")
  
  saveRDS(
    object = all_CpG_complete_with_test.45,
    file =  file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("all_CpG_complete_with_test.45", "rds", sep = ".")
    )
  )
}
# done df_state_sample  -----------------------------------------------------------------------
if (compute_df_state_sample) {
  
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      "all_CpG_complete_with_test.45.rds" 
    )
  )
  
  #df_state_sample <- summary_df <- aggregate(. ~ state, data =  all_CpG_complete_with_test.45[,c("state",as.character(meta$sample))], FUN = mean)
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
  df_state_sample <-
    data.frame(matrix(
      NA,
      nrow = length(state_names) * nrow(meta),
      ncol = length(state_summery_column_names)
    ))
  
  # Rename the columns
  colnames(df_state_sample) <- state_summery_column_names
  
  r <- 1
  for (state in  state_names) {
    print(state)
    for (sample in as.character(meta$sample)) {
      df_state_sample$state[r] <- state
      df_state_sample$sample[r] <- sample
      
      index <- all_CpG_complete_with_test.45$state == state
      meth_column <- all_CpG_complete_with_test.45[, sample]
      meth <- meth_column[index]
      
      df_state_sample$mean_methylation[r] <- mean(meth)
      df_state_sample$sd_methylation[r] <- sd(meth)
      
      df_state_sample$methylation_H3K27AC[r] <-
        mean(meth_column[index &
                           all_CpG_complete_with_test.45$name != "NO_CHIP"])
      df_state_sample$methylation_NO_CHIP[r] <-
        mean(meth_column[index &
                           all_CpG_complete_with_test.45$name == "NO_CHIP"])
      
      r <- r + 1
    }
  }
  
  # change state to factor
  df_state_sample$state <-
    factor(x = df_state_sample$state, levels = state_names)
  
  # annotate states
  matching_indices <-
    match(
      df_state_sample$state,
      paste(
        df_universal_annotation$state_order_by_group,
        df_universal_annotation$`States  presented in paper`,
        sep = "_"
      )
    )
  df_state_sample$state_color <-
    df_universal_annotation$color[matching_indices]
  df_state_sample$state_group <-
    factor(df_universal_annotation$Group[matching_indices],
           levels = unique(df_universal_annotation$Group))
  # annotate type
  matching_indices <- match(df_state_sample$sample, meta$sample)
  df_state_sample$type <- meta$Type[matching_indices]
  
  # match color
  matching_indices <-
    match(df_state_sample$type, mean_age_by_type$Type)
  df_state_sample$type_color <-
    mean_age_by_type$color[matching_indices]
  df_state_sample$type <-
    factor(x = df_state_sample$type , levels = mean_age_by_type$Type)
  
  saveRDS(
    object = df_state_sample,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("df_state_sample", "rds", sep = ".")
    )
  )
}
# done df_state_group_sample  -----------------------------------------------------------------------
if (compute_df_state_group_sample) {
  
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("all_CpG_complete_with_test.45", "rds", sep = ".")
    )
  )
  
  state_group_sample_summery_column_names <-
    c(
      "state_group",
      "sample",
      "mean_methylation",
      "sd_methylation",
      "methylation_H3K27AC",
      "methylation_NO_CHIP"
    )
  
  # Create an empty data frame with the specified column names
  df_state_group_sample <-
    data.frame(matrix(
      NA,
      nrow = length(unique(df_universal_annotation$Group)) * nrow(meta),
      ncol = length(state_group_sample_summery_column_names)
    ))
  
  # Rename the columns
  colnames(df_state_group_sample) <- state_group_sample_summery_column_names
  
  r <- 1
  for (sg in unique(df_universal_annotation$Group)) {
    print(sg)
    for (sample in as.character(meta$sample)) {
      print(sample)
      df_state_group_sample$state_group[r] <- sg
      df_state_group_sample$sample[r] <- sample
      
      index <- all_CpG_complete_with_test.45$state_group == sg
      meth_column <- all_CpG_complete_with_test.45[, sample]
      meth <- meth_column[index]
      
      df_state_group_sample$mean_methylation[r] <- mean(meth)
      df_state_group_sample$sd_methylation[r] <- sd(meth)
      
      df_state_group_sample$methylation_H3K27AC[r] <-
        mean(meth_column[index &
                           all_CpG_complete_with_test.45$name != "NO_CHIP"])
      df_state_group_sample$methylation_NO_CHIP[r] <-
        mean(meth_column[index &
                           all_CpG_complete_with_test.45$name == "NO_CHIP"])
      
      r <- r + 1
    }
  }
  
  # change state to factor
  df_state_group_sample$state_group <-
    factor(x = df_state_group_sample$state_group, levels = unique(df_universal_annotation$Group))
  
  # annotate states
  matching_indices <-
    match(
      df_state_group_sample$state_group,
      unique(df_universal_annotation$Group)
    )
  df_state_group_sample$state_group_color <-
    df_universal_annotation$color[matching_indices]
 
  # annotate type
  matching_indices <- match(df_state_group_sample$sample, meta$sample)
  df_state_group_sample$type <- meta$Type[matching_indices]
  
  # match color
  matching_indices <-
    match(df_state_group_sample$type, mean_age_by_type$Type)
  df_state_group_sample$type_color <-
    mean_age_by_type$color[matching_indices]
  df_state_group_sample$type <-
    factor(x = df_state_group_sample$type , levels = mean_age_by_type$Type)
  
  saveRDS(
    object = df_state_group_sample,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("df_state_group_sample", "rds", sep = ".")
    )
  )
}
# done df_state  -----------------------------------------------------------------------
if (compute_df_state) {
  
    load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("all_CpG_complete_with_test.45", "rds", sep = ".")
    )
  )
  
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
  df_state <-
    data.frame(matrix(
      NA,
      nrow = length(state_names),
      ncol = length(state_column_names)
    ))
  
  # Rename the columns
  colnames(df_state) <- state_column_names
  
  r <- 1
  for (state in  state_names) {
    print(state)
    
    df_state$state[r] <- state
    
    index <- all_CpG_complete_with_test.45$state == state
    meth <-
      all_CpG_complete_with_test.45[index, as.character(meta$sample)]
    df_state$mean_delta[r] <-
      mean(all_CpG_complete_with_test.45$pearson.delta[index])
    df_state$sd_delta[r] <-
      sd(all_CpG_complete_with_test.45$pearson.delta[index])
    df_state$mean_minus_log_p[r] <-
      mean(-log(all_CpG_complete_with_test.45$pearson.p_val[index]), na.rm = TRUE)
    df_state$sd_minus_log_p[r] <-
      sd(-log(all_CpG_complete_with_test.45$pearson.p_val[index]), na.rm = TRUE)
    df_state$mean_methylation[r] <-
      mean(as.matrix(meth))
    df_state$sd_methylation[r] <-
      sd(as.matrix(meth))
    df_state$n_total[r] <- sum(index)
    df_state$n_H3K27AC[r] <-
      sum(index & all_CpG_complete_with_test.45$name != "NO_CHIP")
    r <- r + 1
  }
  
  n_all_CpG <- sum(df_state$n_total)
  n_H3K27AC_CpG <- sum(df_state$n_H3K27AC)
  expected_fraction <- n_H3K27AC_CpG / n_all_CpG
  
  df_state$expected <-
    df_state$n_total * expected_fraction
  
  df_state$state <-
    factor(x = df_state$state,
           levels = state_names)
  
  
  matching_indices <-
    match(
      df_state$state,
      paste(
        df_universal_annotation$state_order_by_group,
        df_universal_annotation$`States  presented in paper`,
        sep = "_"
      )
    )
  df_state$state <-
    factor(df_state$state, levels = df_state$state)
  df_state$state_color <-
    factor(df_universal_annotation$color[matching_indices],
           levels = unique(df_universal_annotation$color))
  df_state$state_group <-
    factor(df_universal_annotation$Group[matching_indices],
           levels = unique(df_universal_annotation$Group))
  
  df_state$p_val_test_all_methylation_KW <- NA
  
  for (st in state_names) {
    print(st)
    meth <-
      data.table::setDT(all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$state == st, as.character(meta$sample)])
    meth <-
      data.table::melt(
        data = meth,
        id.vars = NULL,
        variable.name = "sample",
        value.name = "methylation"
      )
    # # Perform Kruskal-Wallis test
    kw_state <- kruskal.test(methylation ~ sample, data = meth)
    df_state$p_val_test_all_methylation_KW[df_state$state == st] <-
      format.pval(kw_state$p.value, digits = 3)
  }
  
  saveRDS(
    object = df_state,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("df_state", "rds", sep = ".")
    )
  )
}
# done df_state_group  -----------------------------------------------------------------------
if (compute_df_state_group) {
  
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("all_CpG_complete_with_test.45", "rds", sep = ".")
    )
  )
  
  state_group_column_names <-
    c(
      "state_group",
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
  df_state_group <-
    data.frame(matrix(
      NA,
      nrow = length(unique(df_universal_annotation$Group)),
      ncol = length(state_group_column_names)
    ))
  
  # Rename the columns
  colnames(df_state_group) <- state_group_column_names
  
  r <- 1
  for (sg in  unique(df_universal_annotation$Group)) {
    print(sg)
    
    df_state_group$state_group[r] <- sg

    index <- all_CpG_complete_with_test.45$state_group == sg
    meth <-
      all_CpG_complete_with_test.45[index, as.character(meta$sample)]
    
    df_state_group$mean_delta[r] <-
      mean(all_CpG_complete_with_test.45$pearson.delta[index])
    df_state_group$sd_delta[r] <-
      sd(all_CpG_complete_with_test.45$pearson.delta[index])
    df_state_group$mean_minus_log_p[r] <-
      mean(-log(all_CpG_complete_with_test.45$pearson.p_val[index]), na.rm = TRUE)
    df_state_group$sd_minus_log_p[r] <-
      sd(-log(all_CpG_complete_with_test.45$pearson.p_val[index]), na.rm = TRUE)
    df_state_group$mean_methylation[r] <-
      mean(as.matrix(meth))
    df_state_group$sd_methylation[r] <-
      sd(as.matrix(meth))
    df_state_group$n_total[r] <- sum(index)
    df_state_group$n_H3K27AC[r] <-
      sum(index & all_CpG_complete_with_test.45$name != "NO_CHIP")
    r <- r + 1
  }

  n_all_CpG <- sum(df_state_group$n_total)
  n_H3K27AC_CpG <- sum(df_state_group$n_H3K27AC)
  expected_fraction <- n_H3K27AC_CpG / n_all_CpG
  
  df_state_group$expected_n_H3K27AC <-
    df_state_group$n_total * expected_fraction
  
  df_state_group$state_group <- factor(x = df_state_group$state_group,levels = df_state_group$state_group)
  
  df_state_group$color <- unique(as.character(df_universal_annotation$color))
  
 saveRDS(
    object = df_state_group,
    file = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("df_state_group", "rds", sep = ".")
    )
  )
}
# done make_plots  -----------------------------------------------------------------------
if (compute_make_plots) {
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("all_CpG_complete_with_test.45", "rds", sep = ".")
    )
  )
  
  load_variable_if_not_exists(
    variable_name = "df_state_sample",
    file_path = "12.pipeline/results/WholeGenome/df_state_sample.rds"
  )
  
  load_variable_if_not_exists(
    variable_name = "df_state",
    file_path = "12.pipeline/results/WholeGenome/df_state.rds"
  )
  
  load_variable_if_not_exists(
    variable_name = "df_state_group",
    file_path = "12.pipeline/results/WholeGenome/df_state_group.rds"
  )
  
  load_variable_if_not_exists(
    variable_name = "df_state_group_sample",
    file_path = "12.pipeline/results/WholeGenome/df_state_group_sample.rds"
  )
  
  # # Convert your dataframe to a data.table
  # # Melt the dataframe
  # all_CpG_meth <- data.table::melt(data = meth, id.vars = NULL, variable.name = "sample", value.name = "methylation")

  
  library(grid)
# create rectanlges 
    muh_grob <- grid::rectGrob(
    x=1:(length(unique(df_universal_annotation$Group))-1), y=0, gp=gpar(
      color='black', fill=unique(df_universal_annotation$color)[-1], alpha=1))
  

  # Figure 1B    
    p <- ggplot(data = df_state_group_sample[! as.character(df_state_group_sample$state_group) == "others", ],
                mapping = aes(x = state_group,y = mean_methylation,color = type))+
      geom_point(position = position_jitter(width = 0.2, height = 0), size = 1, alpha = 0.8)+
    theme_classic(base_size = 14)+
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          face = "bold"
        ))+
      coord_cartesian(clip='off')+
      annotation_custom(
        grob=muh_grob, xmin = 0, xmax = 1, ymin = -0.16, ymax=0.15
      )+
      labs(x=NULL,y="mean sample methylation")+
      scale_color_discrete(name = "lifestyle/diet")+
      theme(
        text = element_text(size = 14),  # Adjust size as needed
        axis.text = element_text(face = "bold"),  # Make all axis text bold
        axis.title = element_text(face = "bold"),  # Make axis titles bold
        legend.text = element_text(face = "bold")  # Make legend text bold
      )+
      scale_x_discrete(labels = function(x) ifelse(x == "quescient", "quiescent", ifelse(x == "HET", "heterochromatin", x)))
  
  
  ggsave(filename = file.path(
    OUTPUT_FOLDER,
    "results",
    "WholeGenome",
    paste("mean_methylation_per_in_groupV2", "png", sep = ".")
  ),plot = p,width = 9,height = 8.7)
  
 if(FALSE){
   # overlay a modern and a new sample
  # Define breaks
  hist_breaks <- seq(0, 1, 0.01)
  
  # Plot first histogram (S1148)
  hist(all_CpG_complete_with_test.45$S1148[all_CpG_complete_with_test.45$state_group == "TSS"],
       breaks = hist_breaks,
       col = rgb(1, 0, 0, 0.5),  # Red with 50% transparency
       xlim = c(0, 1),
       ylim = c(0, 15000),  # Adjust as needed
       main = "Methylation Distributions at TSS",
       xlab = "Methylation Level",
       ylab = "Frequency"
  )
  
  # Add second histogram (I2978) on top
  hist(all_CpG_complete_with_test.45$I2978[all_CpG_complete_with_test.45$state_group == "TSS"],
       breaks = hist_breaks,
       col = rgb(0, 0, 1, 0.5),  # Blue with 50% transparency
       add = TRUE
  )
  
  # Add legend
  legend("topright",
         legend = c("S1148", "I2978"),
         fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)),
         border = NA)
  
}
  
  
  df_plot_no_arf <-
    all_CpG_complete_with_test.45[!grepl(pattern = "GapArtf",
                                      x = all_CpG_complete_with_test.45$state,
                                      ignore.case = TRUE), ]
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
    scale_fill_manual(values = unique(df_universal_annotation$color)[-1])+
    theme(
      text = element_text(size = 14),  # Adjust size as needed
      axis.text = element_text(face = "bold"),  # Make all axis text bold
      axis.title = element_text(face = "bold"),  # Make axis titles bold
      legend.text = element_text(face = "bold")  # Make legend text bold
    )+
    scale_x_discrete(labels = function(x) ifelse(x == "quescient", "quiescent", ifelse(x == "HET", "heterochromatin", x)))
  
  
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
    scale_fill_manual(values = unique(df_universal_annotation$color)[-1])+
    theme(
      text = element_text(size = 14),  # Adjust size as needed
      axis.text = element_text(face = "bold"),  # Make all axis text bold
      axis.title = element_text(face = "bold"),  # Make axis titles bold
      legend.text = element_text(face = "bold")  # Make legend text bold
    )+
    scale_x_discrete(labels = function(x) ifelse(x == "quescient", "quiescent", ifelse(x == "HET", "heterochromatin", x)))
  
  
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
    scale_fill_manual(values = unique(df_universal_annotation$color)[-1])+
    theme(
      text = element_text(size = 14),  # Adjust size as needed
      axis.text = element_text(face = "bold"),  # Make all axis text bold
      axis.title = element_text(face = "bold"),  # Make axis titles bold
      legend.text = element_text(face = "bold")  # Make legend text bold
    )+
    scale_x_discrete(labels = function(x) ifelse(x == "quescient", "quiescent", ifelse(x == "HET", "heterochromatin", x)))
  
  
  
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
    scale_fill_manual(values = unique(df_universal_annotation$color)[-1])+
    theme(
      text = element_text(size = 14),  # Adjust size as needed
      axis.text = element_text(face = "bold"),  # Make all axis text bold
      axis.title = element_text(face = "bold"),  # Make axis titles bold
      legend.text = element_text(face = "bold")  # Make legend text bold
    )+
    scale_x_discrete(labels = function(x) ifelse(x == "quescient", "quiescent", ifelse(x == "HET", "heterochromatin", x)))
  
  
  # ggsave(filename = file.path(
  #   OUTPUT_FOLDER,
  #   "results",
  #   "WholeGenome",
  #   paste("methylation_WG_complete_with_test", "png", sep = "."))
  #   ,plot = p,width = 36,height = 12)
  
  index <- !grepl(pattern = "GapArtf",
                  x = df_state_sample$state,
                  ignore.case = TRUE)
    df_plot <-
    df_state_sample[index, ]
  
  muh_grob <- grid::rectGrob(
    x=1:length(unique(df_plot$state)), y=0, gp=gpar(
      color='black', fill=as.character(df_state$state_color[df_state$state %in% df_plot$state]), alpha=1))
  
  fill_vec <-
    df_universal_annotation$color[(100 - length(unique(df_plot$state))):length(unique(df_plot$state))]
  
  p <-
    ggplot(data = df_plot ,
           mapping = aes(x = state, y = mean_methylation)) +
    #geom_boxplot(outlier.shape = NA)+
    geom_point(position = position_jitter(width = 0.2, height = 0),size = 1, alpha = 0.7, mapping = aes(color = type)) +
    ggplot2::theme(
      axis.text.x = element_blank(),
      # axis.text.x = element_text(
      #   angle = 90,
      #   vjust = 0.5,
      #   hjust = 1,
      #   face = "bold"
      # ),
      axis.text.y = element_text(face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      panel.background = ggplot2::element_blank(),
      axis.ticks.x = element_blank()
    ) +
    ylab("mean sample methylation") +
    #xlab("chromatin state labeling")+
    annotation_custom(
      grob=muh_grob, xmin = 0, xmax = 1, ymin = -0.05, ymax=0.1
    )+xlab("")+
    coord_cartesian(clip='off')
  

  # Assuming df_plot and muh_grob are already defined
  p <- ggplot(data = df_plot, aes(x = state, y = mean_methylation)) +
    geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.8, alpha = 0.7, aes(color = type)) +
    theme(
      axis.text.x = element_blank(),  # Hide default x-axis text
      axis.text.y = element_text(face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      panel.background = element_blank(),
      axis.ticks.x = element_blank()  # Hide default x-axis ticks
    ) +
    labs(y = "mean sample methylation", x = "") +
    annotation_custom(grob = muh_grob, xmin = 0, xmax = 1, ymin = -0.2, ymax = 0.2) +
    coord_cartesian(clip = 'off') +
    geom_text(aes(label = state), y = -0.001, vjust = 0.5,hjust = 1, size = 2,angle = 90)+  # Manually add x-axis labels
    theme(legend.position="none")
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("100methylation_WG_complete_with_test", "png", sep = ".")
    )
    ,
    plot = p,
    width = 17,
    height = 7
  )
  
  
  
  df_plot_summery_by_state <-
    df_state[!grepl(pattern = "GapArtf", x = df_state$state),]
  
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
  
  p <-
    ggplot(data = df_state,
           mapping = aes(
             x = state,
             y = n_H3K27AC / expected,
             fill = state_group
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
    plot = p,
    width = 36,
    height = 12
  )
  
  df_state$p_value_NO_CHIP <- NA
  df_state$p_value_CHIP <- NA
  
  for (s in 1:length(state_names)) {
    print(s)
    data <-
      df_state_sample[state_names[s] == df_state_sample$state, ]
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
    
    df_state$p_value_NO_CHIP[as.character(df_state$state) == state_names[s]] <-
      p_value_NO_CHIP
    
    # compare H3K37Ac vs no chip
    n <-
      df_state$n_H3K27AC[as.character(df_state$state) == state_names[s]]
    
    if (n > 3) {
      # Perform Pearson correlation test
      cor_test_CHIP <-
        cor.test(data$age, data$methylation_H3K27AC, method = "pearson")
      # Extract the p-value from the test results
      p_value_CHIP <- cor_test_CHIP$p.value
      df_state$p_value_CHIP[as.character(df_state$state) == state_names[s]] <-
        p_value_CHIP
    }
    
    
    
    n1 <-
      df_state$n_total[as.character(df_state$state) == state_names[s]] - df_state$n_H3K27AC[as.character(df_state$state) == state_names[s]]
    
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
  
  
  big_labels <- as.character(df_state$state)
  big_labels[!as.character(df_state$state) %in% head(as.character(df_state$state)[order(df_state$n_H3K27AC, decreasing = TRUE)], 10)] <-
    " "
  df_state$big_labels <- big_labels
  
  # Create the pie chart with borders for all 100 states
  p <-
    ggplot(data = df_state, aes(
      x = 1,
      y = n_H3K27AC,
      fill = state,
      label = big_labels
    )) +
    coord_polar(theta = "y", start = 0) +  # Convert the bar plot into a pie chart
    labs(fill = "State", x = NULL) +   # Legend title
    theme_void() +  # Remove unnecessary elements
    geomtextpath::geom_textpath(
      data = df_state,
      position = position_stack(vjust = 0.5),
      aes(x = 1.75, ),
      size = 3.4,
      show.legend = FALSE
    ) +
    geom_bar(
      width = 1,
      stat = "identity",
      color = "black",
      # Add borders to segments
      #mapping = aes(fill = state_group)
      ) +
      #theme(legend.position = "none")+
      scale_fill_manual(values = as.character(df_state$state_color))

  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      "plots",
      paste("100pie_chart H3K27Ac", "png", sep = ".")
    )
    ,
    plot = p,
    width = 12,
    height = 12
  )
  
  df_state_group$state_group <- factor(x = df_state_group$state_group,levels = unique(df_state_group$state_group))
  big_labels <- as.character(df_state_group$state_group)
  big_labels[!big_labels %in% head(big_labels[order(df_state_group$n_H3K27AC, decreasing = TRUE)], 3)] <- " "
  df_state_group$big_labels <- big_labels
  
  
  p <-
    ggplot(data = df_state_group, aes(
      x = 1,
      y = n_H3K27AC,
      fill = state_group,
      label = big_labels
    )) +
    coord_polar(theta = "y", start = 0) +  # Convert the bar plot into a pie chart
    labs(fill = "state_group", x = NULL) +   # Legend title
    theme_void() +  # Remove unnecessary elements
    geom_bar(
      width = 1,
      stat = "identity",
      color = "black",
      # Add borders to segments
      #mapping = aes(fill = state_group)
    ) +
    geomtextpath::geom_textpath(
      data = df_state_group,
      position = position_stack(vjust = 0.5),
      aes(x = 1.1 ),
      size = 9,
      show.legend = FALSE
    )+
    #theme(legend.position = "none")+
    scale_fill_manual(values = as.character(unique(df_universal_annotation$color)))
  
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      "plots",
      paste("pie_chart_state_group_H3K27Ac", "png", sep = ".")
    )
    ,
    plot = p,
    width = 12,
    height = 12
  )
  
  # plot for each sample mean methylation of chromatin state
             ggplot(data = df_state_sample,
                    mapping = aes(x = sample, y = mean_methylation, color = type)) +
               #geom_boxplot()+
               geom_point(mapping = aes(color = state_group))
             
  
             
  #ggplot(data = df_state_group,mapping = aes(x = ))           
             
}

# done state_dependen_analysis --------------------------------------------------------
if (state_dependen_analysis) {
  # load data
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("all_CpG_complete_with_test.45", "rds", sep = ".")
    )
  )
  
  
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

  for (s in 1:length(state_names)) {
    start_time <- Sys.time()
    state <- state_names[s]
    print(state)
    df_CpG_in_state <-
      all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$state == state,]
    
    df_CpG_in_state$pearson.p_val[df_CpG_in_state$pearson.delta == 0] <- 1
    print(paste("# CpG in :", state, nrow(df_CpG_in_state)))
    saveRDS(object = df_CpG_in_state, file = file.path(folder_path, paste(state, "rds", sep = ".")))
    
    pearson.state <-
      list(
        sim_typisation = real_age_typisation,
        permuation_order = state,
        data = df_CpG_in_state[, c("pearson.p_val", "pearson.statistic", "pearson.delta")]
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
  
  #summery_all_state_results <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/results/WholeGenome/analysis_by_state/summery.pearson.all.states.rds")
  
  summery_all_state_results_wide <-
    pivot_wider(data = summery_all_state_results,
                names_from = slope,
                values_from = fraction_direction_signfincant)
  
  dummy_value <- 1
  summery_all_state_results_wide$negative_enrichment <-
    (summery_all_state_results_wide$negative_signfincant_CpG) / (dummy_value + summery_all_state_results_wide$positive_signfincant_CpG)
  
  alpha <- 8
  min_delta <- 0.5
  ggplot(
    data =  summery_all_state_results_wide[summery_all_state_results_wide$min_delta == min_delta &
                                             summery_all_state_results_wide$minus_log_alpha >= 6, ],
    mapping = aes(x = minus_log_alpha,
                  y = negative_enrichment,
                  color = permuation_type)
  ) +
    geom_point(position = "jitter") +
    theme(legend.position = "none")
  
  
  slice_0 <-
    summery_all_state_results_wide[summery_all_state_results_wide$min_delta == min_delta &
                                     summery_all_state_results_wide$minus_log_alpha == 0, ]
  
  slice_alpha <-
    summery_all_state_results_wide[summery_all_state_results_wide$min_delta == min_delta &
                                     summery_all_state_results_wide$minus_log_alpha == alpha, ]
  
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
      paste("fraction of significant CpG ", p[value] <= e ^ (-8), ", ", delta[meth] > 0.5)
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
  #all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$pearson.delta > min_delta]
  
  
  ggplot(
    data = summery_all_state_results[summery_all_state_results$min_delta == min_delta &
                                       summery_all_state_results$minus_log_alpha >= 6, ],
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
                                summery_all_state_results$minus_log_alpha == alpha,]
  slice[order(slice$fraction_direction_signfincant, decreasing = TRUE), c("permuation_type",
                                                                          "slope",
                                                                          "fraction_direction_signfincant")]
  
}
# done genome wide only delta good states --------------------------------------------------
if (FALSE) {
  min_delta <- 0.5
  max_q <- 0.01
  
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(OUTPUT_FOLDER,"results", "WholeGenome", "all_CpG_complete_with_test.45.rds")
  )
  
  #total <- sum(all_CpG_complete_with_test.45$pearson.delta > min_delta  & all_CpG_complete_with_test.45$pearson.p_val <= max_p)
  
  # genome_CpG_with_delta <-
  #   all_CpG_complete_with_test.45[!grepl(pattern = "GapArtf",
  #                                     x = all_CpG_complete_with_test.45$state,
  #                                     ignore.case = TRUE) &
  #                                all_CpG_complete_with_test.45$pearson.delta > min_delta,]
  
  
  #p.multiple <- p.adjust(all_CpG_complete_with_test.45$pearson.p_val, method = "BY")
  #sum(p-adjust < )
  
  genome_CpG_with_delta <-
    all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$pearson.delta > min_delta,]
  
  counts_min_delta_by_state_group <- genome_CpG_with_delta %>%
    count(state_group) 
  
  counts_min_delta <- genome_CpG_with_delta %>%
    count(state)
  
  significant_genome_CpG_with_delta <-
      genome_CpG_with_delta[genome_CpG_with_delta$pearson.q_vals.min_delta_0.5 <= max_q,]
  
  counts_sig_min_delta <- significant_genome_CpG_with_delta %>%
    count(state)
  
  
  #plot_age_correlation(df_row = significant_genome_CpG_with_delta[1,],typisation = real_age_typisation)
  
  significant_genome_CpG_with_delta$slope <-
    ifelse(significant_genome_CpG_with_delta$pearson.statistic <= 0,
           "negative",
           "positive")
  
  # Example calculation
  result <- significant_genome_CpG_with_delta %>%
    group_by(state) %>%
    summarize(
      fraction_negative = mean(slope == "negative", na.rm = TRUE),
      fraction_positive = mean(slope == "positive", na.rm = TRUE)
    )
  
  counts_sig_min_delta_slope <- significant_genome_CpG_with_delta %>%
    mutate(
      state = factor(state),  # Convert to factor if not already
      slope = factor(slope)   # Convert to factor ensuring all possible slopes are levels
    ) %>%
    group_by(state, slope) %>%
    count() %>%
    ungroup() %>%
    complete(state, slope, fill = list(n = 0))  # Fill missing combinations with zero counts
  
  counts_sig_min_delta_slope_by_state_group <- significant_genome_CpG_with_delta %>%
    group_by(state_group) %>%
    count(slope)
  
  counts_sig_min_delta_slope_by_state_group <- counts_sig_min_delta_slope_by_state_group %>%
    left_join(counts_min_delta_by_state_group, by = "state_group", suffix = c("", "_total"))
  
  counts_sig_min_delta_by_state_group <- significant_genome_CpG_with_delta %>%
    count(state_group)
  counts_sig_min_delta_by_state_group <- counts_sig_min_delta_by_state_group %>%
    left_join(counts_min_delta_by_state_group,by = "state_group", suffix = c("", "_total"))
  counts_sig_min_delta_by_state_group <- counts_sig_min_delta_by_state_group %>%
    mutate(proportion = n / n_total)
  counts_sig_min_delta_by_state_group$state_color <- unique(df_universal_annotation$color)
  
  
  # change HET to heterochromatin etc. 
  counts_sig_min_delta_by_state_group$state_group <- ifelse(
      counts_sig_min_delta_by_state_group$state_group  == "quescient",
      "quiescent",
      ifelse(
        counts_sig_min_delta_by_state_group$state_group  == "HET",
        "heterochromatin",
        as.character(counts_sig_min_delta_by_state_group$state_group)
      )
    )
  
  # # Prepare the data for the Chi-square test
  # observed <- with(counts_sig_min_delta_by_state_group, data.frame(
  #   significant = n,
  #   not_significant = n_total - n
  # ))
  # 
  # # Conduct the Chi-square test
  # chi_test_result <- chisq.test(observed)
  # 
  # # Print the results
  # print(chi_test_result)
  # 
  # 
  # # goodnwes of fit test
  # 
  # # Assuming counts_sig_min_delta_by_state_group is already defined with columns n and n_total
  # 
  # # Calculate total observed significant CpGs
  # total_observed <- sum(counts_sig_min_delta_by_state_group$n)
  # 
  # # Total number of CpGs across all groups
  # total_groups <- sum(counts_sig_min_delta_by_state_group$n_total)
  # 
  # # Expected uniform proportion
  # expected_uniform_proportion <- total_observed / total_groups
  # 
  # # Assuming a uniform distribution, calculate expected count for each group
  # counts_sig_min_delta_by_state_group$expected <- 
  #   (total_observed / total_groups) * counts_sig_min_delta_by_state_group$n_total
  # 
  # # Perform the Chi-square Goodness-of-Fit test for the entire dataset
  # observed_counts <- counts_sig_min_delta_by_state_group$n
  # expected_counts <- counts_sig_min_delta_by_state_group$expected
  # 
  # # Correcting the usage of chisq.test for goodness-of-fit rather than test of independence
  # chi_test_result <- chisq.test(x = observed_counts, p = expected_counts / sum(expected_counts), rescale.p = TRUE)
  # 
  # # Obtain standardized residuals
  # counts_sig_min_delta_by_state_group$stdres <- chi_test_result$stdres
  # 
  # # Convert standardized residuals to approximate p-values
  # counts_sig_min_delta_by_state_group$p_value_from_residuals <- 2 * pnorm(abs(counts_sig_min_delta_by_state_group$stdres), lower.tail = FALSE)
  # 
  # # Adjusting for multiple testing using Bonferroni correction
  # counts_sig_min_delta_by_state_group$p_value_adjusted <- p.adjust(counts_sig_min_delta_by_state_group$p_value_from_residuals, method = "bonferroni")
  # 
  # # Update significance column to include significance only for positive residuals
  # counts_sig_min_delta_by_state_group$significance <-
  #   ifelse(
  #     counts_sig_min_delta_by_state_group$stdres > 0,
  #     ifelse(
  #       counts_sig_min_delta_by_state_group$p_value_adjusted <= 0.001,
  #       '***',
  #       ifelse(
  #         counts_sig_min_delta_by_state_group$p_value_adjusted <= 0.01,
  #         '**',
  #         ifelse(
  #           counts_sig_min_delta_by_state_group$p_value_adjusted <= 0.05,
  #           '*',
  #           ''
  #         )
  #       )
  #     ),
  #     ''
  #   )
  
  counts_sig_min_delta_by_state_group$state_group <-
    factor(counts_sig_min_delta_by_state_group$state_group,
           levels = counts_sig_min_delta_by_state_group$state_group)
  #Figure 1c
  ggplot(data = counts_sig_min_delta_by_state_group, aes(x = state_group, y = proportion, fill = state_group)) +
    geom_col() + # Use geom_col for pre-calculated values
    #geom_text(aes(label = significance, y = proportion + 0.0005), position = position_dodge(width = 0.9), vjust = -0.5, check_overlap = TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate X labels for readability
    labs(
      x = "",
      y = expression("Proportion of significant CpGs (q-value <= max_q)" ~ "from CpGs with" ~ delta > min_delta),
      fill = "State Group"
    ) +
    scale_fill_manual(values = unique(df_universal_annotation$color))
  
  
  # ggplot(data = counts_sig_min_delta_by_state_group, aes(x = state_group, y = proportion, fill = state_group)) +
  #   geom_col() + # Use geom_col for pre-calculated values
  #   geom_text(aes(label = significance, y = proportion + 0.0005), position = position_dodge(width = 0.9), vjust = -0.5, check_overlap = TRUE) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate X labels for readability
  #   labs(
  #     x = "",
  #     y = expression("Proportion of significant CpGs (p <= 0.001)" ~ "from CpGs with" ~ delta > 0.5),
  #     fill = "State Group"
  #   ) +
  #   scale_fill_manual(values = unique(df_universal_annotation$color))
  
#Figure 1C
  
  # p <- ggplot(data = counts_sig_min_delta_by_state_group, aes(x = state_group, y = proportion, fill = state_group)) +
  #   geom_col(color = "black", show.legend = FALSE) +  # Use geom_col for pre-calculated values and hide legend
  #   geom_text(
  #     aes(label = significance, y = proportion + 0.001),
  #     position = position_dodge(width = 0.9),
  #     vjust = -0.5, size = 6,  # Increased size for better visibility
  #     fontface = "bold",  # Make in-plot text bold
  #     check_overlap = TRUE
  #   ) +
  #   geom_hline(yintercept = expected_uniform_proportion, linetype = "dashed", color = "blue", size = 1) +
  #   scale_fill_manual(values = unique(df_universal_annotation$color)) +  # Apply your color palette
  #   labs(
  #     x = NULL,  # Removing the x-axis label for clarity
  #     y = expression("Proportion of significant CpGs" ~ (p <= 0.001) ~ "from CpGs with" ~ delta[meth] > 0.5)
  #   ) +
  #   theme_minimal(base_size = 14) +  # Use a minimal theme with a larger base font size
  #   theme(
  #     axis.text = element_text(face = "bold"),  # Make all axis text bold
  #     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjusting x-axis labels for readability
  #     axis.title = element_text(face = "bold"),  # Make axis titles bold
  #     plot.title = element_text(face = "bold"),  # Make plot title bold
  #     plot.subtitle = element_text(face = "bold"),  # Make plot subtitle bold
  #     legend.title = element_text(face = "bold"),  # Make legend title bold
  #     legend.text = element_text(face = "bold"),  # Make legend text bold
  #     panel.grid.major = element_blank(),  # Clean background without major grid lines
  #     panel.grid.minor = element_blank(),  # Clean background without minor grid lines
  #     legend.position = "none"  # Removing the legend to reduce clutter
  #   ) # +
  #   # scale_x_discrete(labels = function(x) ifelse(x == "quescient", "quiescent", ifelse(x == "HET", "heterochromatin", x)))
  # 
  # ggsave(
  #   filename = file.path(
  #     OUTPUT_FOLDER,
  #     "results",
  #     "WholeGenome",
  #     "plots",
  #     paste("fraction significant CpG per state group", "png", sep = ".")
  #   )
  #   ,
  #   plot = p,
  #   width = 10,
  #   height = 8
  # )
  
  
  ## reviewer control figure : major comment
  
  # # Compute mean and variance of methylation for each sample
  # df_methylation_stats <- all_CpG_complete_with_test.45 %>%
  #   select(all_of(as.character(meta$sample))) %>%
  #   summarise(across(everything(), list(Mean_Methylation = ~mean(.x, na.rm = TRUE), 
  #                                       Variance_Methylation = ~var(.x, na.rm = TRUE)))) %>%
  #   pivot_longer(cols = everything(), names_to = c("sample", ".value"), names_sep = "_")
  # 
  # # Merge with metadata
  # df_meta_methylation <- merge(meta, df_methylation_stats, by = "sample")
  # 
  # # Scatter plot for mean methylation with variance as error bars
  # ggplot(df_meta_methylation, aes(x = age_mean_BP, y = Mean,color = df_meta_methylation$Type )) +
  #   geom_point(size = 3) +  # Mean methylation points
  #   geom_errorbar(aes(ymin = Mean - sqrt(Variance), 
  #                     ymax = Mean + sqrt(Variance)), 
  #                 width = 0, color = "gray50") +  # Error bars
  #   scale_x_reverse() +  # Reverse time axis so older samples are on the left
  #   theme_minimal() +
  #   labs(
  #     x = "Age (Years BP)",
  #     y = "Mean Methylation Level",
  #     title = "Mean Methylation Over Time with Variability",
  #     caption = "Error bars represent sqrt(variance) as a measure of spread"
  #   )
  
  
  # # Load required libraries
  # library(dplyr)
  # library(matrixStats)
  # 
  # # Convert state_group to character to ensure correct comparison
  # all_CpG_complete_with_test.45$state_group <- as.character(all_CpG_complete_with_test.45$state_group)
  # 
  # # Identify unique chromatin state groups
  # state_groups <- unique(all_CpG_complete_with_test.45$state_group)
  # 
  # # Initialize an empty long-format data frame to store results
  # df_chromatin_variance_long <- data.frame(state_group = character(), Variance = numeric(), stringsAsFactors = FALSE)
  # 
  # # Loop through each chromatin state group
  # for (state in state_groups) {
  #   print(paste("Processing Chromatin State Group:", state))  # Print current state
  #   
  #   # Filter data to include only rows corresponding to the current chromatin state
  #   df_subset <- all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$state_group == state, ] 
  #   
  #   # Select only numeric sample columns (biosamples)
  #   df_numeric <- df_subset %>% select(all_of(meta$sample))
  #   
  #   # Ensure conversion to numeric matrix
  #   df_numeric <- as.matrix(df_numeric)
  #   df_numeric <- apply(df_numeric, 2, as.numeric)  # Convert all columns to numeric
  #   
  #   # Compute row-wise variance for all CpG sites
  #   variance_values <- rowVars(df_numeric, na.rm = TRUE)
  #   
  #   # Create a long-format data frame with all variance values and chromatin state labels
  #   df_variance_state <- data.frame(state_group = rep(state, length(variance_values)), Variance = variance_values)
  #   
  #   # Append results to the long-format data frame
  #   df_chromatin_variance_long <- rbind(df_chromatin_variance_long, df_variance_state)
  # }
  # 
  # # Print the first rows of the final long-format variance table
  # head(df_chromatin_variance_long)
  # 
  # # Load required libraries
  # library(dplyr)
  # library(ggplot2)
  # 
  # # Compute mean variance and standard error (SEM) for each chromatin state group
  # df_variance_summary <- df_chromatin_variance_long %>%
  #   group_by(state_group) %>%
  #   summarise(
  #     Mean_Variance = mean(Variance, na.rm = TRUE),
  #     SEM = sd(Variance, na.rm = TRUE) / sqrt(n())
  #   ) %>%
  #   arrange(desc(Mean_Variance))
  # 
  # # Print summary data
  # print(df_variance_summary)
  # 
  # # Load required libraries
  # library(dplyr)
  # library(ggplot2)
  # 
  # # Ensure state_group and Group are of the same type
  # df_universal_annotation$Group <- as.character(df_universal_annotation$Group)
  # df_variance_summary$state_group <- as.character(df_variance_summary$state_group)
  # 
  # # Merge color information based on state_group
  # df_variance_summary <- df_variance_summary %>%
  #   left_join(df_universal_annotation, by = c("state_group" = "Group"))
  # 
  # # Print to verify merging
  # head(df_variance_summary)
  # 
  # 
  # # Create a bar plot with error bars (SEM)
  # ggplot(df_variance_summary, aes(x = reorder(state_group, Mean_Variance), y = Mean_Variance, fill = state_group)) +
  #   geom_bar(stat = "identity", show.legend = FALSE) +
  #   geom_errorbar(aes(ymin = Mean_Variance - SEM, ymax = Mean_Variance + SEM), width = 0.3, color = "black") +
  #   coord_flip() +  # Flip axes for better readability
  #   theme_minimal() +
  #   labs(
  #     x = "Chromatin State Group",
  #     y = "Mean Variance of CpG Methylation",
  #     title = "Mean Variance of CpG Methylation Across Chromatin State Groups",
  #     caption = "Error bars represent SEM (Standard Error of the Mean)"
  #   ) +
  #  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12))
  
  # Define the output directory
  subfolder_path <- file.path(this.dir, "review")
  if (!dir.exists(subfolder_path)) {
    dir.create(subfolder_path)
  }
  
  # Get the levels of the factor (faster than unique() for factors)
  state_groups <- levels(all_CpG_complete_with_test.45$state_group)
  
  # Write subsets efficiently
  for (group_name in state_groups) {
    cat("Processing:", group_name, "\n")  # Print progress
    
    # Create file path
    file_path <- file.path(subfolder_path, paste0(group_name, ".csv"))
    
    # Use `which()` to speed up subsetting (faster for large data frames)
    subset_idx <- which(all_CpG_complete_with_test.45$state_group == group_name)
    
    # Write directly to CSV to avoid creating large intermediate objects
    write.csv(all_CpG_complete_with_test.45[subset_idx, , drop = FALSE], file_path, row.names = FALSE)
  }
  
  # Completion message
  cat("Data has been split and saved in:", subfolder_path, "\n")
  
  # File path for a specific state group
  file_path <- file.path(subfolder_path, paste0(group_name, ".csv"))
  
  # Read the file into a data frame
  df_group <- read.csv(file_path)
  
  # Display the first few rows
  head(df_group)
  
  
  ## end reviewer insert 
  
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
  
  matching_indices <-
    match(
      counts_sig_min_delta_slope_by_state_group$state_group,
      unique(df_universal_annotation$Group)
    )
  counts_sig_min_delta_slope_by_state_group$state_color <- unique(df_universal_annotation$color)[matching_indices]
  # sanity check
  # sum(significant_genome_CpG_with_delta$slope[significant_genome_CpG_with_delta$state == "99_TSS2"] == "positive",na.rm = TRUE)
  # sum(significant_genome_CpG_with_delta$slope[significant_genome_CpG_with_delta$state == "99_TSS2"] == "negative",na.rm = TRUE)
  
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
    paste("fraction of significant CpG ", p[value] <= 0.001, ", ", delta[meth] > 0.5)
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
  
  summary_df_by_state_group <- counts_sig_min_delta_slope_by_state_group %>%
    group_by(state_group) %>%
    summarize(total = sum(n))
  
  
  summary_df_wide <- summary_df %>%
    pivot_wider(
      names_from = slope,
      values_from = c(n),
      names_glue = "{slope}_{.value}"
    )
  
  summary_df_wide_by_state_group <- counts_sig_min_delta_slope_by_state_group %>%
    pivot_wider(
      names_from = slope,
      values_from = c(n),
      names_glue = "{slope}_{.value}"
    )
  
  summary_df_wide$log_fraction <-
    log(summary_df_wide$negative_n / summary_df_wide$positive_n)
  
  
  summary_df_wide_by_state_group$log_fraction <-
    log(summary_df_wide_by_state_group$negative_n / summary_df_wide_by_state_group$positive_n)

  # p <- ggplot(
  #   data = summary_df_wide,#summary_df_wide_by_state_group,
  #   mapping = aes(
  #     x = state_group,
  #     y = log_fraction,
  #     fill = state_group,
  #     label = paste("Negative:", negative_n, " ; Positive:", positive_n)
  #   )
  # ) +
  #   geom_bar(stat = "identity", colour = "black") +
  #   geom_text(
  #     position = position_stack(vjust = 0.5),
  #     size = 7,
  #     mapping = aes(y  =  3)
  #   ) +  # Add text annotations
  #   ggplot2::theme(
  #     panel.background = ggplot2::element_blank(),
  #     # axis.text.x = element_text(
  #     #   angle = 90,
  #     #   vjust = 0.5,
  #     #   hjust = 1
  #     #)
  #     ) +
  #     ylab(expression(
  #       paste(log(negative[correlation] / positive[correlation]), " ", p[value] <= 0.001, ", ", delta[meth] > 0.5)
  #     )) +
  #       xlab("Chromatin state groups") +
  #       scale_fill_manual(values = unique(df_universal_annotation$color)[-1]) +
  #       coord_flip()+
  #   theme(
  #     text = element_text(size = 14)  # Adjust size as needed
  #   )
   
  # Figure S4
  
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
    theme(
      panel.background = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.line = element_line(color = "black", size = 1),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_blank(),  # If you want to show the y-axis title in bold, remove this line
      text = element_text(face = "bold")  # This makes all text bold by default, including legend text
    ) +
    ylab(expression(paste(
      log(negative[correlation] / positive[correlation]), " ", q[value] <= 0.01, ", ", delta[meth] > 0.5
    ))) +
    scale_fill_manual(values = unique(df_universal_annotation$color)) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0))
  
  ggsave(
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "WholeGenome",
      paste("enrichment significant CpG per Group", "png", sep = ".")
    ),
    plot = p,
    width = 12,
    height = 12
  )
  
      
      
      summary_df_state <-
        counts_sig_min_delta_slope %>%
        group_by(state, state_color, state_group) %>%
        summarize(
          negative_count = sum(ifelse(slope == "negative", n, 0)),
          positive_count = sum(ifelse(slope == "positive", n, 0)),
          log_fraction = log(negative_count / positive_count)
        )
      
      #Figure S5
      p <- ggplot(
        data = summary_df_state,
        mapping = aes(
          x = state,
          y = log_fraction,
          fill = state_group,
          #label = paste("Negative:", negative_count, " ; Positive:", positive_count)
        )
      ) +
        geom_bar(stat = "identity", colour = "black") +
        # geom_text(
        #   position = position_stack(vjust = 0.5),
        #   size = 3,
        #   mapping = aes(y  =  max(summary_df_state$log_fraction))
        # ) +  # Add text annotations
          theme(
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 12),
            axis.text.y = element_text(face = "bold", size = 12),
            axis.title.x = element_text(face = "bold", size = 14),
            axis.title.y = element_text(face = "bold", size = 14),
            legend.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"),
            text = element_text(face = "bold")
          )+
          ylab(expression(
            paste(log(negative[correlation] / positive[correlation]), " ", p[value] <= 0.001, ", ", delta[meth] > 0.5)
          )) +
            xlab("") +
        scale_fill_manual(values = unique(as.character(summary_df_state$state_color)))#     +
        # scale_x_discrete(labels = gsub(pattern = "_.*",replacement = "",x = summary_df_state$state))
          
          ggsave(
            filename = file.path(
              OUTPUT_FOLDER,
              "results",
              "WholeGenome",
              paste("enrichment significant CpG per state", "png", sep = ".")
            )
            ,
            plot = p,
            width = 18,
            height = 9
          )
          
          openxlsx::write.xlsx(
            x = summary_df_state,
            file.path(
              OUTPUT_FOLDER,
              "results",
              "genome_wide_tests.xlsx"
            ),
            sheetName = paste("pearson <= " ,max_p, " delta >",min_delta ),
            colNames = TRUE
          )
     
  
          
               
}
# done Calculate state specific CpG permuation simulations dependent on delta --------------------------------------------------
if (FALSE) {
  min_delta <- 0.5
  #max_p <- 0.001

  
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(this.dir,"12.pipeline/results/WholeGenome/all_CpG_complete_with_test.45.rds")
  )
  
  # drop kw.p_val kw.statistic  kw.delta, pearson.p_val_BH
  all_CpG_complete_with_test.45 <- all_CpG_complete_with_test.45[, !colnames(all_CpG_complete_with_test.45) %in% c("kw.p_val", "kw.statistic", "kw.delta", "pearson.p_val_BH")]
  
  
  OUTPUT_FOLDER_pearson<- file.path(OUTPUT_FOLDER, paste("45.pearson","CpG_permutation","min_delta",min_delta,sep = "."))
  #OUTPUT_FOLDER_kw<- file.path(OUTPUT_FOLDER, paste("45.KW","CpG_permutation","min_delta",min_delta,sep = "."))
  print("create_CpG_permutations_vertical")
  dir.create(path = OUTPUT_FOLDER_pearson,
             showWarnings = FALSE)
  
  # dir.create(path = OUTPUT_FOLDER_kw,
  #            showWarnings = FALSE)
  
  target_CpG_permutations <- 10000000  ### CHANGED: Adjust as needed
  
    for(chromatin_group in levels(all_CpG_complete_with_test.45$state_group)){
      
      
      print(paste("Processing", chromatin_group))
      # chromatin_group <- "TSS"
      df_temp <- all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$state_group == chromatin_group & 
                                                 all_CpG_complete_with_test.45$pearson.delta > min_delta,]

    
      num_CpGs_in_state <- nrow(df_temp) 
      
      ### CHANGED: keep generating permutations until target CpG permutations reached
      while (TRUE) {
        # Count how many permutations already exist
      
      ###  Check how many permutation files already exist for this group
      existing_files <- list.files(
        path = OUTPUT_FOLDER_pearson,
        pattern = paste0("^", chromatin_group, "\\.pearson\\.CpG_permutation\\..*\\.rds$"),
        full.names = TRUE
      )  
      
      num_existing_files <- length(existing_files)
      total_CpG_permutations <- num_existing_files * num_CpGs_in_state
      
      if (total_CpG_permutations >= target_CpG_permutations) {
        message(paste("Reached", total_CpG_permutations, "permuted CpGs for", chromatin_group, "- skipping."))
        break  #Exit the while-loop
      }
      start_time <- Sys.time()
      time_string <- format(start_time, "%Y_%m_%d_%H_%M_%S")
      
      print(paste("Permutation", num_existing_files + 1, "for", chromatin_group))
    permutation_CpG <-
      create_CpG_permution(df_temp, real_age_typisation)
    
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
    print("parallel_testing_pearson_cor complete")

    # save age permutation
    sim_age_file_name <-
      paste(chromatin_group,
            "pearson",
            "CpG_permutation",
            time_string,
            "rds",
            sep = ".")
    saveRDS(
      object = permutation_age,
      file = file.path(OUTPUT_FOLDER_pearson, sim_age_file_name)
    )
    
    # permutation_food$data <-
    #   parallel_testing_kruskall_valis(df = df_permuted, food.typisation = real_food_typisation)
    # print("parallel_testing_kruskall_valis")
    # # save food permutation
    # sim_food_file_name <-
    #   paste(chromatin_group,
    #         "KW",
    #         "CpG_permutation",
    #         time_string,
    #         "rds",
    #         sep = ".")
    # saveRDS(
    #   object = permutation_food,
    #   file = file.path(OUTPUT_FOLDER_kw, sim_food_file_name)
    # )
    
    end_time <- Sys.time()
    print(paste("Completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
    
    }
  }
  
  OUTPUT_FOLDER_pearson_plot <- file.path( OUTPUT_FOLDER_pearson,"plot")
  dir.create(path = OUTPUT_FOLDER_pearson_plot,
             showWarnings = FALSE)
  print(OUTPUT_FOLDER_pearson_plot)
  for(chromatin_group in levels(all_CpG_complete_with_test.45$state_group)){
  print(chromatin_group)
  true.pval <- all_CpG_complete_with_test.45$pearson.p_val[
    all_CpG_complete_with_test.45$state_group == chromatin_group & 
      all_CpG_complete_with_test.45$pearson.delta > min_delta]
  
  # List of matching RDS files
  state_group_permutation_age_files <- list.files(
    path = OUTPUT_FOLDER_pearson,
    pattern = paste0("^", chromatin_group, "\\.pearson\\.CpG_permutation\\."),
    full.names = TRUE
  )
  
  # Load and combine all p-values into one data.table
  all_pvals <- rbindlist(lapply(state_group_permutation_age_files, function(f) {
    df <- readRDS(f)$data
    data.table(pval = df$pearson.p_val, file = basename(f))
  }))
  
  # 1. Add true p-values as their own "file"
  true_pvals_dt <- data.table(
    pval = true.pval,
    file = "true_data"
  )
  
  # 2. Combine with permutation p-values
  all_pvals_combined <- rbind(all_pvals, true_pvals_dt)
  
  # 3. Plot all p-values including true data
  p <- ggplot(all_pvals_combined, aes(x = pval, fill = file)) +
    geom_histogram(binwidth = 0.01, alpha = 0.5, position = "identity", color = NA) +
    theme_minimal() +
    labs(
      title = paste("Permutation vs. True P-Values -", chromatin_group),
      x = "Pearson p-value",
      y = "Count",
      fill = "Source"
    ) +
    xlim(0, 1)
  
  ggsave(
    plot = p,
    filename = file.path(OUTPUT_FOLDER_pearson_plot,
      paste0(chromatin_group,".permutation_vs_true_data_p_value_destribution.png")),
    width = 10,
    height =10
  )
  
  library(data.table)
  library(ggplot2)
  library(patchwork)  ### CHANGED

  # positive vs negative

  # 1. TRUE DATA
  true_filter <- all_CpG_complete_with_test.45$state_group == chromatin_group &
    all_CpG_complete_with_test.45$pearson.delta > min_delta

  true_pvals_dt <- data.table(
    pval = all_CpG_complete_with_test.45$pearson.p_val[true_filter],
    direction = ifelse(all_CpG_complete_with_test.45$pearson.statistic[true_filter] >= 0, "positive", "negative"),
    file = "true_data"
  )

  # 2. PERMUTATION FILES
  state_group_permutation_age_files <- list.files(
    path = OUTPUT_FOLDER_pearson,
    pattern = paste0("^", chromatin_group, "\\.pearson\\.CpG_permutation\\."),
    full.names = TRUE
  )

  # 3. Load permutation p-values + direction
  all_pvals <- rbindlist(lapply(state_group_permutation_age_files, function(f) {
    df <- readRDS(f)$data
    data.table(
      pval = df$pearson.p_val,
      direction = ifelse(df$pearson.statistic >= 0, "positive", "negative"),
      file = basename(f)
    )
  }))

  true_pvals_dt$data_type <- "true_data"          # ✅ ADD this line
  all_pvals$data_type <- "permutation_data"       # ✅ ADD this line
  
  
  # 4. Combine TRUE and PERMUTATION data
  all_pvals_combined <- rbind(all_pvals, true_pvals_dt)

  # --- Negative and Positive p-value cutoffs
  pval_cutoff_neg <- max(
    all_CpG_complete_with_test.45$pearson.p_val[
      true_filter &
        all_CpG_complete_with_test.45$pearson.q_vals.min_delta_0.5.negative < max_q
    ],
    na.rm = TRUE
  )

  pval_cutoff_pos <- max(
    all_CpG_complete_with_test.45$pearson.p_val[
      true_filter &
        all_CpG_complete_with_test.45$pearson.q_vals.min_delta_0.5.positive < max_q
    ],
    na.rm = TRUE
  )

  cat("P-value cutoff for negative:", pval_cutoff_neg, "\n")
  cat("P-value cutoff for positive:", pval_cutoff_pos, "\n")

  cutoffs_df <- data.frame(
    direction = c("negative", "positive"),
    cutoff = c(pval_cutoff_neg, pval_cutoff_pos)
  )

  # --- 5. Plot 1 (linear scale)
  p_linear <- ggplot(all_pvals_combined, aes(x = pval, fill = file)) +
    geom_histogram(binwidth = 0.01, alpha = 0.5, position = "identity", color = NA) +
    facet_wrap(~direction) +
    geom_vline(
      data = cutoffs_df,
      aes(xintercept = cutoff),
      color = "red",
      linetype = "dashed",
      size = 0.8
    ) +
    theme_minimal() +
    labs(
      title = paste("Permutation vs. True P-Values (Linear Scale) -", chromatin_group),
      x = "Pearson p-value",
      y = "Count",
      fill = "Source"
    ) +
    xlim(0, 1)+
    scale_fill_manual(
      values = c("true_data" = "blue", "permutation_data" = "gray50"),   # ✅ ADD
      name = "Data Source"                                               # ✅ ADD
    )
  

  # --- 6. Plot 2 (log scale)  ### CHANGED
  p_log <- ggplot(all_pvals_combined, aes(x = pval, fill = file)) +
    geom_histogram(binwidth = 0.5, alpha = 0.5, position = "identity", color = NA) +  ### CHANGED
    facet_wrap(~direction) +
    geom_vline(
      data = cutoffs_df,
      aes(xintercept = cutoff),
      color = "red",
      linetype = "dashed",
      size = 0.8
    ) +
    scale_x_log10(
      # If you have p-values down to 1e-13, set lower bound slightly smaller
      limits = c(1e-14, 1),  ### CHANGED
      breaks = c(1e-13, 1e-10, 1e-7, 1e-4, 1e-1, 1)
    ) +  ### CHANGED
    theme_minimal() +
    labs(
      title = paste("Permutation vs. True P-Values (Log Scale) -", chromatin_group),
      x = "Pearson p-value (log scale)",
      y = "Count",
      fill = "Source"
    )+
    scale_fill_manual(
      values = c("true_data" = "blue", "permutation_data" = "gray50"),   # ✅ ADD
      name = "Data Source"                                               # ✅ ADD
    )
  

  # --- 7. Combine and save  ### CHANGED
  final_plot <- p_linear / p_log  # stack them vertically

  ggsave(
    plot = final_plot,
    filename = file.path(
      OUTPUT_FOLDER_pearson_plot,
      paste0(chromatin_group, ".histogram_split_by_direction_combined.png")
    ),
    width = 12, height = 14  # bigger to accommodate two plots
  )

  
  ######### FDR analyis
  
  # Define cutoff thresholds
  cutoffs <- seq(0, 0.01, by = 0.00001)
  
  # Separate true p-values
  true_vals <- all_pvals_combined[file == "true_data", pval]
  
  # Subset to permutation values
  perm_vals_dt <- all_pvals_combined[file != "true_data"]
  
  # Get list of permutations
  perm_files <- unique(perm_vals_dt$file)
  
  # Initialize list to store results
  fdr_list <- lapply(perm_files, function(perm) {
    this_perm_vals <- perm_vals_dt[file == perm, pval]
    data.table(
      cutoff = cutoffs,
      true_hits = sapply(cutoffs, function(cut) sum(true_vals <= cut, na.rm = TRUE)),
      null_hits = sapply(cutoffs, function(cut) sum(this_perm_vals <= cut, na.rm = TRUE)),
      permutation = perm
    )
  })
  
  # Combine all into one long data.table
  fdr_dt <- rbindlist(fdr_list)
  
  # Compute FDR for each permutation
  fdr_dt[, FDR := ifelse(true_hits > 0, null_hits / true_hits, NA_real_)]
  
  # Plot FDR curves for each permutation
  library(ggplot2)
  
  p <- ggplot(fdr_dt, aes(x = cutoff, y = FDR, group = permutation)) +
    geom_line(alpha = 0.4, color = "gray40") +
    theme_minimal() +
    labs(
      title = paste("FDR Curves by Permutation -", chromatin_group),
      x = "p-value cutoff",
      y = "False Discovery Rate"
    ) +
    ylim(0, 1)
  
  ggsave(
    plot = p,
    filename = file.path(OUTPUT_FOLDER_pearson_plot,
                         paste0(chromatin_group,".FDR Curves by Permutation.png")),
    width = 10,
    height =10
  )
  
  
  ######### 1. Setup
  # Define log10 cutoffs: 10^-1 to 10^-5
  cutoffs <- cutoffs <- 10^seq(-2, -8, by = -0.1)


  
  # TRUE p-values
  true_vals <- all_pvals_combined[file == "true_data", pval]
  
  # PERMUTATION p-values
  perm_vals_dt <- all_pvals_combined[file != "true_data"]
  perm_files <- unique(perm_vals_dt$file)
  
  ######### 2. Compute true & null hits + FDR per permutation
  fdr_list <- lapply(perm_files, function(perm) {
    this_perm_vals <- perm_vals_dt[file == perm, pval]
    data.table(
      cutoff = cutoffs,
      true_hits = sapply(cutoffs, function(cut) sum(true_vals <= cut, na.rm = TRUE)),
      null_hits = sapply(cutoffs, function(cut) sum(this_perm_vals <= cut, na.rm = TRUE)),
      permutation = perm
    )
  })
  
  fdr_dt <- rbindlist(fdr_list)
  fdr_dt[, FDR := ifelse(true_hits > 0, null_hits / true_hits, NA_real_)]
  
  ######### 3. Reshape for plotting
  # FDR curves
  fdr_dt[, type := "FDR"]
  fdr_plot <- fdr_dt[, .(cutoff, value = FDR, permutation, type)]
  
  # True hit counts
  true_plot <- fdr_dt[, .(cutoff, value = true_hits, permutation, type = "True Hits")]
  
  # Permuted hit counts
  null_plot <- fdr_dt[, .(cutoff, value = null_hits, permutation, type = "Permuted Hits")]
  
  # Combine all for plotting
  plot_df <- rbindlist(list(fdr_plot, true_plot, null_plot))
  plot_df[, type := factor(type, levels = c("FDR", "True Hits", "Permuted Hits"))]
  
  ######### 4. Plot
  p <- ggplot(plot_df, aes(x = cutoff, y = value, group = permutation)) +
    geom_line(alpha = 0.4) +
    scale_x_log10(
      breaks = cutoffs,
      labels = format(cutoffs, scientific = FALSE)
    ) +
    facet_wrap(~type, scales = "free_y", ncol = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste("FDR and Hit Counts Across Permutations -", chromatin_group),
      x = "p-value cutoff (log10 scale)",
      y = "Value"
    )
  
  ggsave(
    plot = p,
    filename = file.path(OUTPUT_FOLDER_pearson_plot,
                         paste0(chromatin_group,".FDR and Hit Counts Across Permutations.png")),
    width = 10,
    height =10
  )
  
  }
}
# done q-value additon pearson #####################################
if(FALSE) {
  min_delta <- 0.5
  
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(
      this.dir,
      "12.pipeline/results/WholeGenome/all_CpG_complete_with_test.45.rds"
    )
  )
  
  # Suppose your data frame is called df
  all_CpG_complete_with_test.45 <- all_CpG_complete_with_test.45[, !colnames(all_CpG_complete_with_test.45) %in% c("kw.p_val", "kw.statistic", "kw.delta", "pearson.p_val_BH")]
  
  
  OUTPUT_FOLDER_pearson <- file.path(
    OUTPUT_FOLDER,
    paste(
      "45.pearson",
      "CpG_permutation",
      "min_delta",
      min_delta,
      sep = "."
    )
  )
  dir.create(path = OUTPUT_FOLDER_pearson, showWarnings = FALSE)
  
  all_CpG_complete_with_test.45[, paste0("pearson.q_vals.min_delta_", min_delta)] <- NA
  n_permutations <- c()
  
  for (chromatin_group in levels(all_CpG_complete_with_test.45$state_group)) {
    # for debuging: chromatin_group <- levels(all_CpG_complete_with_test.45$state_group)[1]
    start_time <- Sys.time()
    print(chromatin_group)
    
    # Get real p-values in group with delta > threshold
    idx_real <- which(
      all_CpG_complete_with_test.45$state_group == chromatin_group &
        all_CpG_complete_with_test.45$pearson.delta > min_delta
    )
    real_pvals <- all_CpG_complete_with_test.45$pearson.p_val[idx_real]
    
    # Load permutation file
    sim_age_file_name <- list.files(
      path = OUTPUT_FOLDER_pearson,
      pattern = paste0("^", chromatin_group, ".pearson.CpG_permutation."),
      full.names = TRUE
    )
    
    print(sim_age_file_name)
    
    # Combine all permuted p-values from all files  
    perm_pvals <- unlist(lapply(sim_age_file_name, function(f) {
      readRDS(f)$data$pearson.p_val
    }))  
    
    print(length(perm_pvals))
    n_permutations[chromatin_group] <- length(perm_pvals)
    # Pre-sort permuted p-values once  ### CHANGED
    perm_pvals_sorted <- sort(perm_pvals)
    
    # Compute empirical FDR-adjusted q-values
    real_order <- order(real_pvals)
    sorted_real <- real_pvals[real_order]
    emp_fdr <- numeric(length(sorted_real))
    
    #  use findInterval to count permuted p-values efficiently
    num_perm_below <- findInterval(sorted_real, perm_pvals_sorted)  
    
    #  compute number of real p-values ≤ threshold by index
    num_real_below <- seq_along(sorted_real)  
    
    #vectorized FDR calculation
    emp_fdr <- ifelse(num_real_below == 0, 0, num_perm_below / num_real_below)  ### CHANGED
    
    
    # Monotonic adjustment cap value at 1
    emp_qvals_sorted <- pmin(1, rev(cummin(rev(emp_fdr))))  
    emp_qvals <- numeric(length(real_pvals))
    emp_qvals[real_order] <- emp_qvals_sorted
    
    # Insert back into main dataframe
    all_CpG_complete_with_test.45[idx_real, paste0("pearson.q_vals.min_delta_", min_delta)] <- emp_qvals
    duration <- difftime(Sys.time(), start_time, units = "secs")
    print(paste(
      "Finished",
      chromatin_group,
      "in",
      format(.POSIXct(duration, tz = "GMT"), "%H:%M:%S")
    ))
    
  }
  saveRDS(
    object =  all_CpG_complete_with_test.45,
    file = file.path(
      this.dir,
      "12.pipeline/results/WholeGenome/all_CpG_complete_with_test.45.qval.rds"
    )
  )

}

# directional q-value additon pearson #####################################
if(FALSE){
  for (chromatin_group in levels(all_CpG_complete_with_test.45$state_group)) {
    # for debuging: chromatin_group <- levels(all_CpG_complete_with_test.45$state_group)[1]
    start_time <- Sys.time()
    print(chromatin_group)
    
    # Get real p-values in group with delta > threshold
    idx_real <- which(
      all_CpG_complete_with_test.45$state_group == chromatin_group &
        all_CpG_complete_with_test.45$pearson.delta > min_delta
    )
    
    # Get real test statistics (direction info)
    real_stats <- all_CpG_complete_with_test.45$pearson.statistic[idx_real]
    # Get real p-values for the same subset:
    real_pvals <- all_CpG_complete_with_test.45$pearson.p_val[idx_real]
    
    
    # Get direction: TRUE = positive, FALSE = negative
    is_pos <- real_stats > 0
    is_neg <- real_stats < 0
    
    # Load corresponding permuted test statistics (for direction filtering)
    perm_pvals <- unlist(lapply(sim_age_file_name, function(f) {
      readRDS(f)$data$pearson.p_val
    }))
    
    # Load corresponding permuted test statistics (for direction filtering)
    perm_stats <- unlist(lapply(sim_age_file_name, function(f) {
      readRDS(f)$data$pearson.statistic
    }))
    
    # ---- Positive Direction FDR ----
    real_pos <- real_pvals[is_pos]
    perm_pvals_pos <- perm_pvals[perm_stats > 0]
    
    if (length(perm_pvals_pos) > 0 && length(real_pos) > 0) {
      perm_pvals_pos <- sort(perm_pvals_pos)
      real_order_pos <- order(real_pos)
      sorted_real_pos <- real_pos[real_order_pos]
      
      num_perm_below_pos <- findInterval(sorted_real_pos, perm_pvals_pos)
      num_perm_below_pos <- pmin(num_perm_below_pos, length(perm_pvals_pos))
      
      num_real_below_pos <- seq_along(sorted_real_pos)
      emp_fdr_pos <- ifelse(num_real_below_pos == 0, 0, num_perm_below_pos / num_real_below_pos)
      emp_qvals_sorted_pos <- pmin(1, rev(cummin(rev(emp_fdr_pos))))
      emp_qvals_pos <- numeric(length(real_pos))
      emp_qvals_pos[real_order_pos] <- emp_qvals_sorted_pos
      
      col_pos <- paste0("pearson.q_vals.min_delta_", min_delta, ".positive")
      if (!(col_pos %in% colnames(all_CpG_complete_with_test.45))) {
        all_CpG_complete_with_test.45[[col_pos]] <- NA
      }
      all_CpG_complete_with_test.45[idx_real[is_pos], col_pos] <- emp_qvals_pos
    }
    
    
    
    # ---- Negative Direction FDR ----
    real_neg <- real_pvals[is_neg]
    perm_pvals_neg <- perm_pvals[perm_stats < 0]
    
    if (length(perm_pvals_neg) > 0 && length(real_neg) > 0) {
      perm_pvals_neg <- sort(perm_pvals_neg)
      real_order_neg <- order(real_neg)
      sorted_real_neg <- real_neg[real_order_neg]
      
      num_perm_below_neg <- findInterval(sorted_real_neg, perm_pvals_neg)
      num_perm_below_neg <- pmin(num_perm_below_neg, length(perm_pvals_neg))  
      
      num_real_below_neg <- seq_along(sorted_real_neg)
      emp_fdr_neg <- ifelse(num_real_below_neg == 0, 0, num_perm_below_neg / num_real_below_neg)
      emp_qvals_sorted_neg <- pmin(1, rev(cummin(rev(emp_fdr_neg))))
      emp_qvals_neg <- numeric(length(real_neg))
      emp_qvals_neg[real_order_neg] <- emp_qvals_sorted_neg
      
      col_neg <- paste0("pearson.q_vals.min_delta_", min_delta, ".negative")
      if (!(col_neg %in% colnames(all_CpG_complete_with_test.45))) {
        all_CpG_complete_with_test.45[[col_neg]] <- NA
      }
      all_CpG_complete_with_test.45[idx_real[is_neg], col_neg] <- emp_qvals_neg
    }
    
  }  
# plot direction distributes:


min_delta <- 0.5
q_thresh <- 0.1

# Prepare data
df_plot <- all_CpG_complete_with_test.45 %>%
  filter(pearson.delta >= min_delta, !is.na(pearson.p_val), !is.na(pearson.statistic)) %>%
  mutate(
    direction = case_when(
      pearson.statistic > 0 ~ "Positive",
      pearson.statistic < 0 ~ "Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(direction), pearson.p_val > 0)

# Plot
p <- ggplot(df_plot, aes(x = pearson.p_val, fill = direction)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.5) +
  scale_x_log10() +
  facet_wrap(~ state_group, scales = "free_y") +
  geom_vline(xintercept = q_thresh, linetype = "dashed", color = "black") +  # <- added line
  labs(
    title = paste("Histogram of p-values (log scale) by State Group — Δmeth >=", min_delta),
    x = "Pearson p-value (log scale)",
    y = "Count",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(size = 10))

ggsave(plot = p,filename = file.path(OUTPUT_FOLDER_pearson_plot,"q_value Destribtuion by state.png"),width = 16,height = 16)

# Prepare plotting data
plot_data <- all_CpG_complete_with_test.45 %>%
  filter(pearson.delta >= min_delta) %>%
  mutate(
    direction = case_when(
      pearson.statistic > 0 ~ "Positive",
      pearson.statistic < 0 ~ "Negative",
      TRUE ~ NA_character_
    ),
    is_significant_pos = get(paste0("pearson.q_vals.min_delta_", min_delta, ".positive")) < q_thresh,
    is_significant_neg = get(paste0("pearson.q_vals.min_delta_", min_delta, ".negative")) < q_thresh
  ) %>%
  filter(!is.na(direction)) %>%
  group_by(state_group, direction) %>%
  summarise(
    total = n(),
    significant = if (unique(direction) == "Positive") {
      sum(is_significant_pos, na.rm = TRUE)
    } else {
      sum(is_significant_neg, na.rm = TRUE)
    },
    frac_significant = significant / total,
    .groups = "drop"
  )


# Plot
ggplot(plot_data, aes(x = state_group, y = frac_significant, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = paste("Fraction of Significant CpGs (q <", q_thresh, ", delta >=", min_delta, ")"),
    x = "Chromatin State Group",
    y = "Fraction Significant",
    fill = "Correlation Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# diplay examples of pearson correlation ########################################
if(FALSE){
  min_delta <- 0.5
  max_q <- 0.01
  
  load_variable_if_not_exists(
    variable_name = "all_CpG_complete_with_test.45",
    file_path = file.path(
      this.dir,
      "12.pipeline/results/WholeGenome/all_CpG_complete_with_test.45.qval.rds"
    )
  )

  
  # Create a new plotting column where zeros are replaced by a small value
  df_plot <- all_CpG_complete_with_test.45
  df_plot$qval_no_zero <- ifelse(df_plot$pearson.q_vals.min_delta_0.5 == 0, 1e-10, df_plot$pearson.q_vals.min_delta_0.5)
  
  indx <- which(
    df_plot$pearson.q_vals.min_delta_0.5 <=  max_q  &
      df_plot$pearson.delta > min_delta
    
  )
  
  df_plot <- df_plot[indx,]
  
  ggplot(df_plot, aes(x = qval_no_zero)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black", alpha = 0.8) +
    scale_x_log10() +
    labs(
      title = "Histogram of Empirical q-values (log10 scale, zeros replaced)",
      x = "Empirical q-value (log10)",
      y = "Count"
    ) +
    theme_minimal(base_size = 14)
  
  

  
  
  significant_CpGs <- df_plot
  
  significant_CpGs <- significant_CpGs[order(significant_CpGs$pearson.p_val),]
  
  plot_age_correlation(df_row = significant_CpGs[547064,],typisation = real_age_typisation)
  
  plot_age_correlation(df_row = significant_CpGs[10000,],typisation = real_age_typisation)
  
  plot_age_correlation_new(df_row = significant_CpGs[10000,],typisation = real_age_typisation)
  
  plot_age_correlation_type(df_row = significant_CpGs[1500,],typisation = real_age_typisation)
  
  
  # # Create a plotting-safe q-value column
  # df_plot <- all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$pearson.delta > min_delta,]
  # df_plot$qval_no_zero <- ifelse(df_plot$pearson.q_vals.min_delta_0.5 == 0, 1e-20, df_plot$pearson.q_vals.min_delta_0.5)
  # 
  # # Compute significance flag for coloring
  # df_plot$significant <- with(df_plot,
  #                             ifelse(pearson.q_vals.min_delta_0.5 <= max_q & pearson.delta > min_delta, "Significant", "Not Significant")
  # )
  # 
  # # Make state_group a factor to preserve order
  # df_plot$state_group <- factor(df_plot$state_group, levels = unique(df_plot$state_group))
  # 
  # # Manhattan-style plot
  # ggplot(df_plot, aes(x = state_group, y = -log10(qval_no_zero), color = significant)) +
  #   geom_jitter(width = 0.3, alpha = 0.5, size = 0.6) +
  #   scale_color_manual(values = c("Significant" = "firebrick", "Not Significant" = "gray60")) +
  #   labs(
  #     title = paste0("Manhattan-style plot of q-values per Chromatin State (min_delta = ", min_delta, ")"),
  #     x = "Chromatin State Group",
  #     y = expression(-log[10](q["value"])),
  #     color = "Status"
  #   ) +
  #   theme_minimal(base_size = 14) +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     plot.title = element_text(hjust = 0.5)
  #   )
  
  
  min_delta <- 0.5
  max_q <- 0.01
  
  df_plot <- all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$pearson.delta > min_delta,]
  
  # Prepare data: compute % significant per state_group
  bar_data <- df_plot  %>%
    mutate(significant = pearson.q_vals.min_delta_0.5 <= max_q) %>%
    group_by(state_group) %>%
    summarise(
      total = n(),
      significant_count = sum(significant, na.rm = TRUE),
      percent_significant = 100 * significant_count / total
    )
  
  # Make bar plot
  ggplot(bar_data, aes(x = state_group, y = percent_significant)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(
      title = paste0("% of Significant CpGs per Chromatin State (q_value ≤ ", max_q, ", Δ > ", min_delta, ")"),
      x = "Chromatin State Group",
      y = "% Significant CpGs"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  
}

# Done Calculate true KW  Data -----------------------------------------------------------
if (create_true_stat) {
 
  start_time <- Sys.time()
  df_peak_CpG <- all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$name != "NO_CHIP",]
  print(paste("# CpG in H3K27Ac peaks:", nrow(df_peak_CpG)))
  
  print(paste(
    "# CpG in H3K27Ac peaks with all samples :",
    sum(df_peak_CpG$score == nrow(real_age_typisation))
  ))
  print(paste("discard ", 100 * (
    1 - sum(df_peak_CpG$score == nrow(real_age_typisation)) / nrow(df_peak_CpG)
  ), "% of CpG positions"))
  # discard CpG with not full data
  df_peak_CpG_complete <-
    df_peak_CpG[nrow(real_age_typisation) == df_peak_CpG$score, ]
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # compute real data value kruskall_valis
  start_time <- Sys.time()
  df_kruskall_valis <-
    parallel_testing_kruskall_valis(df_peak_CpG_complete, real_food_typisation)
  end_time <- Sys.time()
  print(end_time - start_time)

  # save results
  df_peak_CpG_complete_with_test <-
    cbind(df_peak_CpG_complete,
          df_kruskall_valis)
  
   # df_peak_CpG_complete_with_test <- all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$name != "NO_CHIP",]
  
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
  
  for (s in as.character(meta$sample)) {
    print(s)
    meth <- as.numeric(df_peak_CpG_complete_with_test[, s])
    sw_normality <-
      shapiro.test(sample(
        x = meth,
        size = 1000,
        replace = FALSE
      ))
    print(sw_normality)
  }
  
  
} else{
  df_peak_CpG_complete_with_test <-
    readRDS(file = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds"))
}
# Explore true KW  Data ###########################
if(FALSE){
  
  df_peak_CpG_complete_with_test_min_delta <- df_peak_CpG_complete_with_test[
    df_peak_CpG_complete_with_test$pearson.delta > min_delta,]
  
  print(paste0(nrow(df_peak_CpG_complete_with_test_min_delta ),
               " CpG in peak with min_delta ",
               min_delta))
  
  print(paste0(sum(df_peak_CpG_complete_with_test_min_delta$pearson.q_vals.min_delta_0.5 < max_q ),
               " CpG in peak with q-value smaller  ",
               max_q))
  
  print(paste0(sum(df_peak_CpG_complete_with_test_min_delta$pearson.q_vals.min_delta_0.5.positive < max_q ,na.rm = TRUE),
               " CpG in peak with positive directional q-value smaller  ",
               max_q))
  
  
  print(paste0(sum(df_peak_CpG_complete_with_test_min_delta$pearson.q_vals.min_delta_0.5.negative < max_q ,na.rm = TRUE),
               " CpG in peak with negtaive directional q-value smaller  ",
               max_q))
  
  
  
  df_peak_CpG_complete_with_test_sorted <- df_peak_CpG_complete_with_test[
    order(df_peak_CpG_complete_with_test$kw.p_val),]
  
  plot_food_correlation(df_row = df_peak_CpG_complete_with_test_sorted[1000,] ,
                        typisation = real_food_typisation)
  
  
  
}

# check other publications: describe LOCUS ---------------------------------------------------------------
if (FALSE) {
  
  study_comparison <- as.data.frame(matrix(ncol = 5,nrow = 5))
  colnames(study_comparison) <- c("publication",
                                  "origin",
                                  "# CpG", 
                                  "p.val comparison with tempral shift pearson test",
                                  "p.val comparison with Kruskal Wallis lifestly test")
  
  # all_CpG_complete_with_test.45 <- readRDS("12.pipeline/results/WholeGenome/all_CpG_complete_with_test.45.rds")
  
  #https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-019-0680-7/tables/2
  library(dplyr)
  library(GenomicRanges)
  library(openxlsx)
  
  # Load the data barbados
  df_barbados <- read.table(file = file.path("other_pubs", "Barbados", "barbadosLiftToHg19.bed"), header = FALSE, sep = "\t")
  colnames(df_barbados) <- c("chr", "start", "end")
  
  tests_barbados <- test_genomic_subsets(df_barbados, all_CpG_complete_with_test.45)
  
  study_comparison[1,] <- c("Peter et al., 2016","malnutrition on Barbados",418,tests_barbados$pearson_test_result$p.value,tests_barbados$kw_test_result$p.value)
  
  # Load the data dutch
  df_dutch <- read.table(file = file.path("other_pubs", "Dutch hunger", "dutch.hunger.liftTo.hg19.bed"), header = FALSE, sep = "\t")
  colnames(df_dutch) <- c("chr", "start", "end","name")
  tests_dutch <- test_genomic_subsets(df_dutch, all_CpG_complete_with_test.45)
  
  study_comparison[2,] <- c("Tobi et al., 2014","Dutch hunger",11697,tests_dutch$pearson_test_result$p.value,tests_dutch$kw_test_result$p.value)
  
  
  # Load the data bangladesh
  df_bangladesh <-
    read.xlsx(
      xlsxFile = file.path(
        "other_pubs",
        "rural Bangladesh",
        "bmjopen-2016-November-6-11 - inline-supplementary-material-1.xlsx"
      ),
      sheet = "significant"
    )
  
  df_bangladesh$start <- df_bangladesh$MAPINFO -1 
  df_bangladesh$end <- df_bangladesh$MAPINFO +1 
  names(df_bangladesh)[names(df_bangladesh) == "chr"] <- "CHR"
  
  tests_bangladesh <- test_genomic_subsets(df_bangladesh, all_CpG_complete_with_test.45)
  
  study_comparison[3,] <- c("Finer et al., 2016","hunger bangladesh",57,tests_bangladesh $pearson_test_result$p.value,tests_bangladesh$kw_test_result$p.value)
  
  # Load the data chinese
  df_chineses <-
    read.xlsx(
      xlsxFile = file.path(
        "other_pubs",
        "Chinese famine",
        "13148_2019_680_MOESM2_ESM (1).xlsx"
      ),
      sheet = "cpg rank Chinese"
    )
  
  df_chineses$start <- df_chineses$pos -1 
  df_chineses$end <- df_chineses$pos +1 
  df_chineses$CHR <- paste0("chr",df_chineses$CHR)
  
  #names(df_chineses)[names(df_chineses) == "chr"] <- "CHR"
  
  test_chinese <- test_genomic_subsets(df_chineses, all_CpG_complete_with_test.45)
  
  study_comparison[4, ] <-
    c(
      "He et al., 2019a",
      "chinese great hunger",
      18778,
      test_chinese$pearson_test_result$p.value,
      test_chinese$kw_test_result$p.value
    )
  
  
  # Load the data fibrolast
  df_fibrolast <-
    read.xlsx(
      xlsxFile = file.path(
        "other_pubs",
        "Chinese famine",
        "13148_2019_680_MOESM2_ESM (1).xlsx"
      ),
      sheet = "cpg rank Fibroblasts"
    )
  
  df_fibrolast$start <- df_fibrolast$pos -1 
  df_fibrolast$end <- df_fibrolast$pos +1 
  df_fibrolast$CHR <- paste0("chr",df_fibrolast$CHR)
  
  test_fibrolast <- test_genomic_subsets(df_fibrolast, all_CpG_complete_with_test.45)
  
  study_comparison[5, ] <-
    c("He et al., 2019a",
      "deprived fibrolast",
      56151,
      test_fibrolast$pearson_test_result$p.value,
      test_fibrolast$kw_test_result$p.value
    )
  
  study_comparison$bonferroni.correction <-  pmin(as.numeric(study_comparison$`p.val comparison with tempral shift pearson test`) * nrow(study_comparison), 1)
  study_comparison$`# CpG` <- as.numeric(study_comparison$`# CpG`)
  study_comparison$`p.val comparison with tempral shift pearson test` <- as.numeric(  study_comparison$`p.val comparison with tempral shift pearson test`)
  write.xlsx(study_comparison, file.path(
    "other_pubs","TableS4.xlsx"))
  
  # Display the first few rows of the DataFrame
  head(df_bangladesh)

  GH <-
    list(
      ID = "GH05J143394",
      chr = "chr5",
      start = 142774324,
      end = 142785910
    )
  Gene <- "NR3C1"
  H3K27ac <- "chr5.142784176.142785514.H3K27ac"
  
  
  df <-
    all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$chrom == GH$chr &
                                 GH$start - 1000 <= all_CpG_complete_with_test.45$start &
                                 all_CpG_complete_with_test.45$end <= GH$end + 1000, ]
  
  unique(df$name)
  
  loc1 <-
    data.frame(name = "up Risk of being in a poorly regulated neurobehavioural profile",
               start =  142783501,
               end = 142783640)
  
  loc2 <-
    data.frame(name = "up Choline",
               start = 142783501,
               end = 142783908)
  
  loc3 <-
    data.frame(name = "up Adult waist circumference, up adult BMI",
               start = 142782759,
               end = 142783164)
  
  loc4 <-
    data.frame(name = "up Meat/fish and vegetable intake, down bread/potato intake in late pregnancy",
               start = 142783579,
               end = 142783714)
  
  loc5 <-
    data.frame(name = "up Adult blood pressure",
               start = 142783578,
               end = 142783714)
  
  xmin <- 142782200
  xmax <- 142785600
  
  p <-
    ggplot(data = df, mapping = aes(x = start, y = -log(pearson.p_val))) +
    geom_rect(
      aes(
        xmin = max(GH$start, xmin),
        xmax = min(xmax, GH$end),
        ymin = -1,
        ymax = -0.5
      ),
      fill = "yellow",
      alpha = 0.9
    ) +
    geom_segment(
      aes(
        x = xmax,
        y = -0.7,
        xend = xmax + 100 ,
        yend = -0.7
      ),
      arrow = arrow(length = unit(0.3, "cm")),
      color = "yellow"
    ) +
    geom_segment(
      aes(
        x = xmin,
        y = -0.7,
        xend = xmin + 100 ,
        yend = -0.7
      ),
      arrow = arrow(length = unit(0.3, "cm")),
      color = "yellow"
    ) +
    geom_rect(
      aes(
        xmin = 142782230,
        xmax = 142783990,
        ymin = -0.5,
        ymax = -0.1
      ),
      fill = "LightSkyBlue",
      alpha = 0.9
    ) +
    
    geom_rect(
      aes(
        xmin = 142784176,
        xmax = 142785514,
        ymin = -0.5,
        ymax = -0.1
      ),
      fill = "LightSkyBlue",
      alpha = 0.9
    ) +
    geom_bar(stat = "identity") +
    geom_text(aes(
      x = (xmin + xmax) / 2,
      y = -0.75,
      label = GH$ID
    ), color = "black") +
    geom_text(aes(
      x = (142782230 + 142783990) / 2,
      y = -0.3,
      label = "H3K27Ac"
    ), color = "black") +
    geom_text(aes(
      x = (142784176 + 142785514) / 2,
      y = -0.3,
      label = "H3K27Ac"
    ), color = "black") +
    geom_bar(data = df[which.min(df$pearson.p_val), ],
             stat = "identity",
             color = "red") +
    ylab(expression(-log(pearson_p[value]))) +
    xlab("Chr5: Genomic position hg19") +
    xlim(c(xmin, xmax)) +
    theme_minimal() +
    # geom_rect(aes(xmin = loc1$start, xmax = loc1$end, ymin = -1.5, ymax = -1.4), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc2$start, xmax = loc2$end, ymin = -1.4, ymax = -1.3), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc3$start, xmax = loc3$end, ymin = -1.3, ymax = -1.2), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc4$start, xmax = loc4$end, ymin = -1.2, ymax = -1.1), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc5$start, xmax = loc5$end, ymin = -1.1, ymax = -1.0), fill = "green", alpha = 0.9)+
    coord_cartesian(expand = c(0, 0))
  
  ggsave(
    plot = p,
    filename = file.path(
      OUTPUT_FOLDER,
      "results",
      "pearson",
      paste(Gene, GH$ID, "png", sep = ".")
    )
  )
  
  # https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-019-0680-7/tables/2
  DMR <-
    list(
      ID = "ENO2",
      chr = "chr12",
      start = 7023752,
      end = 7024121
    )
  DMR <-
    list(
      ID = "ZNF226",
      chr = "chr19",
      start = 44669146,
      end = 44669354
    )
  DMR <-
    list(
      ID = "CCDC51/TMA7",
      chr = "chr3",
      start = 48481268,
      end = 48481793
    )
  
  df <-
    all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$chrom == DMR$chr &
                                 DMR$start - 1000 <= all_CpG_complete_with_test.45$start &
                                 all_CpG_complete_with_test.45$end <= DMR$end + 1000, ]
  
  xmin <- DMR$start - 1000
  xmax <- DMR$end + 1000
  
  unique(df$name)
  
  p <-
    ggplot(data = df, mapping = aes(x = start, y = -log(pearson.p_val))) +
    geom_rect(
      aes(
        xmin = max(GH$start, xmin),
        xmax = min(xmax, GH$end),
        ymin = -1,
        ymax = -0.5
      ),
      fill = "yellow",
      alpha = 0.9
    ) +
    # geom_segment(aes(x = xmax, y = -0.7, xend = xmax + 100 , yend = -0.7), arrow = arrow(length = unit(0.3, "cm")), color = "yellow")+
    # geom_segment(aes(x = xmin, y = -0.7, xend = xmin + 100 , yend = -0.7), arrow = arrow(length = unit(0.3, "cm")), color = "yellow")+
    geom_rect(
      aes(
        xmin = 48481218,
        xmax = 48482176,
        ymin = -0.5,
        ymax = -0.1
      ),
      fill = "LightSkyBlue",
      alpha = 0.9
    ) +
    #geom_rect(aes(xmin = 142784176, xmax = 142785514, ymin = -0.5, ymax = -0.1), fill = "LightSkyBlue", alpha = 0.9)+
    geom_bar(stat = "identity") +
    #geom_text(aes(x = (xmin + xmax) / 2, y = -0.75, label = GH$ID), color = "black") +
    #geom_text(aes(x = (142782230 + 142783990) / 2, y = -0.3, label = "H3K27Ac"), color = "black") +
    #geom_text(aes(x = (142784176 + 142785514) / 2, y = -0.3, label = "H3K27Ac"), color = "black")+
    #geom_bar(data = df[which.min(df$pearson.p_val),],stat = "identity",color = "red")+
    ylab(expression(-log(pearson_p[value]))) +
    xlab(paste(DMR$chr, ": Genomic position hg19")) +
    xlim(c(xmin, xmax)) +
    theme_minimal() +
    # geom_rect(aes(xmin = loc1$start, xmax = loc1$end, ymin = -1.5, ymax = -1.4), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc2$start, xmax = loc2$end, ymin = -1.4, ymax = -1.3), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc3$start, xmax = loc3$end, ymin = -1.3, ymax = -1.2), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc4$start, xmax = loc4$end, ymin = -1.2, ymax = -1.1), fill = "green", alpha = 0.9)+
    # geom_rect(aes(xmin = loc5$start, xmax = loc5$end, ymin = -1.1, ymax = -1.0), fill = "green", alpha = 0.9)+
    coord_cartesian(expand = c(0, 0))
  
}


# Done describe Data ---------------------------------------------------------------
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
  min_delta = 0.5
  paste(
    "number of CpG in peaks with data with min_delta:",
    sum(df_peak_CpG_complete_with_test$pearson.delta > min_delta)
  )
  
  # p <-
  #   ggplot(data = meta,
  #          mapping = aes(x = sample, y = age_mean_BP, color = Type)) +
  #   geom_point() +
  #   geom_errorbar(aes(ymin = age_mean_BP - age_std_BP, ymax = age_mean_BP +
  #                       age_std_BP),
  #                 width = .1) +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(
  #     angle = 0,
  #     vjust = 0.5,
  #     hjust = 1
  #   )) +
  #   scale_y_continuous(
  #     breaks = seq(0, 11000, by = 1000),
  #     # labels = rep("",13),
  #     limits = c(0, 11500),
  #     expand = c(0, 0)
  #   ) +
  #   ylab("years before present (B.P.)") +
  #   xlab("sample ID") +
  #   scale_color_discrete(name = "lifestyle/diet") +
  #   theme(
  #     strip.text.x = element_blank(),
  #     strip.background = element_rect(colour = "white", fill = "white"),
  #     legend.position = c(.9, .4)
  #   ) +
  #   coord_flip()
  # 
  # ggsave(plot = p,filename = "C:/Users/Batyrev/Dropbox/paleo epigenetics R10/FigureS1.png")
  
  
  
  
  # Figure S1
  p <- ggplot(data = meta, mapping = aes(x = sample, y = age_mean_BP, color = Type)) +
    geom_point(size = 3, alpha = 0.6) +  # Enhanced points for better visibility
    geom_errorbar(aes(ymin = age_mean_BP - age_std_BP, ymax = age_mean_BP + age_std_BP),
                  width = 0.1, size = 0.5) +  # Adjusted error bars for clarity
    theme_minimal(base_size = 14) +  # Base size for readability
    theme(
      axis.text.x = element_text(vjust = 1, hjust = 1),  # Adjusting x-axis labels for better readability
      legend.position = c(.9, .4),
      legend.title = element_text(face = "bold"),  # Bold legend title for emphasis
      plot.title = element_text(size = 16, face = "bold"),  # Enhanced plot title
      axis.title.x = element_text(size = 14, face = "bold"),  # Bold and larger x-axis label
      axis.title.y = element_text(size = 14, face = "bold")  # Bold and larger y-axis label
    ) +
    scale_y_continuous(
      breaks = seq(0, 11000, by = 1000),  # Adjusted breaks for y-axis
      limits = c(-200, 11500),
      expand = c(0, 0)
    ) +
    ylab("Years Before Present (B.P.)") +
    xlab("Sample ID") +
    scale_color_discrete(name = "Lifestyle") +  # Clarified legend title
    coord_flip()  # Flipping coordinates for horizontal layout
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "plots",
                         "Age of samples.png"),
    plot = p,
    width = 12,
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
  
  x_max <- 11000
  
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
                       limits = c(0, x_max),) +
    scale_y_continuous(limits = c(-3, 3)) +
    theme(
      panel.background = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ) + xlab("years before present (B.P.)") +
    scale_color_discrete(name = "lifestyle/diet")
  
  p <-
    ggplot(data = meta , aes(x = age_mean_BP, y = 0, fill = Type)) +
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
                       limits = c(0, x_max),) +
    scale_y_continuous(limits = c(-3, 3)) +
    theme(
      panel.background = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12)
    ) + xlab("years before present (B.P.)") +
    scale_fill_manual(name = "lifestyle/diet",
                      values = mean_age_by_type$color,
                      labels = mean_age_by_type$Type)
  
  ggsave(
    filename = file.path(OUTPUT_FOLDER,
                         "plots",
                         "Age of samples 1D.png"),
    plot = p,
    width = 12,
    height = 4
  )

}
#"#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
# Calculate column permutation simulations --------------------------------------------------

if (create_permutations_horizonal_columns) {
  print("create_permutations_horizonal_columns")
  dir.create(path = file.path(OUTPUT_FOLDER_pearson_max_q, "simulation"),
             showWarnings = FALSE)
  
  #df_peak_CpG_complete_with_test_min_delta
  # best_CpGs
  
  n_repetitions <- 100
  for (rep in 1:n_repetitions) {
    # rep <- 1
    start_time <- Sys.time()
    print(rep)
    
    permutation_age  <-
      create_horizontal_permutated_typisation(real_age_typisation)
    permutation_food  <-
      create_horizontal_permutated_typisation(real_food_typisation)
    
    sim_age_typisation <- permutation_age$sim_typisation
    sim_food_typisation <- permutation_food$sim_typisation
    
    permutation_age$data <-
      parallel_testing_pearson_cor(df = df_peak_CpG_complete_with_test_min_delta, age.typisation = sim_age_typisation)
    
    permutation_food$data <-
      parallel_testing_kruskall_valis(df = df_peak_CpG_complete_with_test_min_delta, food.typisation = sim_food_typisation)
    
    # save age permutation
    sim_age_file_name <-
      paste(
        "pearson",
        "horizontal",
        paste(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), collapse = "_"),
        "rds",
        sep = "."
      )
    saveRDS(
      object = permutation_age,
      file = file.path(OUTPUT_FOLDER_pearson_max_q, "simulation", sim_age_file_name)
    )
    
    # save food permutation
    sim_food_file_name <-
      paste(
        "KW",
        "horizontal",
        paste(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), collapse = "_"),
        "rds",
        sep = "."
      )
    saveRDS(
      object = permutation_food,
      file = file.path(OUTPUT_FOLDER_pearson_max_q, "simulation", sim_food_file_name)
    )
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
}
if(FALSE){
  # --- Load libraries
  library(parallel)
  library(pbapply)
  
  # --- Parameters
  total_permutations <- 1e6
  chunk_size <- 1e5
  n_chunks <- total_permutations / chunk_size
  n_cores <- detectCores() - 1
  set.seed(42)
  
  # --- Prepare input
  meth_matrix <- best_CpGs[, as.character(meta$sample), drop = FALSE]
  age_original <- meta$age_mean_BP[match(colnames(meth_matrix), meta$sample)]
  
  # --- Loop over chunks
  for (chunk_id in 1:n_chunks) {
    cat("Running chunk", chunk_id, "of", n_chunks, "\n")
    
    # Precompute permuted ages for this chunk
    age_permutations <- replicate(chunk_size, sample(age_original), simplify = FALSE)
    
    # Setup parallel cluster
    cl <- makeCluster(n_cores)
    clusterExport(cl, varlist = c("meth_matrix", "age_permutations"), envir = environment())
    
    # Run label permutations with progress
    perm_pval_list <- pbapply::pblapply(seq_len(nrow(meth_matrix)), cl = cl, FUN = function(i) {
      meth <- as.numeric(meth_matrix[i, ])
      if (sum(!is.na(meth)) < 3) return(rep(NA, length(age_permutations)))
      
      sapply(age_permutations, function(age_perm) {
        res <- suppressWarnings(cor.test(meth, age_perm, method = "pearson"))
        res$p.value
      })
    })
    
    stopCluster(cl)
    
    # Convert to data frame
    perm_pvals_df <- do.call(rbind, perm_pval_list)
    rownames(perm_pvals_df) <- best_CpGs$name
    
    # Save this chunk
    chunk_file <- file.path(
      OUTPUT_FOLDER_pearson_max_q,
      paste0("CpG_label_permutation_chunk_", chunk_id, "_of_", n_chunks, ".rds")
    )
    
    saveRDS(perm_pvals_df, chunk_file)
    cat("✅ Saved:", chunk_file, "\n\n")
    
    # Clean up memory
    rm(age_permutations, perm_pval_list, perm_pvals_df)
    gc()
  }
  
  
  # # --- Input
  # i <- 1  # Choose CpG row index
  # real_pval <- best_CpGs$pearson.p_val[i]
  # null_pvals <- perm_pvals_df[i, ]
  # 
  # # --- Plot
  # library(ggplot2)
  # 
  # df_plot <- data.frame(pval = as.numeric(null_pvals))
  # 
  # ggplot(df_plot, aes(x = pval)) +
  #   geom_histogram(binwidth = 0.2, fill = "gray70", color = "black", alpha = 0.7) +
  #   geom_vline(xintercept = real_pval, color = "red", linetype = "dashed", size = 1) +
  #   scale_x_log10(
  #     limits = c(1e-14, 1),
  #     breaks = c(1e-12, 1e-9, 1e-6, 1e-3, 1)
  #   ) +
  #   theme_minimal(base_size = 14) +
  #   labs(
  #     title = paste("Null Distribution of Permuted p-values (CpG:", best_CpGs$name[i], ")"),
  #     subtitle = paste("Real p-value =", signif(real_pval, 3)),
  #     x = "Permuted p-values (log10)",
  #     y = "Count"
  #   )
  
  # --- Prepare
  chunk_files <- list.files(
    path = OUTPUT_FOLDER_pearson_max_q,
    pattern = "^CpG_label_permutation_chunk_\\d+_of_\\d+\\.rds$",
    full.names = TRUE
  )
  
  n_chunks <- length(chunk_files)
  n_CpGs <- nrow(best_CpGs)
  real_pvals <- best_CpGs$pearson.p_val
  
  # --- Initialize counts
  count_below <- numeric(n_CpGs)
  total_permuted <- 0
  
  # --- Process each chunk
  for (file in chunk_files) {
    cat("Processing:", file, "\n")
    chunk_df <- readRDS(file)
    
    stopifnot(nrow(chunk_df) == n_CpGs)  # sanity check
    total_permuted <- total_permuted + ncol(chunk_df)
    
    # Count how many permuted p-values are ≤ real p-value
    for (i in seq_len(n_CpGs)) {
      count_below[i] <- count_below[i] + sum(chunk_df[i, ] <= real_pvals[i], na.rm = TRUE)
    }
    
    rm(chunk_df)
    gc()
  }
  
  # --- Compute empirical p-values
  empirical_pvals <- (count_below + 1) / (total_permuted + 1)
  
  # --- Add to best_CpGs
  best_CpGs$empirical_p_val <- empirical_pvals
  
  # --- Add FDR-corrected q-values for empirical p-values
  best_CpGs$empirical_q_val <- p.adjust(best_CpGs$empirical_p_val, method = "BH")
  
  max_q <- 1e-05
  
  ggplot(best_CpGs, aes(x = pearson.p_val, y = empirical_p_val)) +
    geom_point(alpha = 0.5, size = 1.8, color = "steelblue") +
    scale_x_log10(
      limits = c(1e-14, 1),
      breaks = c(1e-12, 1e-9, 1e-6, 1e-3, 1)
    ) +
    scale_y_log10(
      limits = c(1e-14, 1),
      breaks = c(1e-12, 1e-9, 1e-6, 1e-3, 1)
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "orange") +
    geom_vline(xintercept = max_q,linetype = "dashed", color = "green")+
    geom_hline(yintercept = max_q, linetype = "dashed", color = "red")+
    theme_minimal(base_size = 14) +
    labs(
      title = "Empirical p-values vs Pearson p-values",
      x = "Pearson p-value",
      y = "Empirical p-value"
    )
  
  
  ggplot(best_CpGs, aes(x = pearson.p_val, y = empirical_q_val)) +
    geom_point(alpha = 0.5, size = 1.8, color = "steelblue") +
    scale_x_log10(
      limits = c(1e-14, 1),
      breaks = c(1e-12, 1e-9, 1e-6, 1e-3, 1)
    ) +
    scale_y_log10(
      limits = c(1e-14, 1),
      breaks = c(1e-12, 1e-9, 1e-6, 1e-3, 1)
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "orange") +
    geom_vline(xintercept = max_q,linetype = "dashed", color = "green")+
    geom_hline(yintercept = max_q, linetype = "dashed", color = "red")+
    theme_minimal(base_size = 14) +
    labs(
      title = "Empirical q-values vs Pearson p-values",
      x = "Pearson p-value",
      y = "Empirical p_val"
    )
  

  best_significant_CpGs <- best_CpGs[best_CpGs$empirical_q_val <max_q,]
  
  best_significant_CpGs <- best_significant_CpGs[order(best_significant_CpGs$pearson.p_val),]
  
  plot_age_correlation_type(df_row = best_significant_CpGs[1,],typisation = real_age_typisation)
  
  plot_age_correlation_type(df_row = best_significant_CpGs[nrow(best_significant_CpGs),],typisation = real_age_typisation)
  
  
}
# Calculate CpG permuation simulations --------------------------------------------------
if (create_CpG_permutations_vertical) {
  print("create_CpG_permutations_vertical")
  dir.create(path = file.path(OUTPUT_FOLDER, paste("45.pearson","CpG_permutation",sep = ".")),
             showWarnings = FALSE)
  
  dir.create(path = file.path(OUTPUT_FOLDER, paste("45.KW","CpG_permutation",sep = ".")),
             showWarnings = FALSE)
  
  n_repetitions <- 2
  for (rep in 1:n_repetitions) {
    # rep <- 1
    start_time <- Sys.time()
    print(rep)
    
    time_string <- format(start_time, "%Y_%m_%d_%H_%M_%S")
    
    #peak cpg
    permutation_CpG <-
      create_CpG_permution(df_peak_CpG_complete_with_test, real_age_typisation)
    
    
    # permutation_CpG <-
    #   create_CpG_permution(all_CpG_complete_with_test.45, real_age_typisation)
    
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
      file = file.path(OUTPUT_FOLDER, paste("pearson","CpG_permutation",sep = "."), sim_age_file_name)
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
      file = file.path(OUTPUT_FOLDER, paste("KW","CpG_permutation",sep = "."), sim_food_file_name)
    )
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
}
# Evaluate CpG simulations ----------------------------------------------------
if (summerize_CpG_columns_permutations) {
  print("summerize_CpG_columns_permutations")
  # list.files(path =  file.path(OUTPUT_FOLDER, "simulation"))
  
  test <- "pearson"
  #test <- "KW"
  
  file_names <-
    list.files(
      path = file.path(OUTPUT_FOLDER, paste(test,"CpG_permutation",sep = ".")),
      pattern = paste(test, "*.CpG_permutation.*", sep = ".")
    )
  print(length(file_names))
  
  #print(file.path(OUTPUT_FOLDER, paste(test,"CpG_permutation",sep = "."), file_names))
  
  start_time <- Sys.time()
  sim_results <-
    parallel_summerize_permutations(sim_file_names = file.path(OUTPUT_FOLDER, paste(test,"CpG_permutation",sep = "."), file_names))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  saveRDS(
    object = sim_results,
    file = file.path(OUTPUT_FOLDER, paste(test,"CpG_permutation.sim_results.rds",sep="."))
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
  
  sim_results_pearson_CpG <-
    readRDS(
      "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/pearson.CpG_permutation.sim_results.cluster.rds"
    )
  
  sim_results_pearson_horizontal <-
    readRDS(
      "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/pearson.horizontal.sim_results_cluster.rds"
    )
  
  
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
    sim_results_pearson_CpG$minus_log_alpha,-sim_results_pearson_CpG$minus_log_alpha
  )
  df_plot_horizontal$alpha <-
    c(
      sim_results_pearson_horizontal$minus_log_alpha,-sim_results_pearson_horizontal$minus_log_alpha
    )
  
  # prepare real resluts
  real_results_pearson <-
    real_results[real_results$test == "pearson", ]
  
  
  df_plot_real <-
    rbind(real_results_pearson[, c("min_delta" , "n_signfincant_CpG", "name")],
          real_results_pearson[, c("min_delta" , "n_signfincant_CpG", "name")])
  
  df_plot_real$significant <-
    c(
      real_results_pearson$positive_signfincant_CpG,
      real_results_pearson$negative_signfincant_CpG
    )
  
  df_plot_real$alpha <-
    c(real_results_pearson$minus_log_alpha,-real_results_pearson$minus_log_alpha)
  
  df_plot_real_1 <- df_plot_real[df_plot_real$min_delta == delta, ]
  df_plot_1 <- df_plot[df_plot$min_delta == delta, ]
  df_plot_horizontal_1 <-
    df_plot_horizontal[df_plot_horizontal$min_delta == delta, ]
  
  
  p <-
    ggplot(data = df_plot_1[abs(as.numeric(df_plot$alpha)) > cut_off, ] ,
           mapping = aes(x = alpha, y = significant, group = alpha)) +
    geom_boxplot(color = "green", outlier.shape = NA) +
    geom_boxplot(
      data = df_plot_horizontal_1[abs(as.numeric(df_plot_horizontal_1$alpha)) > cut_off, ],
      mapping = aes(x = alpha, y = significant, group = alpha),
      color = "blue",
      outlier.shape = NA
    ) +
    geom_point(data = df_plot_real_1[abs(as.numeric(df_plot_real_1$alpha)) > cut_off, ], color = "red") +
    theme_minimal() +
    scale_color_discrete(#name = "Dose",
      labels = c("CpG permuatation", "age permuation", "real data"))
  
  folder_name <-
    paste("delta", delta, "minLogP", alpha_p, sep = "_")
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
                              as.numeric(df_plot_all$alpha) > cut_off, ],
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
      "# of significant CpGs:  ", p[value] <= e ^ (-9.3) , ", ", delta[meth] > 0.39
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
  
  
# Adjusting the plot to differentiate the unpermuted data point
p <- ggplot(data = df_plot_single, aes(x = slope, y = significant, colour = name)) +
  geom_point(data = df_plot_single[df_plot_single$name != "pearson.real", ], # Plotting only permuted data with jitter
             aes(x = slope, y = significant, colour = name),
             position = position_jitterdodge(jitter.width = 0.25), size = 1, stroke = 0.5) +
  geom_point(data = df_plot_single[df_plot_single$name == "pearson.real", ], # Plotting unpermuted data without jitter
             aes(x = slope, y = significant, colour = name),
             size = 3, stroke = 1) + # Making the unpermuted data point larger
  theme_minimal(base_size = 12) +
  ylab(expression(paste("# of significant CpGs:  ", p[value] <= e ^ (-9.3) , ", ", delta[meth] > 0.39))) +
  xlab("Pearson correlation slope") +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1"),
                     name = "Dataset Type",
                     labels = c("CpG Position Permutation", "Sample Age Permutation", "Unpermuted Data")) +
  theme(legend.position = c(0.95, 0.95), # Coordinates for top right corner
        legend.justification = c("right", "top"), # Adjusts the anchor point of the legend
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.background = element_blank(), # Optionally, make legend background transparent
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14))


# Assuming OUTPUT_FOLDER and folder_name are predefined
file_name <- "count_distribution_age_and_CpG_permutation.png"

# Saving the plot
ggsave(filename = file.path(OUTPUT_FOLDER, "results", "pearson", folder_name, file_name),
       plot = p, width = 6, height = 6)
  
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
  
  #' p <-
  #'   ggplot(
  #'     data = df_plot_simmulations[df_plot_simmulations$name != "pearson.CpG_permutation", ],
  #'     mapping = aes(
  #'       x = alpha_string,
  #'       y = significant,
  #'       color = interaction(name, slope)
  #'     )
  #'   ) +
  #'   geom_boxplot(outlier.shape = NA) +
  #'   xlab(label = expression("significance cut off":alpha)) +
  #'   ylab(label = "number of CpG passing the significance threshold") +
  #'   scale_color_manual(
  #'     values = interaction_colors,
  #'     name = expression("Number of p-values" < e ^ -alpha),
  #'     #"number of p values < exp(-alpha)",
  #'     labels = c(
  #'       #'CpG permutation with negative pearson correlation',
  #'       'Age permutation with negative pearson correlation',
  #'       'data with negative pearson correlation',
  #'       #'CpG permutation with positive pearson correlation',
  #'       'Age permutation with positive pearson correlation',
  #'       'data with positive pearson correlation'
  #'     )
  #'   ) +
  #'   theme_minimal() +
  #'   theme(legend.position = c(.8, .5))+
    
    ggplot(
      data = df_plot_simmulations[!df_plot_simmulations$name == "pearson.CpG_permutation", ],
      mapping = aes(
        x = alpha_string,
        y = significant,
        color = interaction(name, slope)  # Facilitates legend display
      )
    ) +
    geom_boxplot(outlier.shape = NA) +  # Hides outliers
    xlab(label = expression(alpha ~ " significance cut-off")) +  # Correct expression usage
    ylab(label = "Number of CpG passing the significance threshold") +  # Improved label description
    scale_color_manual(
      values = interaction_colors,
      name = expression("Number of p-values" < e^-alpha),  # Correct expression formatting for consistency
      labels = c(
        'Age permutation with negative Pearson correlation',
        'Data with negative Pearson correlation',
        'Age permutation with positive Pearson correlation',
        'Data with positive Pearson correlation'
      )
    ) +
    theme_minimal() +  # Clean minimalistic theme
    theme(
      legend.position = c(0.8, 0.5)  # Adjust legend position to be inside the plot area
    )+
    ylim(c(0,800))
  
  
    # Base colors from RColorBrewer
    base_colors <- brewer.pal(3, "Set1")
    
    # Create extended palette: adjust transparency or use `scales::darken/lighten` for variations
    extended_palette <- c(
      scales::alpha(base_colors[1], 0.5), base_colors[1],
      scales::alpha(base_colors[2], 0.5), base_colors[2],
      scales::alpha(base_colors[3], 0.5), base_colors[3]
    )
    
    # Adjust names in the extended_palette to match the unique_interaction_values
    names(extended_palette) <- c(
      "pearson.CpG_permutation.negative", "pearson.CpG_permutation.positive",
      "pearson.horizontal.negative", "pearson.horizontal.positive",
      "pearson.real.negative", "pearson.real.positive"
    )
    
    # Updated ggplot code using the corrected color palette
    ggplot(df_plot_simmulations, aes(x = alpha_string, y = significant, color = interaction(name, slope))) +
      geom_boxplot(
        data = df_plot_simmulations[df_plot_simmulations$name != "pearson.real", ],
        aes(color = interaction(name, slope)),  # Ensure this interaction correctly matches palette names
        outlier.shape = NA  # Hides outliers
      ) +
      geom_point(
        data = df_plot_simmulations[df_plot_simmulations$name == "pearson.real", ],
        aes(color = interaction(name, slope)),
        size = 2,  # Adjust size for visibility
        shape = 19  # Solid circle
      ) +
      xlab(label = expression(alpha ~ " significance cut-off")) +
      ylab(label = "Number of CpG passing the significance threshold") +
      scale_color_manual(
        values = extended_palette,
        name = "Dataset Type",
        labels = names(extended_palette)
      ) +
      theme_minimal() +
      theme(legend.position = c(0.8, 0.5)) +
      ylim(c(0, 800))
    
  
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
    rbind(real_CpGs_with_delta_cutoff[real_CpGs_with_delta_cutoff$test == "pearson",],
          CpGs_with_delta_cutoff[CpGs_with_delta_cutoff$test == "pearson",])
  
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
  #   #ylab(bquote("Number of CpG with methylation variance " ~ delta[meth] >)) +
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
    #  ylab(expression(paste("Number of CpG with methylation variance ", delta[meth] >)))+
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
    ggplot(data = count_CpGs_with_delta_histogramm[count_CpGs_with_delta_histogramm$test == "pearson", ],
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
  test <- "KW"
  
  print("summerize_horizonal_columns_permutations")
  # list.files(path =  file.path(OUTPUT_FOLDER, "simulation"))
  
  # file_names <-
  #   list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
  #              pattern = "*.horizontal.*")
  
  file_names <-
    list.files(
      path = file.path(OUTPUT_FOLDER, "simulation"),
      pattern = paste(test, "horizontal.*", sep = ".")
    )
  
  start_time <- Sys.time()
  sim_results <-
    parallel_summerize_permutations(sim_file_names = file.path(OUTPUT_FOLDER, "simulation", file_names))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  # saveRDS(
  #   object = sim_results,
  #   file = file.path(OUTPUT_FOLDER, "horizontal.sim_results.rds")
  # )
  
  saveRDS(object = sim_results,
          file = file.path(
            OUTPUT_FOLDER,
            paste(test, "horizontal.sim_results.rds", sep = ".")
          ))
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
  
  
  
  #pearson.real <- readRDS("12.pipeline/real_results.rds")
  
}

# aggregate CpG results -------------------------------------------------------
if (aggregate_CpG_results) {
  
  #   test <- "pearson"
  
  real_results <-
    readRDS(file = file.path(OUTPUT_FOLDER, "real_results.rds"))
  
  test.real_results <- real_results[real_results$test == test,]
  
  #sim_results <-
  #  readRDS(file = file.path(OUTPUT_FOLDER, "horizontal.sim_results.rds"))
  # sim_results <-
  #   readRDS(file = file.path(OUTPUT_FOLDER, "CpG_permutation.sim_results.rds"))
  
  
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
  
  # df_empirical_means$real_count <-
  #   real_results$n_signfincant_CpG[order(real_results$min_delta,
  #                                        real_results$minus_log_alpha,
  #                                        real_results$test)]
  # 
  df_empirical_means$real_count <-
    test.real_results$n_signfincant_CpG[order(test.real_results$min_delta,
                                              test.real_results$minus_log_alpha,
                                              test.real_results$test)]
  
  df_empirical_means$FDR <-
    df_empirical_means$n_signfincant_CpG / df_empirical_means$real_count
  
  df_empirical_means$observed_Vs_expected <-  df_empirical_means$real_count / df_empirical_means$n_signfincant_CpG
  
  saveRDS(
    object = df_empirical_means,
    file = file.path(OUTPUT_FOLDER, paste(test,"df_CpG_vertical_empirical_means.rds",sep = "."))
  )
  #saveRDS(object = df_empirical_means,file = file.path(OUTPUT_FOLDER,"df_horizontal_column_empirical_means.rds"))
  
  ggplot(
    data = sim_results[sim_results$min_delta == 0.39, ],
    mapping = aes(x = minus_log_alpha,
                  y = n_signfincant_CpG,
                  color = name)
  ) +
    geom_point() +
    theme(legend.position = "none")
  
  
  for (delta in unique(sim_results$min_delta)) {
    print(delta)
    sims_delta <- sim_results[sim_results$min_delta == delta,]
    
    for (minus_log_alpha in seq(from = 5, to = 10, by = 1)) {
      print(head(real_results[real_results$min_delta == delta &
                                real_results$minus_log_alpha == minus_log_alpha,]))
      for (t in unique(sims_delta$test)) {
        print(t)
        sims <-
          sims_delta[sims_delta$minus_log_alpha == minus_log_alpha &
                       sims_delta$test == t,]
        sims <-
          sims[order(sims$n_signfincant_CpG, decreasing = TRUE),]
        print(head(sims, 3))
        print(nrow(sims))
      }
    }
  }
  
  
}
# aggregate horizontal results -------------------------------------------------------
if (aggregate_horizontal_results) {
  test <- "KW"
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
    ),]
  
  
  df_empirical_means <-
    df_empirical_means[df_empirical_means$test == test,]
  
  real_results_test <- real_results[test ==  real_results$test,]
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
  df_empirical_means$observed_Vs_expected <-  df_empirical_means$real_count / df_empirical_means$n_signfincant_CpG
  
  #saveRDS(object = df_empirical_means,file = file.path(OUTPUT_FOLDER,"df_CpG_vertical_empirical_means_cluster.rds"))
  saveRDS(object = df_empirical_means,
          file = file.path(
            OUTPUT_FOLDER,
            paste(test, "df_horizontal_column_empirical_means.rds", sep = ".")
          ))
}
# plot landscape -------------------------------------------------------
if (plot_landscape) {
  library(plotly)
  
  load_variable_if_not_exists(
    variable_name =  "df_peak_CpG_complete_with_test" ,
    file_path = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds")
  )
  
  min_alpha <- 6
  max_alpha <- 9.9
  min_delta <- 0.1
  test <- "KW"
  min_real_count <- 2
  
  load_variable_if_not_exists(
    variable_name =  "df_empirical_means" ,
    file_path = file.path(OUTPUT_FOLDER, paste(test,"df_horizontal_column_empirical_means.rds",sep = "."))
  )
  # df_horizontal_column_empirical_means.rds "df_CpG_vertical_empirical_means.rds"
  
  # min(best_CpGs$start)
  # [1] 29633535
  # > max(best_CpGs$start)
  # [1] 29633767
  
  search_df <-
    df_empirical_means[df_empirical_means$minus_log_alpha >= min_alpha &
                         df_empirical_means$minus_log_alpha <= max_alpha &
                         df_empirical_means$min_delta > min_delta &
                         df_empirical_means$test == test &
                         min_real_count <= df_empirical_means$real_count
                       # &  df_empirical_means$n_signfincant_CpG > 0
                       ,]
  
  plot <-plotly::plot_ly(
    data = search_df,
    x = ~min_delta,
    y = ~minus_log_alpha,
    z = ~observed_Vs_expected,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      color = ~observed_Vs_expected,
      colorscale = c('#683531', '#FFE1A1'),
      showscale = TRUE
    )
  )
  
  plot <- layout(plot, scene = list(
    xaxis = list(
      title = '\u03B4<sub>meth</sub>',
      titlefont = list(size = 30, color = "black", family = "Arial, bold")  # Specify bold font for X axis title
    ),
    yaxis = list(
      title = '-log(p<sub>value</sub>)',
      titlefont = list(size = 30, color = "black", family = "Arial, bold")  # Specify bold font for Y axis title
    ),
    zaxis = list(
      title = 'observed versus null',
      titlefont = list(size = 30, color = "black", family = "Arial, bold")  # Specify bold font for Z axis title
    )
  ))
  
  plot
  
  
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
  
  best_values <- search_df[which.min(search_df$FDR),]
  best_CpGs <-
    df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$kw.delta > best_values$min_delta &
                                     df_peak_CpG_complete_with_test$kw.p_val <= exp(-best_values$minus_log_alpha), ]
  
  folder_name <- paste("delta",
                       best_values$min_delta,
                       "minLogP",
                       best_values$minus_log_alpha,
                       sep = "_")
  
  dir.create(
    path = file.path(OUTPUT_FOLDER, "results", test, folder_name),
    showWarnings = FALSE
  )
  
  #figure 3A Figure S
  saveRDS(
    object = best_CpGs,
    file = file.path(OUTPUT_FOLDER,
                     "results",
                     test,
                     folder_name,
                     "best_CpGs.rds")
  )
  
  
  # chnaged folder t new locatrion
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = nrow(best_CpGs), style = 3)
  for (r in 1:nrow(best_CpGs)) {
    p <-
      plot_food_correlation(df_row = best_CpGs[r,], typisation = real_food_typisation)
    ggsave(
      plot = p,
      filename = file.path(OUTPUT_FOLDER_pearson_max_q,
        "plots",
        paste(best_CpGs$chrom[r],
              best_CpGs$start[r],
              "png",
              sep = ".")
      ),
      width = 6,
      height = 5
    )
    
    p <-
      plot_age_correlation_type(df_row = best_CpGs[r,], typisation = real_age_typisation)
    ggsave(
      plot = p,
      filename = file.path(OUTPUT_FOLDER_pearson_max_q,
                           "plots",
        paste("linear",
              best_CpGs$chrom[r],
              best_CpGs$start[r],
              "png",
              sep = ".")
      ),
      width = 6,
      height = 5
    )
    setTxtProgressBar(pb, r)
  }
  close(pb)
}
# plot pearson landscape --------------------------------------------------
if (pearson_landscape) {
  
  load_variable_if_not_exists(
    variable_name =  "df_peak_CpG_complete_with_test" ,
    file_path = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds")
  )
  
  min_alpha <- 6.5
  max_alpha <- 9.9
  min_delta <- 0.2
  test <- "pearson"
  min_real_count <- 1
  
  load_variable_if_not_exists(
    variable_name =  "df_empirical_means" ,
    file_path = file.path(OUTPUT_FOLDER, paste(test,"df_CpG_vertical_empirical_means.rds",sep = "."))
  )
  rm(df_empirical_means)
  gc()
  load_variable_if_not_exists(
    variable_name =  "df_empirical_means" ,
    file_path = file.path(OUTPUT_FOLDER, paste(test,"df_horizontal_column_empirical_means.rds",sep = "."))
  )
  
  
  
    search_df <-
    df_empirical_means[df_empirical_means$minus_log_alpha >= min_alpha &
                         df_empirical_means$minus_log_alpha <= max_alpha &
                         df_empirical_means$min_delta > min_delta &
                         df_empirical_means$test == test &
                         min_real_count <= df_empirical_means$real_count,]
  
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
  
  
  # index <- all_CpG_complete_with_test.45$name != "NO_CHIP" &
  #   all_CpG_complete_with_test.45$pearson.delta > 0.39 &
  #   all_CpG_complete_with_test.45$pearson.p_val <= exp(-9.3)
  #
  # best_CpGs <-
  #   all_CpG_complete_with_test.45[index,]
  
  # best_CpGs2 <-
  #   df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$pearson.delta > 0.3 &
  #                                    df_peak_CpG_complete_with_test$pearson.p_val <= exp(-9.3), ]
  
  
  best_values <- search_df[which.min(search_df$FDR), ]
  best_CpGs <-
    df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$pearson.delta > best_values$min_delta &
                                     df_peak_CpG_complete_with_test$pearson.p_val <= exp(-best_values$minus_log_alpha),]
  
  
  
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
      plot_age_correlation(df_row = best_CpGs[r,], typisation = real_age_typisation)
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
  all_CpG_complete_with_test.45 <-
    readRDS(
      "12.pipeline/results/WholeGenome/all_CpG_complete_with_test.45.rds"
    )
  
  p_significance_cut_off <- exp(-9.3)
  print(paste("p_significance_cut_off", p_significance_cut_off))
  print(paste(
    "-log(p_significance_cut_off)",
    -log(p_significance_cut_off)
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
    all_CpG_complete_with_test.45[all_CpG_complete_with_test.45$name != "NO_CHIP" &
                                 all_CpG_complete_with_test.45$pearson.delta > min_delta  &
                                 all_CpG_complete_with_test.45$pearson.p_val <= p_significance_cut_off, ]
  print(nrow(significant_CpGs))
  
  
  # Create a data frame with counts of each state
  state_counts <- table(significant_CpGs$state)
  state_counts_df <-
    data.frame(state = names(state_counts),
               count = as.numeric(state_counts))
  state_counts_df$state_color <-
    factor(x = df_universal_annotation$color,
           levels = unique(df_universal_annotation$color))
  state_counts_df$state_group <-
    factor(x =  df_universal_annotation$Group,
           levels = unique(df_universal_annotation$Group))
  state_counts_df$state <-
    factor(x = state_counts_df$state,
           levels = unique(state_counts_df$state))
  big_labels <- as.character(state_counts_df$state)
  big_labels[state_counts_df$count == 0] <-  " "
  state_counts_df$big_labels <- big_labels
  count_label <- as.character(state_counts_df$count)
  count_label[state_counts_df$count == 0] <-  " "
  # Create the pie chart with borders
  p <-
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
             color = "black",
             # Add borders to segments
             #mapping = aes(fill = state_group)
             ) +
             geomtextpath::geom_textpath(
               data = state_counts_df,
               position = position_stack(vjust = 0.5),
               aes(x = 1.75,),
               size = 3.4#,
               #show.legend = FALSE
             ) +
               geomtextpath::geom_textpath(
                 data = state_counts_df,
                 position = position_stack(vjust = 0.5),
                 aes(x = 1.2,),
                 size = 3.4,
                 label = count_label
                 #show.legend = FALSE
               ) +
               #theme(legend.position = "none")+
               scale_fill_manual(values = as.character(df_universal_annotation$color))
             
             ggsave(
               filename = file.path(
                 OUTPUT_FOLDER,
                 folder_name,
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
               temp <- significant_CpGs[significant_CpGs$name == peak, ]
               
               bed_for_meme$start[r] <- min(temp$start) - add_nucleotides
               bed_for_meme$end[r] <- max(temp$start) + add_nucleotides
             }
             
             bed_for_meme$score <- 1000.
             bed_for_meme$strand <- "."
             
             write.table(
               x = bed_for_meme[, ],
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
                         significant_CpGs$end[r] <= df_GH$end,]
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
             
             print(paste("#of positive correlation CpGs : ",   nrow(significant_CpGs[significant_CpGs$pearson.statistic > 0 ,])))
             
             
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
             
             # List to store each plot
             plot_list <- list()
             
             for (r in 1:nrow(significant_CpGs)) {
               p <-
                 plot_age_correlation(df_row = significant_CpGs[r,], typisation = real_age_typisation)
               
               # Append the plot to the list
               plot_list[[r]] <- p
               
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
                 ),width = 8,height = 8
               )
             }
}

combined_plot <- wrap_plots(plot_list[!significant_CpGs$GeneAnnotation %in% c("CIDEB","FMO5","SIRT1","NR3C1")], ncol = 9)
print(combined_plot)

ggsave(
  filename = file.path(
    OUTPUT_FOLDER,
    "results",
    test,
    folder_name,
    "plots",
    paste("rest_combined",
      ".png",
      sep = "."
    )
  ),width = 35,
  height = 20
)


# https://maayanlab.cloud/Enrichr/enrich?dataset=47ed0059768b38095dca92d0748c5055

# plot sim vs real --------------------------------------------------------
if (OTHER) {
  library(plotly)
  print("OTHER")
  plot_age_correlation(df_row = df_peak_CpG_complete_with_test[5,], typisation = real_age_typisation)
  
  plot_food_correlation(df_row = df_peak_CpG_complete_with_test[5,], typisation = real_food_typisation)
  
  search_df[which.min(search_df$FDR), ]
  
  best_ex <-
    df_CpG[!is.na(df_CpG$KW_p_val) &
             df_CpG$type_delta > 0.25 &
             df_CpG$KW_p_val < exp(-9.3), ]
  
  ex1 <- best_ex[1, ]
  ex2 <- best_ex[2, ]
  
  
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
           df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3), ]
  
  # min_delta minus_log_alpha n_signfincant_CpG real_count   FDR
  #    0.08             6.4              93.6        300 0.312
  
  length(unique(df_CpG$name[!is.na(df_CpG$KW_p_val) &
                              df_CpG$type_delta > 0.08 &
                              df_CpG$KW_p_val < exp(-6.4)]))
  
  # min_delta minus_log_alpha n_signfincant_CpG real_count   FDR
  #    0.25             9.3              0.17          2    0.085
  df_CpG[!is.na(df_CpG$KW_p_val) &
           df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3), ]
  
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
  min_delta <- 0.5
  max_q <- 0.01
  
  OUTPUT_FOLDER_pearson_max_q <- file.path(OUTPUT_FOLDER_pearson,paste0("max_q_",max_q))
  dir.create(OUTPUT_FOLDER_pearson_max_q, showWarnings = FALSE)
  
  best_CpGs <- df_peak_CpG_complete_with_test[
    df_peak_CpG_complete_with_test$pearson.delta > min_delta &
      df_peak_CpG_complete_with_test$pearson.q_vals.min_delta_0.5 < max_q ,]
  
  nrow(best_CpGs)
 
  bed4 <- best_CpGs[, c("chrom", "start", "end", "name")]
  
  write.table(
    x = bed4,
    file = file.path(OUTPUT_FOLDER_pearson_max_q,"best_CpGs.bed"),
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
    temp <- best_CpGs[best_CpGs$name == peak, ]
    
    bed_for_meme$start[r] <- min(temp$start) - add_nucleotides
    bed_for_meme$end[r] <- max(temp$start) + add_nucleotides
  }
  
  bed_for_meme$score <- 1000.
  bed_for_meme$strand <- "."
  
  write.table(
    x = bed_for_meme[, c("chrom","start","end","name")],
    file = file.path(OUTPUT_FOLDER_pearson_max_q,
      paste("bed_for_meme",
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
  
  # Split into chunks of max 999 rows
  chunk_size <- 999
  num_chunks <- ceiling(nrow(bed4) / chunk_size)
  
  for (i in seq_len(num_chunks)) {
    start_row <- (i - 1) * chunk_size + 1
    end_row <- min(i * chunk_size, nrow(bed4))
    
    bed_chunk <- bed4[start_row:end_row, ]
    
    write.table(
      x = bed_chunk,
      file = file.path(OUTPUT_FOLDER_pearson_max_q, paste0("best_CpGs_part", i, ".bed")),  # ✅ CHANGED
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  
  ### make Genhancer query with files https://genome.ucsc.edu/cgi-bin/hgTables
  # https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2502956011_yAfYHsfwhLRkg8cQImu8jfZaSnG8&clade=mammal&org=Human&db=hg19&hgta_group=regulation&hgta_track=geneHancer&hgta_table=geneHancerInteractionsDoubleElite&hgta_regionType=userRegions&position=chr9%3A1%2C332%2C600-1%2C332%2C609&hgta_outputType=bed&hgta_outFileName=
  # GH Interactio double elite hg19
  
  # Database: hg19    Primary Table: geneHancerInteractionsDoubleElite Data last updated: 2019-01-15
  # Big Bed File Download: /gbdb/hg19/geneHancer/geneHancerInteractionsDoubleElite.v2.hg19.bb
  
  # Read all hgTables_part*.BED files into one data frame
  hg_files <- list.files(
    path = OUTPUT_FOLDER_pearson_max_q,
    pattern = "^hgTables_part\\d+\\.BED$",  # Match hgTables_part1.BED, part2.BED, etc.
    full.names = TRUE
  )
  
  df_GH_list <- lapply(hg_files, function(f) {
    read.delim(
      file = f,
      header = FALSE,
      skip = 1,
      col.names = colnames(best_CpGs)[1:5]
    )
  })
  
  df_GH <- do.call(rbind, df_GH_list)
  
  
  best_CpGs$GeneAnnotation <- NA
  best_CpGs$EnhancerAnnotation <- NA
  
  for (r in 1:nrow(best_CpGs)) {
    # r <- 1
    annotated_GH <-
      df_GH[best_CpGs$chrom[r] == df_GH$chrom &
              df_GH$start <= best_CpGs$start[r] &
              best_CpGs$end[r] <= df_GH$end,]
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
    file.path(OUTPUT_FOLDER_pearson_max_q,
      "best_CpGs.xlsx"
    ),
    sheetName = "significant_CpGS",
    colNames = TRUE
  )
  
  
  write.csv2(
    x = unique(best_CpGs$GeneAnnotation),
    file = file.path(OUTPUT_FOLDER_pearson_max_q,
      paste(
        "unique.genes",
        min_delta,
        "delta.GH_Interactions",
        "q",
        max_q,
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
  normalize_PCA <- FALSE
  if (normalize_PCA) {
    normalized_string <- "normalized"
  } else{
    normalized_string <- ""
  }

  # meta37 <-
  #   meta[as.character(meta$sample) %in% as.character(real_food_typisation$sample),]
  # 
  # meta37 <-
  #   meta
  # 
  # truncate to only three groups
  # trauncate min delta pearson.delta
  data_meth <-
    df_peak_CpG_complete_with_test_min_delta[, as.character(real_age_typisation$sample)]
  #print("erase non complete rows")
  data_meth <- data_meth[complete.cases(data_meth),]
  data_meth <-
    t(data_meth)
  
  #print(sum(complete.cases(df_peak_CpG_complete_with_test)/nrow(data_meth)))
  #data_meth <- as.numeric(data_meth)
  
  data_meth_with_variance <-
    data_meth[, which(apply(data_meth, 2, var) != 0)]
  pca_prcomp_res <-
    prcomp(x = data_meth_with_variance, scale = normalize_PCA)
  
  pca_summery <- summary(pca_prcomp_res)$importance[2,]
  
  saveRDS(object = pca_prcomp_res, file =  file.path(
    OUTPUT_FOLDER,
    paste(
      "is_NORMALIZED",
      normalize_PCA,
      "pca_prcomp_res",
      "N_SAMPLES",
      "45",
      "rds",
      sep = "."
    )
  ))
  
  df_pca <- cbind(meta, pca_prcomp_res$x)
  
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
  
  # p_pca <-
  #   ggplot(data = df_pca , aes(
  #     x = PC1,
  #     y = PC2,
  #     label = sample,
  #     color = Type
  #   )) +
  #   geom_point() +
  #   ggrepel::geom_text_repel() +
  #   #labs(subtitle = "PCA_1v2 & Type") +
  #   theme_minimal() +
  #   xlab(paste("PC1 : ", round(pca_summery[1] * 100, 1), "%")) +
  #   ylab(paste("PC2 : ", round(pca_summery[2] * 100, 1), "%"))+
  #   scale_color_manual(values = setNames(mean_age_by_type$color, mean_age_by_type$Type))
  
  p_pca <-
    ggplot(data = df_pca, aes(
      x = PC1,
      y = PC2,
      label = sample,
      color = Type
    )) +
    geom_point(size = 3, alpha = 0.8) +  # Larger points with slight transparency
    ggrepel::geom_text_repel(
      size = 3.5,  # Slightly larger text size for better readability
      box.padding = unit(0.35, "lines"),  # Add padding around text
      point.padding = unit(0.5, "lines"),  # Avoid clashing text with points
      segment.color = 'grey50'  # Use a softer line color for text pointers
    ) +
    scale_color_manual(values = setNames(mean_age_by_type$color, mean_age_by_type$Type)) +
    labs(
      title = "Principal Component Analysis (PCA)",  # Add a main title
      subtitle = "PCA Scatter Plot of PC1 vs PC2 by diet",  # Include a descriptive subtitle
      caption = "Each point represents a sample colored by its diet."  # Caption at the bottom
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),  # Bold and larger title font
      plot.subtitle = element_text(size = 14),  # Slightly smaller subtitle font
      plot.caption = element_text(size = 12, margin = margin(t = 10, b = 10)),  # Caption with margins for space
      axis.title = element_text(size = 14),  # Axis titles sizing
      axis.text = element_text(size = 12),  # Axis labels sizing
      legend.position = "right",  # Ensure the legend is positioned well
      legend.title = element_blank(),  # Remove the legend title for cleaner look
      legend.text = element_text(size = 12)  # Legend text sizing
    ) +
    xlab(paste("PC1: ", round(pca_summery[1] * 100, 1), "% Variance Explained")) +  # More descriptive axis label
    ylab(paste("PC2: ", round(pca_summery[2] * 100, 1), "% Variance Explained"))  # More descriptive axis label
  
  
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


# DMP analysis 
if(FALSE){
  # criterion: methylation difference greater than a threshold
threshold_delta <- 0.3  # Adjust based on your specific criteria
threshold_p <- 0.001
df_peak_CpG_complete_with_test$DMP <-
  df_peak_CpG_complete_with_test$pearson.delta > threshold_delta &
  df_peak_CpG_complete_with_test$pearson.p_val <= threshold_p # Add a logical column for DMP status

library(dplyr)

# Assuming df_peak_CpG_complete_with_test is your data frame
df_peak_CpG_complete_with_test <- df_peak_CpG_complete_with_test %>%
  mutate(DMP = (pearson.delta > threshold_delta) & (pearson.p_val <= threshold_p))


# Sort by chromosome, start position, and name to ensure grouping works correctly
df_sorted <- df_peak_CpG_complete_with_test %>%
  arrange(chrom, start, name)

# Function to find clusters of DMPs and calculate additional metrics
find_dmp_clusters <- function(df) {
  df <- df %>%
    mutate(cluster = cumsum(!DMP | lag(name, default = first(name)) != name))
  
  grouped <- df %>%
    group_by(chrom, name, cluster) %>%
    filter(any(DMP)) %>%
    summarise(
      start = min(start),
      end = max(end),
      avg_delta = mean(pearson.delta),
      num_CpGs = n(),  # Count number of CpGs in the cluster
      segment_length = end - start + 1,  # Calculate the length of the segment
      median_p_value = median(pearson.p_val, na.rm = TRUE),  # Calculate median p-value
      .groups = 'drop'
    ) %>%
    arrange(chrom, start)
  
  # Reset cluster ids to be consecutive
  grouped$cluster <- as.integer(factor(grouped$cluster))
  
  return(grouped)
}

# Apply the function
dmr_groups <- find_dmp_clusters(df_sorted)

dmr_groups[dmr_groups$num_CpGs == max(dmr_groups$num_CpGs),]

dmr_groups[dmr_groups$segment_length == max(dmr_groups$segment_length),]
}
