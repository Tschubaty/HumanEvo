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
# library(foreach)
# library(doParallel)
library(ggplot2)
#library(Rtsne)
#library(ggrepel)
library("snow")
library("parallel")

# set constants  ----------------------------------------------------------

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

# set running parameters -------------------------------------------------------
create_true_stat <- FALSE
create_permutations <- FALSE

# load Data ---------------------------------------------------------------

df_peak_CpG <-
  readRDS(
    "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/H3K27ac.only.39.samples.merged.hg19.rds"
  )
meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
meta <- meta[order(meta$age_mean_BP), ]
real_age_typisation <-
  meta[order(meta$age_mean_BP), c("sample", "age_mean_BP")]
real_food_typisation <-
  meta[meta$Type %in% TYPES, c("sample", "Type")]


# define functions ---------------------------------------------------------

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
      as.numeric(df_row[typisation$sample[typisation$Type == ty]])
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
    meth <- as.numeric(df_row[typisation$sample[typisation$Type == s]])
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
  # Export it form base to workspace
  clusterExport(clus, deparse(substitute(age.typisation)))
  
  df_meth <-
    df[, age.typisation$sample[order(age.typisation$age_mean_BP)]]
  age <-
    age.typisation$age_mean_BP[order(age.typisation$age_mean_BP)]
  
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
      x = df[, typisation$sample],
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
create_horizontal_permutated_typisation <- function(typisation){
  permuation_order <- sample(nrow(typisation))
  sim_typisation <-  typisation
  sim_typisation$sample <-sim_typisation$sample[permuation_order]
  return(list(sim_typisation=sim_typisation,permuation_order = permuation_order))
}

# Calculate true Data -----------------------------------------------------------
if (create_true_stat) {
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
    df_peak_CpG[nrow(real_age_typisation) == df_peak_CpG$score, ]
  
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
  saveRDS(
    object = df_peak_CpG_complete_with_test,
    file = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds")
  )
} else{
  df_peak_CpG_complete_with_test <-
    readRDS(file = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds"))
}

# Calculate simulations --------------------------------------------------

if (create_permutations) {
  dir.create(path = file.path(OUTPUT_FOLDER,"simulation"), showWarnings = FALSE)
  n_repetitions <- 99 
  for (rep in 1:n_repetitions) {
    # rep <- 1
    start_time <- Sys.time()
    print(rep)

    permutation_age  <- create_horizontal_permutated_typisation(real_age_typisation)
    permutation_food  <- create_horizontal_permutated_typisation(real_food_typisation)
    
    sim_age_typisation <- permutation_age$sim_typisation
    sim_food_typisation <- permutation_food$sim_typisation
    
    permutation_age$data <- parallel_testing_pearson_cor(df = df_peak_CpG_complete_with_test,age.typisation = sim_age_typisation)
    permutation_food$data <- parallel_testing_kruskall_valis(df = df_peak_CpG_complete_with_test, food.typisation = sim_food_typisation)

    # save age permutation 
    sim_age_file_name <-
      paste("pearson",
            "horizontal",
            paste(permutation_age$permuation_order,collapse = "_"),
            "rds",
            sep = ".")
    saveRDS(object = permutation_age,
            file = file.path(OUTPUT_FOLDER, "simulation", sim_age_file_name))   
    
    # save food permutation 
    sim_food_file_name <-
      paste("KW",
            "horizontal",
            paste(permutation_food$permuation_order,collapse = "_"),
            "rds",
            sep = ".")
    saveRDS(object = permutation_food,
            file = file.path(OUTPUT_FOLDER, "simulation", sim_food_file_name))   

    end_time <- Sys.time()
    print(end_time - start_time)
  }
}
#############################################################################


##################################################################################s

############ all simmulations collected
RESOLUTION_WIDTH <- 0.25


df_simulation.stats <-
  data.frame(
    file_name = list.files(
      path = file.path(OUTPUT_FOLDER, "simulation"),
      pattern = "Kruskal_Wallis.type.label_shuffle.*"
    ),
    hemming.distance_permuation = NA,
    permuation_order = NA
  )

string_permuation_order <-
  gsub(pattern = "Kruskal_Wallis.type.label_shuffle.|.rds",
       replacement = "",
       x = df_simulation.stats$file_name)
df_simulation.stats$permuation_order <- string_permuation_order

permuation_orders <-
  lapply(
    X = string_permuation_order,
    FUN = function(string) {
      numbers <- as.numeric(unlist(strsplit(x = string, split = "\\.")))
      return(as.vector(numbers))
    }
  )

type_real_order1 <-
  as.numeric(factor(
    real_sample_typisation$type,
    levels = c("Farmer", "HG", "Steppe")
  ))
type_real_order2 <-
  as.numeric(factor(
    real_sample_typisation$type,
    levels = c("Steppe", "Farmer", "HG")
  ))
type_real_order3 <-
  as.numeric(factor(
    real_sample_typisation$type,
    levels = c("HG", "Steppe", "Farmer")
  ))

df_simulation.stats$hemming.distance_permuation <-
  sapply(
    permuation_orders,
    FUN =  function(permutation_order) {
      type_permutation_order <-
        as.numeric(factor(real_sample_typisation$type[permutation_order]))
      d1 <- sum(type_real_order1 != type_permutation_order)
      d2 <- sum(type_real_order2 != type_permutation_order)
      d3 <- sum(type_real_order3 != type_permutation_order)
      return(min(c(d1, d2, d3)))
    }
  )

# permutation_order <- permuation_orders[[1]]
# type_permutation_order <- as.numeric(factor(real_sample_typisation$type[permutation_order]))
# d1 <- sum(type_real_order1 != type_permutation_order)
# d2 <- sum(type_real_order2 != type_permutation_order)
# d3 <- sum(type_real_order3 != type_permutation_order)



# print("spaical case of sim number 8 is very significant")
# as.numeric(factor(real_sample_typisation$type[permuation_orders[[8]]]))
# as.numeric(factor(real_sample_typisation$type))
# significant_sample_typisation <- real_sample_typisation
# significant_sample_typisation$group <- as.numeric(factor(real_sample_typisation$type[permuation_orders[[8]]]))
# meta$significant_sample_typisation
#
# saveRDS(object = significant_sample_typisation$group,file = file.path(OUTPUT_FOLDER,"significant_sample_typisation.rds"))

min_deltas <- seq(0, 0.50, 0.01)
minus_log_alphas <- seq(from = 1, to = 10, by = 0.1)
sim_results_colnames <-
  c(
    "name",
    "type" ,
    "min_delta",
    "minus_log_alpha",
    "n_signfincant_CpG",
    "hemming.distance_permuation"
  )
sim_results <- data.frame()

for (n in 1:nrow(df_simulation.stats)) {
  start_time <- Sys.time()
  #n <- 1
  print(n)
  sim <-
    readRDS(file.path(
      OUTPUT_FOLDER,
      "simulation",
      df_simulation.stats$file_name[n]
    ))
  
  get_n_significant <- function(minus_log_alpha, min_delta) {
    p_expr <- paste("exp(-", minus_log_alpha, ")", sep = "")
    alpha <- eval(parse(text = p_expr))
    n_signfincant_CpG <-
      sum(sim$p_val < alpha & sim$delta > min_delta, na.rm = TRUE)
  }
  
  temp_sim <-
    data.frame(matrix(
      nrow = length(min_deltas) * length(minus_log_alphas),
      ncol = length(sim_results_colnames)
    ))
  colnames(temp_sim) <- sim_results_colnames
  
  temp_sim$name <- paste("sim", n, sep = "_")
  temp_sim$type <- "permutation"
  temp_sim$hemming.distance_permuation <-
    df_simulation.stats$hemming.distance_permuation[n]
  temp_sim$permuation_order = df_simulation.stats$permuation_order[n]
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
  temp_sim$n_signfincant_CpG <-
    mapply(FUN = get_n_significant ,  temp_sim$minus_log_alpha, temp_sim$min_delta)
  
  sim_results <- rbind(sim_results, temp_sim)
  
  end_time <- Sys.time()
  print(end_time - start_time)
}

###########################################################################

# real data
################################################################

get_n_significant <- function(minus_log_alpha, min_delta) {
  p_expr <- paste("exp(-", minus_log_alpha, ")", sep = "")
  alpha <- eval(parse(text = p_expr))
  n_signfincant_CpG <-
    sum(df_CpG$KW_p_val < alpha &
          df_CpG$type_delta > min_delta, na.rm = TRUE)
}

real_results <-
  data.frame(matrix(
    nrow = length(min_deltas) * length(minus_log_alphas),
    ncol = length(sim_results_colnames)
  ))
colnames(real_results) <- sim_results_colnames

real_results$name <- "real_data"
real_results$type <- "real"
real_results$hemming.distance_permuation <-
  df_simulation.stats$hemming.distance_permuation[n]
real_results$permuation_order = paste(1:37, collapse = ".")
real_results$minus_log_alpha <-
  rep(
    x = minus_log_alphas,
    times = length(min_deltas),
    length.out = NA,
    each = 1
  )
real_results$min_delta <-
  rep(
    x = min_deltas,
    times = 1,
    length.out = NA,
    each = length(minus_log_alphas)
  )
real_results$n_signfincant_CpG <-
  mapply(FUN = get_n_significant ,  real_results$minus_log_alpha, real_results$min_delta)


# aggirgate means
################################################################

df_empirical_means <-
  aggregate(n_signfincant_CpG ~ min_delta + minus_log_alpha,
            sim_results,
            FUN =
              mean)

df_empirical_means <-
  df_empirical_means[order(df_empirical_means$min_delta,
                           df_empirical_means$minus_log_alpha), ]

df_empirical_means$real_count <-
  real_results$n_signfincant_CpG[order(real_results$min_delta, real_results$minus_log_alpha)]

df_empirical_means$FDR <-
  df_empirical_means$n_signfincant_CpG / df_empirical_means$real_count

library(plotly)

plot_ly(
  df_empirical_means,
  x = ~ min_delta ,
  y = ~ minus_log_alpha,
  z = ~ FDR
)

search_df <-
  df_empirical_means[df_empirical_means$minus_log_alpha < 10, ]

plot_ly(
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

search_df[which.min(search_df$FDR), ]

best_ex <-
  df_CpG[!is.na(df_CpG$KW_p_val) &
           df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3), ]

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
                            df_CpG$type_delta > 0.08 & df_CpG$KW_p_val < exp(-6.4)]))

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