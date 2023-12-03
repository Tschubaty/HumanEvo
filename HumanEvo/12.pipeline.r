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
pearson_test_colnames <- c("pearson.p_val","pearson.statistic", "pearson.delta")
KW_test_colnames <- c("kw.p_val", "kw.statistic", "kw.delta")

# set running parameters -------------------------------------------------------
create_true_stat <- FALSE
create_permutations_horizonal_columns <- FALSE
summerize_CpG_columns_permutations <- FALSE
summerize_horizonal_columns_permutations <- FALSE
summerize_real_values <- TRUE
OTHER <- FALSE
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
    return(sum(df$p_val < exp(-minus_log_alpha) &
                 df$delta > min_delta, na.rm = TRUE))
  }

# parallel_compute_delta_pval_landscape <- function(sim_listfile,
#                                                   min_deltas = seq(0, 0.50, 0.01),
#                                                   minus_log_alphas = seq(from = 1, to = 10, by = 0.1)) {
#   mcmapply
# }
 
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
           minus_log_alphas = seq(from = 1, to = 10, by = 0.1)) {
    sim_results_colnames <-
      c(
        "test",
        "permuation_type" ,
        "min_delta",
        "minus_log_alpha",
        "n_signfincant_CpG",
        "name"
      )

    #create and register cluster
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(n.cores)
    clusterExport(my.cluster, "get_n_significant")
    doParallel::registerDoParallel(cl = my.cluster)
    
    sim_results <- foreach(n = 1:length(sim_file_names),
            .combine = 'rbind') %dopar% {
              
      #load data
      sim_listfile <-
        readRDS(file = sim_file_names[n])
      sim <- sim_listfile$data
      
      sim_file_name <- basename(sim_file_names[n])
      
      
      # allocated df
      temp_sim <-
        data.frame(matrix(
          nrow = length(min_deltas) * length(minus_log_alphas),
          ncol = length(sim_results_colnames)
        ))
      colnames(temp_sim) <- sim_results_colnames
      
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
     
      colnames(sim) <- c("p_val","statistic","delta")
      temp_sim$n_signfincant_CpG <-
        mapply(
          minus_log_alpha = temp_sim$minus_log_alpha,
          min_delta = temp_sim$min_delta,
          FUN = get_n_significant,
          MoreArgs = list(df=sim)
        )
      
      # append results
      return(temp_sim)
            }
    parallel::stopCluster(cl = my.cluster)
    return(sim_results)
  }



#' create CpG permuartion for each column in dataframe 
#' @param df with the  column name df$sample and
#' methylation values that will be permuted inseide 
#' @param typisation define permuted columns 
#' @returns list with permuted dataframe $data and 
#' list of permuation_order in each column in $permuation_order
create_CpG_permution <- function(df,typisation){
  
  permuation_order <- list()

  n_CpG <- nrow(df)
  for(i in 1:nrow(typisation)){
    permuation <- sample(1:n_CpG)
    df[,typisation$sample[i]] <- df[permuation,typisation$sample[i]]
    permuation_order[[typisation$sample[i]]] <- permuation
  }
  df <- df[, !colnames(df) %in% 
       c(pearson_test_colnames,KW_test_colnames)]
  return(list(data = df,permuation_order = permuation_order))
}

#' create a linear correlation plot meth vs age from row  in dataframe row 
#' @param df_row is a row in the CpG dataframe with 
#' columnnames like the typisation$sample 
#' @param typisation age typisation (can be permuted for showing sillulated values)
#' @returns plot with linear correlation line meth cs age
plot_age_correlation <- function(df_row, typisation) {
  df_plot <- typisation
  df_plot$meth <- as.numeric(df_row[, typisation$sample])
  
  p<- ggplot(data = df_plot,
             mapping = aes(x = age_mean_BP, y = meth)) + 
    ggtitle(label = paste(df_row$chrom, df_row$start),
            subtitle = paste("p_val: ",format.pval(df_row$pearson.p_val))) +
    geom_point() +
    ylim(0,1)+
    xlab("sample age estimation in Years before Present")+
    ylab("CpG methylation")+
    geom_smooth(method = "lm",)+
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
  df_plot$meth <- as.numeric(df_row[, typisation$sample])
  
  p<- ggplot(data = df_plot,
             mapping = aes(x = Type, y = meth,color = Type)) + 
    ggtitle(label = paste(df_row$chrom, df_row$start),
            subtitle = paste("p_val: ",format.pval(df_row$kw.p_val))) +
    geom_boxplot() + 
    geom_point(position = "jitter") +
    ylim(0,1)+
    xlab("sample food Type")+
    ylab("CpG methylation")+
    theme_minimal()
  return(p)
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
    file = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds"))
  
  
  saveRDS(object = list(permuation_order = NA ,sim_typisation = real_food_typisation , data = df_peak_CpG_complete_with_test[,KW_test_colnames]),
          file = file.path(OUTPUT_FOLDER, "KW.real.rds"))
  saveRDS(object = list(permuation_order = NA ,sim_typisation = real_age_typisation , data = df_peak_CpG_complete_with_test[,pearson_test_colnames]) ,
          file = file.path(OUTPUT_FOLDER, "pearson.real.rds"))
} else{
  df_peak_CpG_complete_with_test <-
    readRDS(file = file.path(OUTPUT_FOLDER, "df_peak_CpG_complete_with_test.rds"))
}

# Calculate column permutation simulations --------------------------------------------------

if (create_permutations_horizonal_columns) {
  print("create_permutations_horizonal_columns")
  dir.create(path = file.path(OUTPUT_FOLDER,"simulation"), showWarnings = FALSE)
  n_repetitions <- 100
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

# Calculate CpG permuation simulations --------------------------------------------------

if (create_CpG_permutations_vertical) {
  print("create_CpG_permutations_vertical")
  dir.create(path = file.path(OUTPUT_FOLDER,"simulation"), showWarnings = FALSE)
  n_repetitions <- 2
  for (rep in 1:n_repetitions) {
    # rep <- 1
    start_time <- Sys.time()
    print(rep)
    
    time_string <- format(start_time, "%Y_%m_%d_%H_%M_%S")
    
    permutation_CpG <- create_CpG_permution(df_peak_CpG_complete_with_test,real_age_typisation)
    
    df_permuted <- permutation_CpG$data
    
    permutation_age <- list()
    permutation_age$typisation <- real_age_typisation
    permutation_age$permuation_order <- permutation_CpG$permuation_order
    
    permutation_food <- list()
    permutation_food$typisation <- real_food_typisation
    permutation_food$permuation_order <- permutation_CpG$permuation_order
    
    permutation_age$data <- parallel_testing_pearson_cor(df = df_permuted,age.typisation = real_age_typisation)
    permutation_food$data <- parallel_testing_kruskall_valis(df = df_permuted, food.typisation = real_food_typisation)
    
    # save age permutation 
    sim_age_file_name <-
      paste("pearson",
            "CpG_permutation",
            time_string,
            "rds",
            sep = ".")
    saveRDS(object = permutation_age,
            file = file.path(OUTPUT_FOLDER, "simulation", sim_age_file_name))   
    
    # save food permutation 
    sim_food_file_name <-
      paste("KW",
            "CpG_permutation",
            time_string,
            "rds",
            sep = ".")
    saveRDS(object = permutation_food,
            file = file.path(OUTPUT_FOLDER, "simulation", sim_food_file_name))   
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
}


# Evaluate CpG simulations ----------------------------------------------------
if(summerize_CpG_columns_permutations){
  print("summerize_CpG_columns_permutations")
  # list.files(path =  file.path(OUTPUT_FOLDER, "simulation"))
  
  file_names <- list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
                           pattern = "*.CpG_permutation.*")
  start_time <- Sys.time()
  sim_results <-
    parallel_summerize_permutations(sim_file_names = file.path(OUTPUT_FOLDER, "simulation", file_names))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  saveRDS(object = sim_results,file = file.path(OUTPUT_FOLDER,"CpG_permutation.sim_results.rds"))
}
# Evaluate horizontal simulations ----------------------------------------------------
if(summerize_horizonal_columns_permutations){
  print("summerize_horizonal_columns_permutations")
  # list.files(path =  file.path(OUTPUT_FOLDER, "simulation"))

  file_names <- list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
                           pattern = "*.horizontal.*")
  start_time <- Sys.time()
  sim_results <-
    parallel_summerize_permutations(sim_file_names = file.path(OUTPUT_FOLDER, "simulation", file_names))
  end_time <- Sys.time()
  print(end_time - start_time)
 
saveRDS(object = sim_results,file = file.path(OUTPUT_FOLDER,"horizontal.sim_results.rds"))
}

# Evaluate real values ----------------------------------------------------
if(summerize_real_values){
  print("summerize_real_values")
file_names <- list.files(path = file.path(OUTPUT_FOLDER),
                         pattern = "*.real.*")
start_time <- Sys.time()
real_results <-
  parallel_summerize_permutations(sim_file_names = file.path(OUTPUT_FOLDER,file_names))
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(object = real_results,file = file.path(OUTPUT_FOLDER,"real_results.rds"))
}

# aggregate results -------------------------------------------------------

real_results <- readRDS(file = file.path(OUTPUT_FOLDER,"real_results.rds"))
#sim_results <- readRDS(file = file.path(OUTPUT_FOLDER,"horizontal.sim_results.rds"))
sim_results <- readRDS(file = file.path(OUTPUT_FOLDER,"CpG_permutation.sim_results.rds"))


df_empirical_means <-
  aggregate(n_signfincant_CpG ~ min_delta + minus_log_alpha + permuation_type + test,
            sim_results,
            FUN =
              mean)

df_empirical_means <-
  df_empirical_means[order(df_empirical_means$min_delta,
                           df_empirical_means$minus_log_alpha,
                           df_empirical_means$test), ]

df_empirical_means$real_count <-
  real_results$n_signfincant_CpG[order(real_results$min_delta, real_results$minus_log_alpha,real_results$test)]

df_empirical_means$FDR <-
  df_empirical_means$n_signfincant_CpG / df_empirical_means$real_count

saveRDS(object = df_empirical_means,file = file.path(OUTPUT_FOLDER,"df_empirical_means.rds"))
# plot kw landscape -------------------------------------------------------

min_alpha <- 5
max_alpha <- 10
min_delta <- 0.1
test <- "KW"
min_real_count <- 2

search_df <-
  df_empirical_means[df_empirical_means$minus_log_alpha >= min_alpha &
                       df_empirical_means$minus_log_alpha <= max_alpha &
                       df_empirical_means$min_delta > min_delta &
                       df_empirical_means$test == test & 
                       min_real_count <= df_empirical_means$real_count,]

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
  df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$kw.delta >= best_values$min_delta &
                                   df_peak_CpG_complete_with_test$kw.p_val <= exp(-best_values$minus_log_alpha), ]

folder_name <- paste(
  "delta",
  best_values$min_delta,
  "minLogP",
  best_values$minus_log_alpha,
  sep = "_"
)

dir.create(path = file.path(OUTPUT_FOLDER,"results","kw",folder_name), showWarnings = FALSE)


saveRDS(object = best_CpGs,
        file = file.path(
          OUTPUT_FOLDER,
          "results",
          "kw",
          folder_name,
          "best_CpGs.rds"
        ))
for(r in 1:nrow(best_CpGs)){
 p <-  plot_food_correlation(df_row = best_CpGs[r,],typisation = real_food_typisation)
 ggsave(plot = p,filename = file.path(
   OUTPUT_FOLDER,
   "results",
   "kw",
   folder_name,
   paste(
     best_CpGs$chrom[r],
     best_CpGs$start[r],
     ".png",
     sep = "."
   )))
}
# plot pearson landscape --------------------------------------------------

min_alpha <- 5
max_alpha <- 10
min_delta <- 0.2
test <- "pearson"
min_real_count <- 2

search_df <-
  df_empirical_means[df_empirical_means$minus_log_alpha >= min_alpha &
                       df_empirical_means$minus_log_alpha <= max_alpha &
                       df_empirical_means$min_delta > min_delta &
                       df_empirical_means$test == test & 
                       min_real_count <= df_empirical_means$real_count,]

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
best_CpGs <- df_peak_CpG_complete_with_test[df_peak_CpG_complete_with_test$pearson.delta >= best_values$min_delta &
                                              df_peak_CpG_complete_with_test$pearson.p_val <= exp(-best_values$minus_log_alpha),]

folder_name <- paste(
  "delta",
  best_values$min_delta,
  "minLogP",
  best_values$minus_log_alpha,
  sep = "_"
)

dir.create(path = file.path(OUTPUT_FOLDER,"results","pearson",folder_name), showWarnings = FALSE)


saveRDS(object = best_CpGs,
        file = file.path(
          OUTPUT_FOLDER,
          "results",
          "pearson",
          folder_name,
          "best_CpGs.rds"
        ))
for(r in 1:nrow(best_CpGs)){
  p <-  plot_age_correlation(df_row = best_CpGs[r,],typisation = real_age_typisation)
  ggsave(plot = p,filename = file.path(
    OUTPUT_FOLDER,
    "results",
    "pearson",
    folder_name,
    paste(
      best_CpGs$chrom[r],
      best_CpGs$start[r],
      ".png",
      sep = "."
    )))
}

# plot sim vs real --------------------------------------------------------


plot_age_correlation(df_row = df_peak_CpG_complete_with_test[5,],typisation = real_age_typisation)

plot_food_correlation(df_row = df_peak_CpG_complete_with_test[5,],typisation = real_food_typisation)


  print("OTHER")
if(OTHER){
library(plotly)


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
}