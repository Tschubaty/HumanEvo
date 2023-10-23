#################################################################
##
##
##  input:
##
##
##  output:
##
##
##  v_01 08.09.2023
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
# library(foreach)
# library(doParallel)
library(ggplot2)
#library(Rtsne)
#library(ggrepel)
library("parallel")
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- "10.optimal_p_val_and_delta"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
INPUT_FOLDER <-"04.methylation_vs_Type"
N_SAMPLES <- 37
ANNOTATION <- "hg19"
create_permutations <- FALSE

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

df_CpG <- readRDS(file =  file.path(INPUT_FOLDER,"peak_CpG_with_kruskal_valis_p_val.rds"))
meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))

TYPES <- unique(meta$Type)
## only big sample groups
TYPES <- c("Farmer", "Steppe", "HG")

sample_types <- list()
for (t in TYPES) {
  print(t)
  sample_types[[t]] <- meta$sample[meta$Type %in% t]
}

real_sample_typisation <- readRDS(file.path(INPUT_FOLDER,"real_sample_typisation.rds"))

delta_mean_type <- function(df_row) {
  df_analyis <- data.frame()
  for (s in names(sample_types)) {
    #  s <- names(sample_types)[1]
    meth <-
      as.numeric(df_row[unlist(sample_types[s], use.names = FALSE)])
    # test prerequisite if at least 5 samples
    if (sum(!is.na(meth)) < 1) {
      return(NA)
    }
    temp <- data.frame(meth = meth , type = s)
    df_analyis <- rbind(df_analyis, temp)
  }
  type_means <-
    aggregate(
      df_analyis$meth,
      list(df_analyis$type),
      FUN = function(x) {
        mean(x, na.rm = TRUE)
      }
    )
  return(max(type_means$x) - min(type_means$x))
}

delta_mean_type2 <- function(df_row, typisation) {
  typisation$meth <- NA
  for (ty in unique(typisation$type)) {
    # ty <- unique(typisation$type)[1]
    meth <-
      as.numeric(df_row[typisation$sample[typisation$type == ty]])
    # test prerequisite if at least 1 sample
    if (sum(!is.na(meth)) < 1) {
      return(NA)
    }
    typisation$meth[typisation$type == ty] <- meth
  }
  type_means <-
    aggregate(
      typisation$meth,
      list(typisation$type),
      FUN = function(x) {
        mean(x, na.rm = TRUE)
      }
    )
  return(max(type_means$x) - min(type_means$x))
}

#delta_mean_type2(df_row = df_CpG[1000,],typisation = real_sample_typisation)

############################## simulation paralle R ######################### horizontal ###

# test function for horizontal column permuation
meth_vs_type_test2 <- function(df_row, typisation) {
  typisation$meth <- NA
  for (s in unique(typisation$type)) {
    meth <- as.numeric(df_row[typisation$sample[typisation$type == s]])
    # test prerequisite if at least 5 samples
    if (sum(!is.na(meth)) < 5) {
      return(NA)
    }
    typisation$meth[typisation$type == s] <- meth
  }
  
  kru_test <- kruskal.test(meth ~ type, data = typisation)
  return(kru_test)
}

parallel_testing_kruskall_valis <- function(df,typisation){
  print("clust")
  clust <- makeCluster(parallel::detectCores() - 1)
  # Export it form base to workspace
  clusterExport(clust, "meth_vs_type_test2")
  #print(deparse(substitute(df)))
  #print(deparse(substitute(typisation)))
  clusterExport(clust, "delta_mean_type2")
  clusterExport(clust, deparse(substitute(typisation)))
  print("test")
  simulated_test <-
    parApply(
      cl = clust,
      X = df,
      MARGIN = 1,
      FUN = function(r) {
        return(meth_vs_type_test2(r, typisation))
      }
    )
  print("delta")
  delta <-
    parApply(
      cl = clust,
      X = df,
      MARGIN = 1,
      FUN = function(r) {return(delta_mean_type2(df_row = r,typisation = typisation))}
    )
  
  print("pval")
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
  print("stat")
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
  
  df_x <- data.frame(p_val = p_val,statistic = statistic,delta = delta)
  return(df_x)
}

# test_reults <- parallel_testing_kruskall_valis(df = df_CpG,typisation = real_sample_typisation)
# 
# df_CpG$KW_p_val <- test_reults[[1]]
# df_CpG$KW_statistic <- test_reults[[2]]
# df_CpG$type_delta <- test_reults[[3]]

#te <- meth_vs_type_test2(df_CpG[1000,],real_sample_typisation)
#df_CpG$p_val_Kruskal_Wallis <- meth_vs_type_test2(df_CpG,real_sample_typisation)

############################## produce sim ###################################
if (create_permutations) {
  n_repetitions <- 99
  library("snow")
  for (rep in 1:n_repetitions) {
    #
    # rep <- 1
    print(rep)
    start_time <- Sys.time()
    
    permuted_sample_typisation <- real_sample_typisation
    permuation_order <- sample(x = 1:37)
        permuted_sample_typisation$sample <-
      permuted_sample_typisation$sample[permuation_order]
    
    df_x <- parallel_testing_kruskall_valis(df = df_CpG,typisation = permuted_sample_typisation)
    
    end_time <- Sys.time()
    st = paste(permuation_order,  collapse = ".")
    print(end_time - start_time)
    
    file_name <-
      paste("Kruskal_Wallis",
            "type",
            "label_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = df_x,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
  }
}
#############################################################################

















##################################################################################s

############ all simmulations collected
RESOLUTION_WIDTH <- 0.25 


df_simulation.stats <- data.frame(file_name = list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
                                                      pattern = "Kruskal_Wallis.type.label_shuffle.*"),hemming.distance_permuation = NA,permuation_order = NA)

string_permuation_order <- gsub(pattern = "Kruskal_Wallis.type.label_shuffle.|.rds",replacement = "",x = df_simulation.stats$file_name)
df_simulation.stats$permuation_order <- string_permuation_order

permuation_orders <- lapply(X = string_permuation_order,FUN = function(string){
  numbers <- as.numeric(unlist(strsplit(x = string, split = "\\.")))
  return(as.vector(numbers ))
})

df_simulation.stats$hemming.distance_permuation <- sapply(permuation_orders,FUN =  function(permutation_order){
  type_permutation_order <- as.numeric(factor(real_sample_typisation$type[permutation_order]))
  type_real_order <- as.numeric(factor(real_sample_typisation$type))
  return(sum(type_real_order != type_permutation_order))
})

print("spaical case of sim number 8 is very significant")
as.numeric(factor(real_sample_typisation$type[permuation_orders[[8]]]))
as.numeric(factor(real_sample_typisation$type))
significant_sample_typisation <- real_sample_typisation
significant_sample_typisation$group <- as.numeric(factor(real_sample_typisation$type[permuation_orders[[8]]]))
meta$significant_sample_typisation 

saveRDS(object = significant_sample_typisation$group,file = file.path(OUTPUT_FOLDER,"significant_sample_typisation.rds"))


sim_results <- data.frame()

#n <- 1
min_delta <- 0

for (n in 1:nrow(df_simulation.stats)) {
  print(n)
  sim <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", df_simulation.stats$file_name[n]))
  
  for(min_delta in seq(0,0.25,0.01)){
    
  adjusted_df <- data.frame(sim$p_val,sim$delta,sim$statistic)
  # put p-val with low delta to 1
  adjusted_df$sim.p_val[adjusted_df$sim.delta < min_delta] <- 1
  
  p <-
    ggplot(data = adjusted_df, mapping =  aes(x = -log(sim.p_val))) +
    geom_histogram(binwidth = RESOLUTION_WIDTH) + # beaks = seq(0,12,0.25)
    theme_minimal() +
    ggtitle("- log hist_sim_Kruskal_Wallis_chi_squared_p_val_CpGs for a CpG")

  bp <- ggplot_build(p)

  hist_data <- bp$data[[1]]
  
  
  temp <- data.frame(count = hist_data$count,
             x = hist_data$x ,
             name = paste("sim", n, sep = "_"),
             type = "permutation",
             delta = min_delta)

 
  sim_results <- rbind(sim_results, temp)
  }
}
###########################################################################

# real data 
################################################################
real_results <- data.frame()
for(min_delta in seq(0,0.25,0.01)){
  
  adjusted_df <- data.frame(p_val = df_CpG$KW_p_val,delta = df_CpG$type_delta,statistic = df_CpG$KW_statistic)
  # put p-val with low delta to 1s
  adjusted_df$p_val[df_CpG$type_delta < min_delta] <- 1
  
  
p <-
  ggplot(data = adjusted_df,
         mapping =  aes(x = -log(p_val))) +
  geom_histogram(binwidth = RESOLUTION_WIDTH) + 
  theme_minimal() +
  ggtitle("- log hist_sim_KW_p_val_CpGs for a CpG")

bp <- ggplot_build(p)

hist_data <- bp$data[[1]]

temp <- data.frame(count = hist_data$count,
                   x = hist_data$x ,
                   name = "original",
                   type = "original",
                   delta = min_delta)


real_results <- rbind(real_results, temp)
}
#############################################################

all_results  <-
  rbind(real_results, sim_results)

saveRDS(
  object = all_results,
  file = file.path(OUTPUT_FOLDER, paste("empirical_p_val_resulition_width_",RESOLUTION_WIDTH,"_all_results.rds",sep = ""))
)

#############################################################


for(min_delta in seq(0, 0.25, 0.01)) {
  max_pval <- 0.01
  min_log_p <- -log(max_pval)
  
  tail_results <- all_results[all_results$delta == min_delta, ]
  tail_results <- tail_results[tail_results$x >= min_log_p,]
  tail_results$x <- as.factor(tail_results$x)
  
  
  p <-  ggplot(
    data = tail_results,
    mapping = aes(
      x = x,
      y = count,
      group = name,
      colour = type,
      fill = type,
      shape = type,
    )
  ) + geom_point(position = "jitter") +
    theme_minimal()
  
  
  
  ggsave(
    plot = p,
    filename = file.path(
      OUTPUT_FOLDER,
      paste(
        "type permutation vs real classification",
        "max_pval",
        max_pval,
        "min_delta",
        min_delta,
        "png",
        sep = "."
      )
    ),
    width = 16,
    height = 10
  )
  
  print(min_delta)
  for (x in unique(tail_results$x)) {
    bin <- tail_results[tail_results$x == 5,]
    print(bin$name[which.max(bin$count)])
  }
}

##############################################################

############################################################## compute ROC ############################

best_values <- data.frame()

for(min_delta in seq(0, 0.25, 0.01)) {
  print(min_delta)
  
  mean_simulation_count <-
    aggregate(x = all_results[all_results$delta == min_delta &
                                all_results$type != "original", "count"], list(all_results$x[all_results$delta == min_delta &
                                                                                               all_results$type != "original"]), mean)
  
  ggplot(data = mean_simulation_count, mapping = aes(x = Group.1, y = x)) +
    geom_line()
  
  mean_simulation_count_CDF <-
    data.frame("minus.log.p" = mean_simulation_count$Group.1, empiric_CDF = NA)
  for (r in 1:nrow(mean_simulation_count)) {
    mean_simulation_count_CDF$empiric_mean_sim_CDF[r] <-
      sum(mean_simulation_count$x[mean_simulation_count$Group.1 >= mean_simulation_count_CDF$minus.log.p[r]])
    
  }
  
  real_data_count_CDF <-
    data.frame(
      "minus.log.p" = real_results$x[real_results$delta == min_delta],
      empiric_CDF = NA,
      min_delta = min_delta
    )
  for (r in 1:nrow(real_data_count_CDF)) {
    df_with_delta <-  real_results[real_results$delta == min_delta, ]
    
    real_data_count_CDF$empiric_CDF[r] <-
      sum(df_with_delta$count[df_with_delta$x >= real_data_count_CDF$minus.log.p[r]])
    
  }
  
  df_ROC <-
    data.frame(
      minus.log.p = real_data_count_CDF$minus.log.p,
      false_positive_rate =   mean_simulation_count_CDF$empiric_mean_sim_CDF[1:nrow(real_data_count_CDF)] / real_data_count_CDF$empiric_CDF,
      real_to_random_rate = real_data_count_CDF$empiric_CDF / mean_simulation_count_CDF$empiric_mean_sim_CDF[1:nrow(real_data_count_CDF)]
    )
  
  
  p <- ggplot(data = df_ROC,
              mapping = aes(x = minus.log.p, y = real_to_random_rate)) + geom_line() +
    geom_point() + geom_point(data = df_ROC[which.max(df_ROC$real_to_random_rate),], color =
                                "red") + theme_minimal() + ggtitle("optimal p value for maximum true positive rate")
  
  ggsave(
    plot = p,
    filename = file.path(
      OUTPUT_FOLDER,
      paste(
        "data to random positive ratio P val cutoff value & min_delta",
        min_delta,
        "png",
        sep = "."
      )
    ),
    width = 16,
    height = 10
  )
  
  best_val <- df_ROC[which.max(df_ROC$real_to_random_rate),]
  best_val$p_val <- exp(-best_val$minus.log.p)
  best_val$min_delta <- min_delta
  best_values <-rbind(best_values,best_val)
}



###########################################################################################################

