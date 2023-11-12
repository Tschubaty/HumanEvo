#################################################################
##
##
##  input:
##
##
##  output:
##
##
##  v_02 12.11.2023
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
real_sample_typisation <- readRDS(file.path(INPUT_FOLDER,"real_sample_typisation.rds"))



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

type_real_order1 <- as.numeric(factor(real_sample_typisation$type,levels = c("Farmer","HG","Steppe")))
type_real_order2 <- as.numeric(factor(real_sample_typisation$type,levels = c("Steppe","Farmer","HG")))
type_real_order3 <- as.numeric(factor(real_sample_typisation$type,levels = c("HG","Steppe","Farmer")))

df_simulation.stats$hemming.distance_permuation <- sapply(permuation_orders,FUN =  function(permutation_order){
  type_permutation_order <- as.numeric(factor(real_sample_typisation$type[permutation_order]))
  d1 <- sum(type_real_order1 != type_permutation_order)
  d2 <- sum(type_real_order2 != type_permutation_order)
  d3 <- sum(type_real_order3 != type_permutation_order)
  return(min(c(d1,d2,d3)))
})

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
sim_results_colnames <- c("name","type" ,"min_delta", "minus_log_alpha", "n_signfincant_CpG", "hemming.distance_permuation")
sim_results <- data.frame()

for (n in 1:nrow(df_simulation.stats)) {
  start_time <- Sys.time()
  #n <- 1
  print(n)
  sim <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", df_simulation.stats$file_name[n]))
  
  get_n_significant <- function(minus_log_alpha,min_delta) {
    p_expr <- paste("exp(-", minus_log_alpha, ")", sep = "")
    alpha <- eval(parse(text = p_expr))
    n_signfincant_CpG <-
      sum(sim$p_val < alpha & sim$delta > min_delta, na.rm = TRUE)
  }
  
  temp_sim <- data.frame(matrix(nrow = length(min_deltas)*length(minus_log_alphas),ncol = length(sim_results_colnames)))
  colnames(temp_sim) <- sim_results_colnames
  
  temp_sim$name <- paste("sim", n, sep = "_")
  temp_sim$type <- "permutation"
  temp_sim$hemming.distance_permuation <- df_simulation.stats$hemming.distance_permuation[n]
  temp_sim$permuation_order = df_simulation.stats$permuation_order[n]
  temp_sim$minus_log_alpha <- rep(x = minus_log_alphas, times = length(min_deltas), length.out = NA, each = 1) 
  temp_sim$min_delta <- rep(x = min_deltas, times = 1, length.out = NA, each = length(minus_log_alphas)) 
  temp_sim$n_signfincant_CpG <-mapply(FUN=get_n_significant ,  temp_sim$minus_log_alpha, temp_sim$min_delta)
  
  sim_results <- rbind(sim_results,temp_sim)
  
  end_time <- Sys.time()
  print(end_time - start_time)
}

###########################################################################

# real data 
################################################################

get_n_significant <- function(minus_log_alpha,min_delta) {
  p_expr <- paste("exp(-", minus_log_alpha, ")", sep = "")
  alpha <- eval(parse(text = p_expr))
  n_signfincant_CpG <-
    sum(df_CpG$KW_p_val < alpha & df_CpG$type_delta > min_delta, na.rm = TRUE)
}

real_results <- data.frame(matrix(nrow = length(min_deltas)*length(minus_log_alphas),ncol = length(sim_results_colnames)))
colnames(real_results) <- sim_results_colnames

real_results$name <- "real_data"
real_results$type <- "real"
real_results$hemming.distance_permuation <- df_simulation.stats$hemming.distance_permuation[n]
real_results$permuation_order = paste(1:37,collapse = ".")
real_results$minus_log_alpha <- rep(x = minus_log_alphas, times = length(min_deltas), length.out = NA, each = 1) 
real_results$min_delta <- rep(x = min_deltas, times = 1, length.out = NA, each = length(minus_log_alphas)) 
real_results$n_signfincant_CpG <-mapply(FUN=get_n_significant ,  real_results$minus_log_alpha, real_results$min_delta)


# aggirgate means 
################################################################

df_empirical_means <-
  aggregate(n_signfincant_CpG ~ min_delta + minus_log_alpha, sim_results, FUN =
              mean)

df_empirical_means <- df_empirical_means[order(df_empirical_means$min_delta,df_empirical_means$minus_log_alpha),]

df_empirical_means$real_count <- real_results$n_signfincant_CpG[order(real_results$min_delta,real_results$minus_log_alpha)]

df_empirical_means$FDR <- df_empirical_means$n_signfincant_CpG / df_empirical_means$real_count

library(plotly)

plot_ly(df_empirical_means, x = ~min_delta , y = ~minus_log_alpha, z = ~FDR )

search_df <- df_empirical_means[df_empirical_means$minus_log_alpha < 10,]

plot_ly(search_df , x = ~min_delta , y = ~minus_log_alpha, z = ~FDR ,marker = list(color = ~FDR, colorscale = c('#683531','#FFE1A1'), showscale = TRUE))

search_df[which.min(search_df$FDR),]

best_ex <- df_CpG[!is.na(df_CpG$KW_p_val) & df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3),]

ex1 <- best_ex[1,]
ex2 <-best_ex[2,]


df_plot <- data.frame(matrix(nrow = N_SAMPLES,ncol = 0))
df_plot$sample <- as.character(real_sample_typisation$sample)
df_plot$meth <- as.numeric(ex1[real_sample_typisation$sample])
df_plot$type <- real_sample_typisation$type

ggplot(data = df_plot,
       mapping = aes(x = type, y = meth, color = type)) + ggtitle(paste(ex1$chrom, ex1$start)) +
  geom_boxplot() + geom_point(position = "jitter")+
  theme_minimal()


df_plot <- data.frame(matrix(nrow = N_SAMPLES,ncol = 0))
df_plot$sample <- as.character(real_sample_typisation$sample)
df_plot$meth <- as.numeric(ex2[real_sample_typisation$sample])
df_plot$type <- real_sample_typisation$type

ggplot(data = df_plot,
       mapping = aes(x = type, y = meth, color = type)) + ggtitle(paste(ex2$chrom, ex2$start)) +
  geom_boxplot() + geom_point(position = "jitter")+
  theme_minimal()



df_CpG[!is.na(df_CpG$KW_p_val) & df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3),]

# min_delta minus_log_alpha n_signfincant_CpG real_count   FDR
#    0.08             6.4              93.6        300 0.312

length(unique(df_CpG$name[!is.na(df_CpG$KW_p_val) & df_CpG$type_delta > 0.08 & df_CpG$KW_p_val < exp(-6.4)]))

# min_delta minus_log_alpha n_signfincant_CpG real_count   FDR
#    0.25             9.3              0.17          2    0.085
df_CpG[!is.na(df_CpG$KW_p_val) & df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3),]

saveRDS(object = df_empirical_means,file = file.path(OUTPUT_FOLDER,"simulation","empirical_means.rds"))

# p <- ggplot(data = df_empirical_means,
#             mapping = aes(x = min_log_p, y = real_count / count)) + geom_point() + theme_minimal() +
#   ggtitle(min_delta, "min_delta")