#################################################################
##
##
##  input:
##
##
##  output:
##
##
##  v_01 13.11.2023
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
OUTPUT_FOLDER <- "10.optimal_p_val_and_delta_AGE"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
INPUT_FOLDER <-"04.methylation_vs_age_CpG"
N_SAMPLES <- 39
ANNOTATION <- "hg19"
recalculate_summery <- FALSE

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

df_CpG <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/04.methylation_vs_age_CpG/H3K27ac.with.pearson.39.samples.hg19.rds")
dim(df_CpG)
meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
#real_sample_typisation <- readRDS(file.path(INPUT_FOLDER,"real_sample_typisation.rds"))

##################################################################################s
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

delta_max <- function(df_row) {
  return(max(df_row) - min(df_row))
}

meth <- df_CpG[,meta$sample[order(meta$age_mean_BP)]]

deltas <- apply(X = meth,MARGIN = 1,FUN = delta_max)

df_CpG$delta <- as.numeric(deltas)
############ all simmulations collected

df_simulation.stats <- data.frame(file_name = list.files(path = file.path(INPUT_FOLDER, "simulation"),
                                                      pattern = "pearson_p_val.methylation.vertical_column_shuffle.*"),hemming.distance_permuation = NA,permuation_order = NA)

string_permuation_order <- gsub(pattern = "pearson_p_val.methylation.vertical_column_shuffle.|.rds",replacement = "",x = df_simulation.stats$file_name)
df_simulation.stats$permuation_order <- string_permuation_order

permuation_orders <- lapply(X = string_permuation_order,FUN = function(string){
  numbers <- as.numeric(unlist(strsplit(x = string, split = "_")))
  return(as.vector(numbers ))
})

real_order <- seq(length(permuation_orders[[1]]))

#distance hemming
df_simulation.stats$hemming.distance_permuation <-
  sapply(
    permuation_orders,
    FUN =  function(permutation_order) {
      return(sum(real_order != permutation_order))
    }
  )

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


if(recalculate_summery){
  


sim_results_colnames <- c("name","type" ,"min_delta", "minus_log_alpha", "n_signfincant_CpG", "hemming.distance_permuation")
sim_results <- data.frame()

for (n in 1:nrow(df_simulation.stats)) {
  start_time <- Sys.time()
  #n <- 1
  print(n)
  sim <-
  readRDS(file.path(INPUT_FOLDER, "simulation", df_simulation.stats$file_name[n]))
  
  get_n_significant <- function(minus_log_alpha,min_delta) {
    p_expr <- paste("exp(-", minus_log_alpha, ")", sep = "")
    alpha <- eval(parse(text = p_expr))
    n_signfincant_CpG <-
      sum(sim < alpha & df_CpG$delta > min_delta, na.rm = TRUE)
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


saveRDS(object = sim_results,file = file.path(OUTPUT_FOLDER,"sim_results.rds"))
###########################################################################

# real data 
################################################################

#sim_results <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/10.optimal_p_val_and_delta_AGE/sim_results.rds")

get_n_significant_real <- function(minus_log_alpha,min_delta) {
  p_expr <- paste("exp(-", minus_log_alpha, ")", sep = "")
  alpha <- eval(parse(text = p_expr))
  n_signfincant_CpG <-
    sum(df_CpG$pearson_p_val < alpha & df_CpG$delta > min_delta, na.rm = TRUE)
}

real_results <- data.frame(matrix(nrow = length(min_deltas)*length(minus_log_alphas),ncol = length(sim_results_colnames)))
colnames(real_results) <- sim_results_colnames

real_results$name <- "real_data"
real_results$type <- "real"
real_results$hemming.distance_permuation <- 0
real_results$permuation_order = paste(1:37,collapse = ".")
real_results$minus_log_alpha <- rep(x = minus_log_alphas, times = length(min_deltas), length.out = NA, each = 1) 
real_results$min_delta <- rep(x = min_deltas, times = 1, length.out = NA, each = length(minus_log_alphas)) 
real_results$n_signfincant_CpG <-mapply(FUN=get_n_significant_real ,  real_results$minus_log_alpha, real_results$min_delta)


# aggirgate means 
################################################################

df_empirical_means <-
  aggregate(n_signfincant_CpG ~ min_delta + minus_log_alpha, sim_results, FUN =
              mean)

df_empirical_means <- df_empirical_means[order(df_empirical_means$min_delta,df_empirical_means$minus_log_alpha),]

df_empirical_means$real_count <- real_results$n_signfincant_CpG[order(real_results$min_delta,real_results$minus_log_alpha)]

df_empirical_means$FDR <- df_empirical_means$n_signfincant_CpG / df_empirical_means$real_count

saveRDS(object = df_empirical_means,file = file.path(OUTPUT_FOLDER,"df_empirical_means.rds"))

}else{
  sim_results <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/10.optimal_p_val_and_delta_AGE/sim_results.rds")
  df_empirical_means <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/10.optimal_p_val_and_delta_AGE/df_empirical_means.rds")
}

library(plotly)

plot_ly(df_empirical_means, x = ~min_delta , y = ~minus_log_alpha, z = ~FDR )

search_df <- df_empirical_means[df_empirical_means$minus_log_alpha > min(minus_log_alphas) & df_empirical_means$min_delta > min(min_deltas),]

################################################################################################################
############################## define which minimum to take in which area ######################################
################################################################################################################
search_df <- df_empirical_means[df_empirical_means$minus_log_alpha < 8.5  & df_empirical_means$minus_log_alpha > 6 & df_empirical_means$min_delta > 0.1,]

plot_ly(search_df , x = ~min_delta , y = ~minus_log_alpha, z = ~FDR ,marker = list(color = ~FDR, colorscale = c('#683531','#FFE1A1'), showscale = TRUE))

optimal_value <- search_df[which.min(search_df$FDR),]

bestCpG <- df_CpG[!is.na(df_CpG$pearson_p_val) & df_CpG$delta > optimal_value$min_delta & df_CpG$pearson_p_val < exp(-optimal_value$minus_log_alpha),]

save.folder <- paste("min_delta",optimal_value$min_delta,"max_minus_log_p",optimal_value$minus_log_alpha,sep = "_")
dir.create(file.path(OUTPUT_FOLDER,save.folder), showWarnings = FALSE)

saveRDS(object = bestCpG,file = file.path(OUTPUT_FOLDER,save.folder,"bestCpG.rds"))

for(r in 1:nrow(bestCpG)){
  current_ex <- bestCpG[r,]
  
  df_plot <- data.frame(matrix(nrow = N_SAMPLES,ncol = 0))
  df_plot$sample <- as.character(meta$sample[order(meta$age_mean_BP)])
  df_plot$meth <- as.numeric(current_ex[meta$sample[order(meta$age_mean_BP)]])
  df_plot$age <- as.numeric(meta$age_mean_BP[order(meta$age_mean_BP)])
  
  saveRDS(object = df_plot,file = file.path(OUTPUT_FOLDER,save.folder,paste(current_ex$chrom, current_ex$start,"rds",sep = ".")))
  fit <- lm(df_plot$meth ~ df_plot$age)
  
  pval <- lmp(fit)
  
  
  p <- ggplot(data = df_plot,
         mapping = aes(x = age, y = meth)) + 
    ggtitle(paste(current_ex$chrom, current_ex$start),subtitle = paste("p-val:",format.pval(pval))) +
    geom_point()+
    # stat_smooth(method = "lm", 
    #             formula = df_plot$meth ~ df_plot$age, 
    #             geom = "smooth")+ 
    ylim(0,1)+
    theme_minimal()+
    geom_smooth(method = "lm")
  
 ggsave(plot = p,filename = file.path(OUTPUT_FOLDER,save.folder,paste(current_ex$chrom, current_ex$start,"png",sep = "."))) 
  # df_plot <- data.frame(matrix(nrow = N_SAMPLES,ncol = 0))
  # df_plot$sample <- as.character(real_sample_typisation$sample)
  # df_plot$meth <- as.numeric(ex2[real_sample_typisation$sample])
  # df_plot$type <- real_sample_typisation$type
  # 
  # ggplot(data = df_plot,
  #        mapping = aes(x = type, y = meth, color = type)) + ggtitle(paste(ex2$chrom, ex2$start)) +
  #   geom_boxplot() + geom_point(position = "jitter")+
  #   theme_minimal()  
  
}


# 
# 
# 
# df_CpG[!is.na(df_CpG$KW_p_val) & df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3),]
# 
# # min_delta minus_log_alpha n_signfincant_CpG real_count   FDR
# #    0.08             6.4              93.6        300 0.312
# 
# length(unique(df_CpG$name[!is.na(df_CpG$KW_p_val) & df_CpG$type_delta > 0.08 & df_CpG$KW_p_val < exp(-6.4)]))
# 
# # min_delta minus_log_alpha n_signfincant_CpG real_count   FDR
# #    0.25             9.3              0.17          2    0.085
# df_CpG[!is.na(df_CpG$KW_p_val) & df_CpG$type_delta > 0.25 & df_CpG$KW_p_val < exp(-9.3),]
# 
# saveRDS(object = df_empirical_means,file = file.path(OUTPUT_FOLDER,"simulation","empirical_means.rds"))
# 
# # p <- ggplot(data = df_empirical_means,
# #             mapping = aes(x = min_log_p, y = real_count / count)) + geom_point() + theme_minimal() +
# #   ggtitle(min_delta, "min_delta")