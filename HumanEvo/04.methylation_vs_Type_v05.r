#################################################################
##
##
##  input:
##
##
##
##  output:
##
##
##  v_05 06.11.2023
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
OUTPUT_FOLDER <- "gatherData"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
dir.create(path = file.path(OUTPUT_FOLDER, "simulation"),
           showWarnings = FALSE)
#INPUT_FOLDER <-
N_SAMPLES <- 37
ANNOTATION <- "hg19"
create_permutations <- TRUE

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

######################################## Input Data CODE ###########################

meta <-
  readRDS(
    "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/gatherData/37meta.rds"
  )

df_CpG <-
  readRDS(
    "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/gatherData/37_completecase_peak_CpG_hg19_cut_peak_annotation.rds"
  )

real_sample_typisation <-
  readRDS(
    "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/gatherData/37real_sample_typisation.rds"
  )


########################## define functions ########################
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

parallel_testing_kruskall_valis <- function(df, typisation) {
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
      FUN = function(r) {
        return(delta_mean_type2(df_row = r, typisation = typisation))
      }
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
  
  df_x <-
    data.frame(p_val = p_val,
               statistic = statistic,
               type_mean_delta = delta)
  return(df_x)
}

############################## produce sim ###################################
real_sample_typisation <- parallel_testing_kruskall_valis(df_CpG, real_sample_typisation) 

df_CpG_results <- cbind(df_CpG,real_sample_typisation)
saveRDS(object = df_CpG_results,file = file.path(OUTPUT_FOLDER,"37_completecase_peak_CpG_hg19_cut_peak_annotation_with_kruskall_valis.rds"))

if (create_permutations) {
  n_repetitions <- 100
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
    
    df_x <-
      parallel_testing_kruskall_valis(df = df_CpG, typisation = permuted_sample_typisation)
    
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
    
    saveRDS(object = permuted_sample_typisation, file.path(
      OUTPUT_FOLDER,
      "simulation",
      paste("permuted_sample_typisation",
            st,
            "rds",
            sep = ".")
    ))
    
    saveRDS(object = df_x,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
  }
}
#############################################################################














#############################################################################

#  sum(!is.na(random_p_val))/nrow(df_CpG)

sim_files <-
  list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
             pattern = "Kruskal_Wallis_chi_squared_p_val.type.label_shuffle.*")

# df <- readRDS(file = file.path(
#   OUTPUT_FOLDER,
#   paste(
#     "4b"
#     ,
#     "H3K27ac",
#     "with",
#     "pearson",
#     N_SAMPLES,
#     "samples",
#     ANNOTATION,
#     "rds",
#     sep = "."
#   )
# ))

sim_counts <- list()
x_val <- list()

#n <- 1

for (n in 1:length(sim_files)) {
  print(n)
  sim_pval <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", sim_files[n]))
  
  p <-
    ggplot(data = data.frame(sim_pval), mapping =  aes(x = -log(sim_pval))) +
    geom_histogram(binwidth = 0.25) + # beaks = seq(0,12,0.25)
    theme_minimal() +
    ggtitle("- log hist_sim_Kruskal_Wallis_chi_squared_p_val_CpGs for a CpG")

  bp <- ggplot_build(p)

  hist_data <- bp$data[[1]]

  sim_counts[[n]] <- hist_data$count
  x_val[[n]] <- hist_data$x
}
###########################################################################

sim_results <- data.frame()

for (n in 1:length(sim_counts)) {
  temp <-
    data.frame(count = sim_counts[[n]],
               x = x_val[[n]] ,
               name = paste("sim", n, sep = "_"),
               file = sim_files[n])

  sim_results <- rbind(sim_results, temp)
}

#sim_files[n]
################################################################

#ggplot(data = sim_results,mapping = aes(x = x,y = count,group=name, colour=name))+geom_point(position="dodge")
p <-
  ggplot(data = data.frame(p_val = df_CpG$KW_p_val),
         mapping =  aes(x = -log(p_val))) +
  geom_histogram(binwidth = 0.25) + # beaks = seq(0,12,0.25)
  theme_minimal() +
  ggtitle("- log hist_sim_KW_p_val_CpGs for a CpG")

bp <- ggplot_build(p)

hist_data <- bp$data[[1]]

real_results <- hist_data[, c("x", "count")]
#############################################################

real_results$label <- "real data"
real_results$file <- "realData"
real_results$name <- "original"

sim_results$label <- "label permuation"
all_results  <-
  rbind(real_results, sim_results)
#all_results$x <- as.factor(all_results$x)
saveRDS(
  object = all_results,
  file = file.path(OUTPUT_FOLDER, "p_val_dataframe_all_results.rds")
)

#############################################################

#all_results  <- rbind(real_results,horizontal_sim_results,sim_results)

uncorrected_alpha <- 0.01
min_log_p <- -log(uncorrected_alpha)

tail_results <- all_results
tail_results <- tail_results[tail_results$x >= min_log_p, ]
tail_results$x <- as.factor(tail_results$x)


p <-  ggplot(
  data = tail_results,
  mapping = aes(
    x = x,
    y = count,
    group = name,
    colour = name,
    fill = name
  )
) + geom_point(position = "jitter") +
  theme_minimal()



ggsave(
  plot = p,
  filename = file.path(
    OUTPUT_FOLDER,
    paste(
      "type permutation vs real classification",
      "png",
      sep = "."
    )
  ),
  width = 16,
  height = 10
)

##################################################################

saveRDS(object = df_CpG,file = file.path(OUTPUT_FOLDER,"peak_CpG_with_kruskal_valis_p_val.rds"))

############################################################## compute ROC ############################

all_results

mean_simulation_count <- aggregate(x = all_results[all_results$label != "real data","count"], list(all_results$x[all_results$label != "real data"]), mean)

ggplot(data = mean_simulation_count,mapping = aes(x = Group.1,y = x))+geom_line()

mean_simulation_count_CDF <- data.frame("minus.log.p" = mean_simulation_count$Group.1, empiric_CDF = NA) 
for(r in 1:nrow(mean_simulation_count)){
  
  mean_simulation_count_CDF$empiric_mean_sim_CDF[r] <-  sum(mean_simulation_count$x[mean_simulation_count$Group.1 >= mean_simulation_count_CDF$minus.log.p[r]])
  
}

real_data_count_CDF <- data.frame("minus.log.p" = real_results$x, empiric_CDF = NA) 
for(r in 1:nrow(real_results)){
  
  real_data_count_CDF$empiric_CDF[r] <-  sum(real_results$count[real_results$x >= real_data_count_CDF$minus.log.p[r]])
  
}



df_ROC <-
  data.frame(
    minus.log.p = real_data_count_CDF$minus.log.p,
    false_positive_rate =   mean_simulation_count_CDF$empiric_mean_sim_CDF[1:nrow(real_data_count_CDF)] / real_data_count_CDF$empiric_CDF,
    real_to_false_rate = real_data_count_CDF$empiric_CDF / mean_simulation_count_CDF$empiric_mean_sim_CDF[1:nrow(real_data_count_CDF)]
  )


p <- ggplot(data = df_ROC,
       mapping = aes(x = minus.log.p, y = real_to_false_rate)) + geom_line() +
  geom_point() + geom_point(data = df_ROC[which.max(df_ROC$real_to_false_rate), ], color =
                              "red")+theme_minimal()+ggtitle("optimal p value for maximum true positive rate")

ggsave(
  plot = p,
  filename = file.path(
    OUTPUT_FOLDER,
    paste(
      "true positive to false positive ratio depending on P val cutoff value",
      "png",
      sep = "."
    )
  ),
  width = 16,
  height = 10
)

##############################################################


