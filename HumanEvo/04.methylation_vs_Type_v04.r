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
OUTPUT_FOLDER <- "04.methylation_vs_Type"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
#INPUT_FOLDER <-
N_SAMPLES <- 39
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


meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
sample_order <- as.character(meta$sample[order(meta$age_mean_BP)])
#age <- meta$age_mean_BP[order(meta$age_mean_BP)]

df_CpG <-
  readRDS(
    "C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/07.methylation_p_val/single_CpG.H3K27ac.hg19.39.samples.p_val.delta_meth.rds"
  )

ggplot(data = df_CpG)

TYPES <- unique(meta$Type)
## only big sample groups
TYPES <- c("Farmer", "Steppe", "HG")

sample_types <- list()
for (t in TYPES) {
  print(t)
  sample_types[[t]] <- meta$sample[meta$Type %in% t]
}


real_sample_typisation <- data.frame()
for (t in TYPES) {
  print(t)
  real_sample_typisation <-
    rbind(real_sample_typisation,
          data.frame(sample = sample_types[[t]] , type = t))
}

saveRDS(object = real_sample_typisation,file = file.path(OUTPUT_FOLDER,"real_sample_typisation.rds"))

CpG_per_type <- data.frame()


for (i in 1:length(sample_types)) {
  print(i)
  
  
  
  for (s in sample_types[[i]]) {
    print(s)
    CpG_per_type <-
      rbind(CpG_per_type , data.frame(meth = df_CpG[, s], type = names(sample_types[i])))
  }
}


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
  return(list(p_val,statistic,delta))
}

test_reults <- parallel_testing_kruskall_valis(df = df_CpG,typisation = real_sample_typisation)

df_CpG$KW_p_val <- test_reults[[1]]
df_CpG$KW_statistic <- test_reults[[2]]
df_CpG$type_delta <- test_reults[[3]]

#te <- meth_vs_type_test2(df_CpG[1000,],real_sample_typisation)
#df_CpG$p_val_Kruskal_Wallis <- meth_vs_type_test2(df_CpG,real_sample_typisation)

############################## produce sim ###################################

if (create_permutations) {
  n_repetitions <- 98
  library("snow")
  #Create cluster
  clust <- makeCluster(parallel::detectCores() - 1)
  # Export it form base to workspace
  clusterExport(clust, "meth_vs_type_test2")
  clusterExport(clust, "sample_types")
  #Apply the declared function
  for (rep in 1:n_repetitions) {
    #
    # rep <- 1
    print(rep)
    start_time <- Sys.time()
    
    permuted_sample_typisation <- real_sample_typisation
    permuted_sample_typisation$sample <-
      sample(permuted_sample_typisation$sample)
    
    clusterExport(clust, "permuted_sample_typisation")

    simulated_test <-
      parApply(
        cl = clust,
        X = df_CpG,
        MARGIN = 1,
        FUN = function(r) {
          return(meth_vs_type_test2(r, permuted_sample_typisation))
        }
      )
    delta <-
      parApply(
        cl = clust,
        X = df_CpG,
        MARGIN = 1,
        FUN = delta_mean_type
      )
    
    
    random_p_val <-
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
    
    random_statistic <-
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
    
    end_time <- Sys.time()
    st = format(end_time, "%Y-%m-%d_%H.%M.%s")
    print(end_time - start_time)
    
    file_name <-
      paste("Kruskal_Wallis_chi_squared_p_val",
            "type",
            "label_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = random_p_val,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
    file_name <-
      paste(
        "Kruskal_Wallis_chi_squared_statistic",
        "type",
        "label_shuffle",
        st,
        "rds",
        sep = "."
      )
    saveRDS(object = random_statistic,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
    file_name <-
      paste("permuted_sample_typisation",
            "type",
            "label_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = permuted_sample_typisation,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
  }
  stopCluster(clust)
}


if(create_permutations) {
  n_repetitions <- 98
  library("snow")
  #Create cluster
  clust <- makeCluster(parallel::detectCores() - 1)
  # Export it form base to workspace
  clusterExport(clust, "meth_vs_type_test2")
  clusterExport(clust, "sample_types")
  #Apply the declared function
  for (rep in 1:n_repetitions) {
    #
    # rep <- 1
    print(rep)
    start_time <- Sys.time()
    
    permuted_sample_typisation <- real_sample_typisation
    permuted_sample_typisation$sample <-
      sample(permuted_sample_typisation$sample)
    
    clusterExport(clust, "permuted_sample_typisation")
    
    
    
    simulated_test <-
      parApply(
        cl = clust,
        X = df_CpG,
        MARGIN = 1,
        FUN = function(r) {
          return(meth_vs_type_test2(r, permuted_sample_typisation))
        }
      )
    delta <-
      parApply(
        cl = clust,
        X = df_CpG,
        MARGIN = 1,
        FUN = delta_mean_type
      )
    
    
    random_p_val <-
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
    
    random_statistic <-
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
    
    end_time <- Sys.time()
    st = format(end_time, "%Y-%m-%d_%H.%M.%s")
    print(end_time - start_time)
    
    file_name <-
      paste("Kruskal_Wallis_chi_squared_p_val",
            "type",
            "label_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = random_p_val,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
    file_name <-
      paste(
        "Kruskal_Wallis_chi_squared_statistic",
        "type",
        "label_shuffle",
        st,
        "rds",
        sep = "."
      )
    saveRDS(object = random_statistic,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
    file_name <-
      paste("permuted_sample_typisation",
            "type",
            "label_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = permuted_sample_typisation,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
  }
  stopCluster(clust)
}
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
               name = paste("sim", n, sep = "_"))

  sim_results <- rbind(sim_results, temp)
}


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

real_results$name <- "real data"
sim_results$name <- "label permuation"
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

