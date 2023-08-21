#################################################################
##
##
##  input:
##           "methylation+chip"
##
##
##
##  output: methylation+chip   
##
##
##  v_01 21.08.2023
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
##################################### CONSTANTS ########################################
INPUT_FOLDER <- "methylation+chip"
OUTPUT_FOLDER <- INPUT_FOLDER 
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

################ CODE ###########################
generate_meta <- FALSE
generate_pearson_cor.methylation.shuffle <- FALSE
generate_pearson_cor.methylation.horizonal_shuffle <- FALSE
generate_pearson_cor.methylation.vertical_column_shuffle <- FALSE

##################################################
if (generate_meta) {
  meta <-
    readRDS(file = file.path("03.plots", paste("meta", "rds", sep = ".")))
  
  df <-
    readRDS(file.path(
      INPUT_FOLDER,
      paste(
        "H3K27ac",
        "only",
        N_SAMPLES,
        "samples",
        "merged",
        ANNOTATION,
        "rds",
        sep = "."
      )
    ))
  # compute in score how many samples
  df$score <-
    apply(
      df[, c(7:ncol(df))],
      1 ,
      FUN = function(r) {
        sum(!is.na(r))
      }
    )
  
  df <- df[df$score == N_SAMPLES, ]
  
  ################################################
  
  age <- meta$age_mean_BP[order(meta$age_mean_BP)]
  df$pearson_p_val <- NA
  df$pearson_cor <- NA
  # for one permutation
  df$radom_p_val <- NA
  #debug
  #df <- head(x = df,n = 1000)
  
  # for(r in 1:nrow(df)){
  #
  #
  #   meth <- df[r,as.character(levels(meta$sample))]
  #
  #   test <- cor.test(x = age,y =  as.numeric(meth),method = "pearson")
  #
  #   df$pearson_p_val[r] <- test$p.value
  #   df$pearson_cor[r] <- test$statistic
  #
  #   permuted.test <- cor.test(x = sample(age),y =  as.numeric(meth),method = "pearson")
  #
  #   df$radom_p_val[r] <- permuted.test$p.value
  #
  # }
  
  meth <- df[, as.character(levels(meta$sample))]
  
  cor_tests <-
    apply(
      meth ,
      MARGIN = 1,
      FUN = function(r) {
        cor.test(x = age,
                 y =  as.numeric(r),
                 method = "pearson")
      }
    )
  
  df$pearson_p_val <- sapply(
    cor_tests ,
    FUN = function(r) {
      r$p.val
    }
  )
  
  df$pearson_cor <- sapply(
    cor_tests ,
    FUN = function(r) {
      r$statistic
    }
  )
  
  
  saveRDS(object = df, file = file.path(
    OUTPUT_FOLDER,
    paste(
      "4b"
      ,
      "H3K27ac",
      "with",
      "pearson",
      N_SAMPLES,
      "samples",
      ANNOTATION,
      "rds",
      sep = "."
    )
  ))
} else{
  meta <-
    readRDS(file = file.path("03.plots", paste("meta", "rds", sep = ".")))
  df <- readRDS(file = file.path(
    OUTPUT_FOLDER,
    paste(
      "4b"
      ,
      "H3K27ac",
      "with",
      "pearson",
      N_SAMPLES,
      "samples",
      ANNOTATION,
      "rds",
      sep = "."
    )
  ))
  meth <- df[, as.character(levels(meta$sample))]
  age <- meta$age_mean_BP[order(meta$age_mean_BP)]
  
}


#############################################################################
# library(parallel)

# mn <- parallel::mclapply(meth[1:10], function(df) {
#   function(r) {
#     cor.test(x = age,
#              y =  as.numeric(r),
#              method = "pearson")
#   }
# }, mc.cores = (parallel::detectCores()-1))

# mn <- parallel::mclapply(meth[1:10], function(df) {
#   function(r) {
#     cor.test(x = age,
#              y =  as.numeric(r),
#              method = "pearson")
#   }
# }, mc.cores = (parallel::detectCores()-1))





# #library("parallel")
# cl <- parallel::makeCluster(parallel::detectCores()-1)
# res <- parallel::parLapply(cl = cl, X = meth[1:10,],
#                            fun = age_corr.test)
# parallel::stopCluster(cl)


n_repetitions <- 100

age_corr.test <- function(r) {
  cor.test(x = age,
           y =  as.numeric(r),
           method = "pearson")
}

if (generate_pearson_cor.methylation.shuffle) {
  for (rep in 1:n_repetitions) {
    # rep <- 1
    print(rep)
    start_time <- Sys.time()
    # vertical permutation
    shuffeld_meth <-
      apply(
        X = meth,
        MARGIN = 2,
        FUN =  function(meth_sample) {
          meth_sample[sample(1:length(meth_sample))]
        }
      )
    
    cor_tests <-
      apply(shuffeld_meth , MARGIN = 1, FUN = age_corr.test)
    
    random_p_val <- sapply(
      cor_tests ,
      FUN = function(r) {
        r$p.val
      }
    )
    
    random_pearson_cor <-
      sapply(
        cor_tests ,
        FUN = function(r) {
          r$statistic
        }
      )
    
    #cor_p.val <- apply(meth , MARGIN = 1, FUN = function(r){cor.test(x = age,y =  as.numeric(r),method = "pearson")$p.val})
    end_time <- Sys.time()
    
    print(end_time - start_time)
    
    st = format(end_time, "%Y-%m-%d_%H.%M")
    file_name <-
      paste("pearson_p_val", "methylation", "shuffle", st, "rds", sep = ".")
    saveRDS(object = random_p_val,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
    file_name <-
      paste("pearson_cor", "methylation", "shuffle", st, "rds", sep = ".")
    saveRDS(object = random_pearson_cor,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
  }
}
#############################################################################
########################### horizonzal permuation ###########################
########################### not parrarelized##############################
#n_repetitions <- 100

# horizontal_random_age_corr.test <- function(r) {
#   cor.test(x = sample(age),
#            y =  as.numeric(r),
#            method = "pearson")}

if (generate_pearson_cor.methylation.horizonal_shuffle) {
  for (rep in 1:n_repetitions) {
    # rep <- 1
    print(rep)
    start_time <- Sys.time()
    # horizontal permutation
    cor_tests <-
      apply(meth , MARGIN = 1, FUN = horizontal_random_age_corr.test)
    
    random_p_val <- sapply(
      cor_tests ,
      FUN = function(r) {
        r$p.val
      }
    )
    
    random_pearson_cor <-
      sapply(
        cor_tests ,
        FUN = function(r) {
          r$statistic
        }
      )
    
    end_time <- Sys.time()
    st = format(end_time, "%Y-%m-%d_%H.%M")
    print(end_time - start_time)
    
    file_name <-
      paste("pearson_p_val",
            "methylation",
            "horizonal_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = random_p_val,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
    file_name <-
      paste("pearson_cor",
            "methylation",
            "horizonal_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = random_pearson_cor,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
  }
  
  
  ############################## parallle R ######################### horizontal ###
  n_repetitions <- 100
  library("snow")
  #Create cluster
  clus <- makeCluster(parallel::detectCores() - 1)
  # Export it form base to workspace
  clusterExport(clus, "age")
  #Apply the declared function
  for (rep in 1:n_repetitions) {
    #
    # rep <- 1
    print(rep)
    start_time <- Sys.time()
    psr_cor_tests <-
      parRapply(
        cl = clus,
        x = meth,
        fun = function(r) {
          cor.test(x = sample(age),
                   y =  as.numeric(r),
                   method = "pearson")
        }
      )
    
    random_p_val <- sapply(
      psr_cor_tests ,
      FUN = function(r) {
        r$p.val
      }
    )
    
    random_pearson_cor <-
      sapply(
        psr_cor_tests ,
        FUN = function(r) {
          r$statistic
        }
      )
    
    end_time <- Sys.time()
    st = format(end_time, "%Y-%m-%d_%H.%M.%s")
    print(end_time - start_time)
    file_name <-
      paste("pearson_p_val",
            "methylation",
            "horizonal_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = random_p_val,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
    file_name <-
      paste("pearson_cor",
            "methylation",
            "horizonal_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = random_pearson_cor,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
  }
  stopCluster(clus)
  #############################################################################
}
############################### Whole COLUMN VEtrical Permuation ############
############################## parallle R #########################
if (generate_pearson_cor.methylation.vertical_column_shuffle) {
  library("snow")
  #Create cluster
  clus <- makeCluster(parallel::detectCores() - 1)
  # Export it form base to workspace
  clusterExport(clus, "age")
  #Apply the declared function
  for (rep in 1:n_repetitions) {
    #
    # rep <- 1
    print(rep)
    start_time <- Sys.time()
    col_permuation <- sample(1:ncol(meth))
    col_permuation_meth <- meth[, col_permuation]
    psr_cor_tests <-
      parRapply(
        cl = clus,
        x = col_permuation_meth,
        fun = function(r) {
          cor.test(x = age,
                   y =  as.numeric(r),
                   method = "pearson")
        }
      )
    
    random_p_val <- sapply(
      psr_cor_tests ,
      FUN = function(r) {
        r$p.val
      }
    )
    
    random_pearson_cor <-
      sapply(
        psr_cor_tests ,
        FUN = function(r) {
          r$statistic
        }
      )
    
    end_time <- Sys.time()
    st =  paste(col_permuation, collapse = "_")
    print(end_time - start_time)
    file_name <-
      paste("pearson_p_val",
            "methylation",
            "vertical_column_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = random_p_val,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
    
    file_name <-
      paste("pearson_cor",
            "methylation",
            "vertical_column_shuffle",
            st,
            "rds",
            sep = ".")
    saveRDS(object = random_pearson_cor,
            file = file.path(OUTPUT_FOLDER, "simulation", file_name))
  }
  stopCluster(clus)
  #############################################################################
}

sim_files <-
  list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
             pattern = "pearson_p_val.methylation.shuffle.*")
horizontal_sim_files <-
  list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
             pattern = "pearson_p_val.methylation.horizonal_shuffle.*")
vertical_column_sim_files <-
  list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
             pattern = "pearson_cor.methylation.vertical_column_shuffle.*")

sim_counts <- list()
horizontal_sim_counts <- list()
vertical_column_sim_counts <- list()

x_val <- list()
horizontal_x_val <- list()
vertical_column_x_val <- list()

#n <- 1

for (n in 1:length(horizontal_sim_files)) {
  print(n)
  sim_pval <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", sim_files[n]))
  horizontal_sim_pval <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", horizontal_sim_files[n]))
  vertical_column_sim_pval <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", vertical_column_sim_files[n]))
  
  
  p <-
    ggplot(data = data.frame(sim_pval), mapping =  aes(x = -log(sim_pval))) +
    geom_histogram(binwidth = 0.25) + # beaks = seq(0,12,0.25)
    theme_minimal() +
    ggtitle("- log hist_sim_pearson_p_val_CpGs for a CpG")
  
  horizontal_p <-
    ggplot(data = data.frame(horizontal_sim_pval),
           mapping =  aes(x = -log(horizontal_sim_pval))) +
    geom_histogram(binwidth = 0.25) + # beaks = seq(0,12,0.25)
    theme_minimal() +
    ggtitle("- log hist_sim_pearson_p_val_CpGs for a CpG")
  
  vertical_column_p <-
    ggplot(data = data.frame(vertical_column_sim_pval),
           mapping =  aes(x = -log(vertical_column_sim_pval))) +
    geom_histogram(binwidth = 0.25) + # beaks = seq(0,12,0.25)
    theme_minimal() +
    ggtitle("- log hist_vertical_column_pearson_p_val_CpGs for a CpG")
  
  ############################
  bp <- ggplot_build(p)
  horizontal_bp <- ggplot_build(horizontal_p)
  vertical_column_bp <- ggplot_build(vertical_column_p)
  
  ############################
  hist_data <- bp$data[[1]]
  horizontal_hist_data <- horizontal_bp$data[[1]]
  vertical_column_hist_data <- vertical_column_bp$data[[1]]
  
  ########################
  sim_counts[[n]] <- hist_data$count
  horizontal_sim_counts[[n]] <- horizontal_hist_data$count
  vertical_column_sim_counts[[n]] <- vertical_column_hist_data$count
  
  ###########################
  x_val[[n]] <- hist_data$x
  horizontal_x_val[[n]] <- horizontal_hist_data$x
  vertical_column_x_val[[n]] <- vertical_column_hist_data$x
}
###########################################################################

sim_results <- data.frame()
horizontal_sim_results <- data.frame()
vertical_column_sim_results <- data.frame()

for (n in 1:length(sim_counts)) {
  temp <-
    data.frame(count = sim_counts[[n]],
               x = x_val[[n]] ,
               name = paste("sim", n, sep = "_"))
  horizontal_temp <-
    data.frame(
      count = horizontal_sim_counts[[n]],
      x = horizontal_x_val[[n]] ,
      name = paste("horizontal_sim", n, sep = "_")
    )
  
  vertical_column_temp <-
    data.frame(
      count = vertical_column_sim_counts[[n]],
      x = vertical_column_x_val[[n]] ,
      name = paste("vertical_column_sim", n, sep = "_")
    )
  
  sim_results <- rbind(sim_results, temp)
  horizontal_sim_results <-
    rbind(horizontal_sim_results, horizontal_temp)
  vertical_column_sim_results <-
    rbind(vertical_column_sim_results, vertical_column_temp)
}
############## sttae of work 

################################################################

#ggplot(data = sim_results,mapping = aes(x = x,y = count,group=name, colour=name))+geom_point(position="dodge")

p <-
  ggplot(data = data.frame(pearson_p_val = df[, c("pearson_p_val")]),
         mapping =  aes(x = -log(pearson_p_val))) +
  geom_histogram(binwidth = 0.25) + # beaks = seq(0,12,0.25)
  theme_minimal() +
  ggtitle("- log hist_sim_pearson_p_val_CpGs for a CpG")

bp <- ggplot_build(p)

hist_data <- bp$data[[1]]

real_results <- hist_data[, c("x", "count")]
#############################################################

real_results$name <- "real data"
horizontal_sim_results$name <- "horizontal (age) permuation"
sim_results$name <- "vertical (methylation) permuation"
vertical_column_sim_results$name <- "whole column (age) permuation"

all_results  <-
  rbind(real_results, horizontal_sim_results, sim_results,vertical_column_sim_results)
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
      "sim horizontal and vertical vs real pearson_p_val_CpGs",
      "png",
      sep = "."
    )
  ),
  width = 16,
  height = 10
)


##################################################################
######################### correlation values #####################
cor_files <-
  list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
             pattern = "pearson_cor.methylation.shuffle.*")
horizontal_cor_files <-
  list.files(path = file.path(OUTPUT_FOLDER, "simulation"),
             pattern = "pearson_cor.methylation.horizonal_shuffle.*")


cor_sim_counts <- list()
cor_x_val <- list()
cor_horizontal_sim_counts <- list()
cor_horizontal_x_val <- list()
cor_sim_results <- data.frame()
cor_horizontal_sim_results <- data.frame()

max_sample <- length(horizontal_sim_files)

for (n in 1:max_sample) {
  # n <- 1
  print(n)
  sim_cor_val <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", cor_files[n]))
  horizontal_sim_cor_val <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", horizontal_cor_files[n]))
  sim_pval <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", sim_files[n]))
  horizontal_sim_pval <-
    readRDS(file.path(OUTPUT_FOLDER, "simulation", horizontal_sim_files[n]))
  
  sim_dataframe <- data.frame(pval = sim_pval, cor = sim_cor_val)
  horizontal_sim_dataframe <-
    data.frame(pval = horizontal_sim_pval, cor = horizontal_sim_cor_val)
  
  
  p <-
    ggplot(data = sim_dataframe, mapping =  aes(x = cor)) + # [sim_dataframe$pval <= uncorrected_alpha,]
    geom_histogram(binwidth = 0.1) + # beaks = seq(0,12,0.25)
    theme_minimal() +
    ggtitle("hist_corr_sim_pearson_CpGs for a CpG")
  
  horizontal_p <-
    ggplot(data = horizontal_sim_dataframe, mapping =  aes(x = cor)) + # [horizontal_sim_dataframe$pval <= uncorrected_alpha,]
    geom_histogram(binwidth = 0.1) + # beaks = seq(0,12,0.25)
    theme_minimal() +
    ggtitle("hist_horizontal_corr_sim_pearson_CpGs for a CpG")
  
  bp <- ggplot_build(p)
  horizontal_bp <- ggplot_build(horizontal_p)
  
  
  cor_hist_data <- bp$data[[1]]
  cor_horizontal_hist_data <- horizontal_bp$data[[1]]
  
  
  cor_sim_counts[[n]] <- cor_hist_data$count
  cor_x_val[[n]] <- cor_hist_data$x
  cor_horizontal_sim_counts[[n]] <- cor_horizontal_hist_data$count
  cor_horizontal_x_val[[n]] <- cor_horizontal_hist_data$x
  
  # put as dataframe
  temp <-
    data.frame(
      count = cor_sim_counts[[n]],
      x = cor_x_val[[n]] ,
      name = paste("cor_sim", n, sep = "_")
    )
  horizontal_temp <-
    data.frame(
      count = cor_horizontal_sim_counts[[n]],
      x = cor_horizontal_x_val[[n]] ,
      name = paste("cor_horizontal_sim", n, sep = "_")
    )
  
  cor_sim_results <- rbind(cor_sim_results, temp)
  cor_horizontal_sim_results <-
    rbind(cor_horizontal_sim_results, horizontal_temp)
  
}
saveRDS(object = cor_sim_results,
        file = file.path(OUTPUT_FOLDER, "cor_sim_results.rds"))
saveRDS(object = cor_horizontal_sim_results,
        file = file.path(OUTPUT_FOLDER, "cor_horizontal_sim_results"))
#cor_sim_results <- readRDS(file.path(OUTPUT_FOLDER,"cor_sim_results.rds"))
#cor_horizontal_sim_results <- readRDS(file.path(OUTPUT_FOLDER,"cor_horizontal_sim_results"))
p <-
  ggplot(data = df, mapping =  aes(x = pearson_cor)) + # [df$pearson_p_val <= uncorrected_alpha,]
  geom_histogram(binwidth = 0.1) + # beaks = seq(0,12,0.25)
  theme_minimal() +
  ggtitle("hist_corr_real_pearson_CpGs for a CpG")


bp <- ggplot_build(p)

hist_data <- bp$data[[1]]

cor_real_results <- hist_data[, c("x", "count")]
#############################################################

cor_real_results$name <- "real data"
cor_horizontal_sim_results$name <- "horizontal (age) permuation"
cor_sim_results$name <- "vertical (methylation) permuation"
cor_all_results  <-
  rbind(cor_real_results,
        cor_horizontal_sim_results,
        cor_sim_results)

file_name <-
  paste("pearson_cor",
        "CpG.methylation",
        "cor_all_results",
        "rds",
        sep = ".")
saveRDS(object = cor_all_results, file = file.path(OUTPUT_FOLDER, file_name))

cor_tail_results <-
  cor_all_results[cor_all_results$x <= -3 | 3 <= cor_all_results$x, ]
#cor_tail_results$x <- as.factor(cor_tail_results$x)

p <-
  ggplot(
    data = cor_tail_results,
    mapping = aes(
      x = x,
      y = count,
      group = name,
      colour = name,
      fill = name
    )
  ) + geom_point(position = "dodge") + theme_minimal()


ggsave(
  filename = file.path(OUTPUT_FOLDER, "cor_tail_results.png"),
  plot = p,
  width = 8,
  height = 5
)



#################################################################
# histogramm # sample coverage
library(tidyr)

df_p_val <- df[, c("pearson_p_val", "radom_p_val")]

df_p_val <- df_p_val %>% gather(key = Type, value = Value)

p <-
  ggplot(data = df_p_val, mapping =  aes(x = -log(Value), fill = Type)) +
  geom_histogram(position = "dodge") +
  theme_minimal() +
  ggtitle("- log hist_all_pearson_p_val_CpGs for a CpG")
ggsave(
  plot = p,
  filename = file.path(
    OUTPUT_FOLDER,
    paste("hist_all_pearson_p_val_CpGs", "png", sep = ".")
  ),
  width = 6,
  height = 4
)


# histogramm # sample coverage
p <-
  ggplot(data = df_p_val, mapping =  aes(x = -log(Value), fill = Type)) +
  geom_histogram(position = "dodge") +
  theme_minimal() +
  ggtitle("- log hist_all_pearson_p_val_CpGs for a CpG")
ggsave(
  plot = p,
  filename = file.path(
    OUTPUT_FOLDER,
    paste("hist_all_pearson_p_val_CpGs", "png", sep = ".")
  ),
  width = 6,
  height = 4
)


# small_pearson <- df_p_val[df_p_val$Type == "pearson_p_val",]
# small_pearson <- small_pearson[order(small_pearson$Value),]
#
# small_random<- df_p_val[df_p_val$Type == "radom_p_val",]
# small_random<- small_random[order(small_random$Value),]

# histogramm # sample coverage
p <-
  ggplot(data = df_p_val[df_p_val$Value <= 0.001, ], mapping =  aes(x = -log(Value), fill =
                                                                      Type)) +
  geom_histogram(position = "dodge") +
  theme_minimal() +
  ggtitle("- log hist_tail_pearson_p_val <= 0.001 for a CpG")
ggsave(
  plot = p,
  filename = file.path(
    OUTPUT_FOLDER,
    paste("hist_tail_pearson_p_val_CpGs", "png", sep = ".")
  ),
  width = 6,
  height = 4
)


r <- sample(which(df$pearson_p_val <= 0.001), size = 1)
example <-
  data.frame(age = age, meth = as.numeric(df[r, as.character(levels(meta$sample))]))

p <- ggplot(data = example, mapping =  aes(x = age, y = meth)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste(df$chrom[r], ".", df$start[r], "-", df$end[r]),
    subtitle = paste("corr", df$pearson_cor[r], "p-val", df$pearson_p_val[r])
  )

p

ggsave(
  plot = p,
  filename = file.path(OUTPUT_FOLDER, paste(
    paste(df$chrom[r], ".", df$start[r], "-", df$end[r]), "png", sep = "."
  )),
  width = 6,
  height = 4
)


# https://huji.zoom.us/j/89527833780?pwd=anJiMWxYYmRpeDNBRmNDRndENU9NUT09

data_meth <- t(df[, as.character(levels(meta$sample))])

pca_prcomp_res <- prcomp(data_meth)

saveRDS(object = pca_prcomp_res,
        file =  file.path(OUTPUT_FOLDER, "pca_prcomp_res.rds"))
#pca_prcomp_res <- readRDS(file = file.path(OUTPUT_FOLDER,"pca_prcomp_res.rds"))

df_pca <- cbind(meta, pca_prcomp_res$x)

eigs <- pca_prcomp_res$sdev ^ 2

pca_explained <- c()

for (i in 1:length(eigs)) {
  pca_explained[i] <- eigs[i] / sum(eigs)
}
# pca1_explained <- eigs[1] / sum(eigs)
#
# pca2_explained <- eigs[2] / sum(eigs)

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
  labs(subtitle = "PCA & sex") + theme_minimal()

ggsave(
  plot = p_pca,
  filename = file.path(OUTPUT_FOLDER, "PCA_vs_sex.png"),
  width = 8,
  height = 6
)

# locality
p_pca <-
  ggplot(data = df_pca , aes(
    x = PC1,
    y = PC2,
    label = sample,
    color = Locality
  )) +
  geom_point() +
  ggrepel::geom_text_repel() +
  labs(subtitle = "PCA & locality")

ggsave(
  plot = p_pca,
  filename = file.path(OUTPUT_FOLDER, "PCA_vs_locality.png"),
  width = 16,
  height = 6
)


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
  labs(subtitle = "PCA & Country")

ggsave(
  plot = p_pca,
  filename = file.path(OUTPUT_FOLDER, "PCA_vs_Country.png"),
  width = 10,
  height = 6
)


##
# "Country"
p_pca <-
  ggplot(data = df_pca , aes(
    x = PC1,
    y = PC2,
    label = sample,
    color = age
  )) +
  geom_point() +
  ggrepel::geom_text_repel() +
  labs(subtitle = "PCA & age")

ggsave(
  plot = p_pca,
  filename = file.path(OUTPUT_FOLDER, "PCA_vs_age.png"),
  width = 10,
  height = 6
)


######################## delta ##############################

df$delta <- as.numeric(apply(X = ,meth,MARGIN = 1,FUN = function(r){max(r)-min(r)}))
nrow(df)
sum(df$pearson_p_val < 0.01,na.rm = TRUE)
sum(df$delta > 0.25,na.rm = TRUE)
sum(df$pearson_p_val < 0.01 & df$delta > 0.25,na.rm = TRUE)

ggplot(data = df,mapping = aes(x = delta))+geom_histogram(binwidth = 0.025)+geom_vline(xintercept = 0.25,color ="red")
hist(df$delta)
hist(-log(df$pearson_p_val))