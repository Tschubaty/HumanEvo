#################################################################
##
##
##  input:
##           "methylation+chip"
##            
##
##
##  output:    04.methylation_vs_age_CpG
##
##
##  v_01 21.04.2023
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
OUTPUT_FOLDER <- "04.methylation_vs_age_CpG"
INPUT_FOLDER <- "methylation+chip"
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
generate_pearson_cor.H3K27ac.peak.shuffle <- FALSE
################ ################################

if(generate_meta){
  
meta <- readRDS(file = file.path("03.plots",paste("meta","rds",sep = ".")))

df <- readRDS(file.path(INPUT_FOLDER,"H3K27ac.peak.mean.methylation.39.samples.merged.hg19.rds"))
# compute in score how many samples 
df$n_samples  <- apply(df[,c(7:ncol(df))], 1 ,FUN = function(r) {sum(!is.na(r))})

df <- df[df$n_samples == N_SAMPLES,]

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

meth <- df[,as.character(levels(meta$sample))]

cor_tests <- apply(meth , MARGIN = 1, FUN = function(r){cor.test(x = age,y =  as.numeric(r),method = "pearson")})

df$pearson_p_val <- sapply(cor_tests , FUN = function(r){r$p.val})

df$pearson_cor <- sapply(cor_tests , FUN = function(r){r$statistic})


saveRDS(object = df,file = file.path(OUTPUT_FOLDER,
                                     paste("4c"
                                       ,"H3K27ac.peak",
                                           "with",
                                           "pearson",
                                           N_SAMPLES,
                                           "samples",
                                           ANNOTATION,
                                           "rds",
                                           sep = ".")))
}else{
  df <- readRDS(file.path(OUTPUT_FOLDER,
                          paste("4c"
                                ,"H3K27ac.peak",
                                "with",
                                "pearson",
                                N_SAMPLES,
                                "samples",
                                ANNOTATION,
                                "rds",
                                sep = ".")))
  meta <- readRDS(file = file.path("03.plots",paste("meta","rds",sep = ".")))
  age <- meta$age_mean_BP[order(meta$age_mean_BP)]
  meth <- df[,as.character(levels(meta$sample))]

}


#############################################################################
n_repetitions <- 100

age_corr.test <- function(r) {
  cor.test(x = age,
           y =  as.numeric(r),
           method = "pearson")}

for(rep in 1:n_repetitions){
  
  print(rep)
  start_time <- Sys.time()
  
  shuffeld_meth <- apply(X = meth, MARGIN = 2,FUN =  function(meth_sample){meth_sample[sample(1:length(meth_sample))]})
  
  cor_tests <- apply(shuffeld_meth , MARGIN = 1, FUN = age_corr.test)
  
  random_p_val <- sapply(cor_tests , FUN = function(r){r$p.val})
  
  random_pearson_cor <- sapply(cor_tests , FUN = function(r){r$statistic})
  
  #cor_p.val <- apply(meth , MARGIN = 1, FUN = function(r){cor.test(x = age,y =  as.numeric(r),method = "pearson")$p.val})
  end_time <- Sys.time()
  
  print(end_time - start_time)
  
  st=format(end_time, "%Y-%m-%d_%H.%M.%s")
  file_name <- paste("pearson_p_val","H3K27ac.peak_methylation", "shuffle",st, "rds", sep = ".")
  saveRDS(object = random_p_val,file = file.path(OUTPUT_FOLDER,"simulation",file_name))

  file_name <- paste("pearson_cor","H3K27ac.peak", "shuffle",st, "rds", sep = ".")
  saveRDS(object = random_pearson_cor,file = file.path(OUTPUT_FOLDER,"simulation",file_name))
}

#############################################################################

sim_files <- list.files(path = file.path(OUTPUT_FOLDER,"simulation"),pattern = "H3K27ac.peak.*")

# df <- readRDS(file = file.path(OUTPUT_FOLDER,
#                                      paste("4c"
#                                            ,"H3K27ac.peak",
#                                            "with",
#                                            "pearson",
#                                            N_SAMPLES,
#                                            "samples",
#                                            ANNOTATION,
#                                            "rds",
#                                            sep = ".")))

sim_counts <- list()
x_val <-list()


#n <- 1 

for(n in 1:length(sim_files)){

print(n)
sim_pval <- readRDS( file.path(OUTPUT_FOLDER,"simulation",sim_files[n]))

p <- ggplot(data = data.frame(sim_pval),mapping =  aes(x = -log(sim_pval)) )+
  geom_histogram(binwidth = 0.25)+ # beaks = seq(0,12,0.25)
  theme_minimal()+
  ggtitle("- log hist_sim_pearson_p_val_CpGs for a CpG")

bp <- ggplot_build(p)

hist_data <- bp$data[[1]]

sim_counts[[n]] <- hist_data$count
x_val[[n]] <- hist_data$x

}

###########################################################################

sim_results <- data.frame()

for(n in 1:length(sim_counts)){
  
  
  temp <- data.frame(count = sim_counts[[n]], x = x_val[[n]] , name = paste("sim",n,sep = "_"))
  
  sim_results <- rbind(sim_results,temp)
}




################################################################

#ggplot(data = sim_results,mapping = aes(x = x,y = count,group=name, colour=name))+geom_point(position="dodge")

p <- ggplot(data = data.frame(pearson_p_val = df[,c("pearson_p_val")]),mapping =  aes(x = -log(pearson_p_val)) )+
  geom_histogram(binwidth = 0.25)+ # beaks = seq(0,12,0.25)
  theme_minimal()+
  ggtitle("- log hist_sim_pearson_p_val_CpGs for a CpG")

bp <- ggplot_build(p)

hist_data <- bp$data[[1]]

real_results <- hist_data[,c("x","count")]
#############################################################

min_log_p <- 5

p <- ggplot(data = sim_results[sim_results$x > min_log_p,])+geom_point(position="dodge",mapping = aes(x = x,y = count,group=name, colour=name))+
  geom_point(data=real_results[real_results$x > min_log_p,],mapping = aes(x=x,y=count),size = 4,shape = 1)


ggsave(plot = p,
       filename = file.path(OUTPUT_FOLDER,paste("sim vs real ChiP peak values  pearson_p_val_CpGs","png",sep = ".")),
       width = 20, height = 16)

##################################################################

# histogramm # sample coverage 
library(tidyr)

df_p_val<- df[,c("pearson_p_val", "radom_p_val")]

df_p_val<- df_p_val %>% gather(key=Type, value=Value) 

p <- ggplot(data = df_p_val,mapping =  aes(x = -log(Value),fill=Type) )+
  geom_histogram(position="dodge")+
  theme_minimal()+
  ggtitle("- log hist_all_pearson_p_val_CpGs for a CpG")
ggsave(plot = p,
       filename = file.path(OUTPUT_FOLDER,paste("hist_all_pearson_p_val_CpGs","png",sep = ".")),
       width = 6, height = 4)


# histogramm # sample coverage 
p <- ggplot(data = df_p_val,mapping =  aes(x = -log(Value),fill=Type) )+
  geom_histogram(position="dodge")+
  theme_minimal()+
  ggtitle("- log hist_all_pearson_p_val_CpGs for a CpG")
ggsave(plot = p,
       filename = file.path(OUTPUT_FOLDER,paste("hist_all_pearson_p_val_CpGs","png",sep = ".")),
       width = 6, height = 4)


# small_pearson <- df_p_val[df_p_val$Type == "pearson_p_val",]
# small_pearson <- small_pearson[order(small_pearson$Value),]
# 
# small_random<- df_p_val[df_p_val$Type == "radom_p_val",]
# small_random<- small_random[order(small_random$Value),]

# histogramm # sample coverage 
p <- ggplot(data = df_p_val[df_p_val$Value <= 0.001,],mapping =  aes(x = -log(Value),fill=Type) )+
  geom_histogram(position="dodge")+
  theme_minimal()+
  ggtitle("- log hist_tail_pearson_p_val <= 0.001 for a CpG")
ggsave(plot = p,
       filename = file.path(OUTPUT_FOLDER,paste("hist_tail_pearson_p_val_CpGs","png",sep = ".")),
       width = 6, height = 4)


r <- sample(which(df$pearson_p_val <= 0.001),size = 1)
example <- data.frame(age = age,meth =as.numeric(df[r,as.character(levels(meta$sample))]))

p <- ggplot(data = example,mapping =  aes(x = age,y = meth ))+
  geom_point()+
  theme_minimal()+
  labs(title= paste(df$chrom[r],".",df$start[r],"-",df$end[r]),
       subtitle=paste("corr",df$pearson_cor[r],"p-val",df$pearson_p_val[r]))

p

ggsave(plot = p,
       filename = file.path(OUTPUT_FOLDER,paste(paste(df$chrom[r],".",df$start[r],"-",df$end[r]),"png",sep = ".")),
       width = 6, height = 4)


data_meth <- t(df[,as.character(levels(meta$sample))])

pca_prcomp_res <- prcomp(data_meth)

saveRDS(object = pca_prcomp_res,file =  file.path(OUTPUT_FOLDER,"pca_prcomp_res.rds"))


df_pca <- cbind(meta,pca_prcomp_res$x)

# sex
p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = sex))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA & sex")

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_vs_sex.png"),width=8, height=6)

# locality
p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = Locality))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA & locality")

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_vs_locality.png"),width=16, height=6)


##
# "Country" 
p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = Country))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA & Country")

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_vs_Country.png"),width=10, height=6)


##
# "Country" 
p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = age))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA & age")

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_vs_age.png"),width=10, height=6)