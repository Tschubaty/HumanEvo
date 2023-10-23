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
#library("parallel")
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- "11.PCA"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
INPUT_FOLDER <-"04.methylation_vs_Type"
N_SAMPLES <- 37
ANNOTATION <- "hg19"

df_CpG <- readRDS(file =  file.path(INPUT_FOLDER,"peak_CpG_with_kruskal_valis_p_val.rds"))
meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
real_sample_typisation <- readRDS(file.path(INPUT_FOLDER,"real_sample_typisation.rds"))
meta37 <- meta[meta$sample %in% real_sample_typisation$sample,]

data_meth <- t(df_CpG[complete.cases(df_CpG),real_sample_typisation$sample])
print("erase non complete rows")
print(sum(complete.cases(df_CpG)/nrow(df_CpG)))

pca_prcomp_res <- prcomp(data_meth)

saveRDS(object = pca_prcomp_res,file =  file.path(OUTPUT_FOLDER,paste("pca_prcomp_res","N_SAMPLES","37","rds",sep = ".")))


df_pca <- cbind(meta37,pca_prcomp_res$x)

# sex
p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = sex))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA_1v2 & sex")+
  theme_minimal()

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_1v2_vs_sex_N_SAMPLES_37.png"),width=8, height=6)

# sex
p_pca <- ggplot(data = df_pca ,aes(x=PC3, y=PC4,label = sample,color = sex))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA_3v4 & sex")+
  theme_minimal()

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_3v4_vs_sex_N_SAMPLES_37.png"),width=8, height=6)


# locality
# p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = Locality))+
#   geom_point()+
#   ggrepel::geom_text_repel()+
#   labs(subtitle="PCA_1v2 & locality")
# 
# ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_1v2_vs_locality_N_SAMPLES_37.png"),width=16, height=6)


##
# "Country" 
p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = Country))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA_1v2 & Country")+
  theme_minimal()

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_1v2_vs_Country.png"),width=10, height=6)

p_pca <- ggplot(data = df_pca ,aes(x=PC3, y=PC4,label = sample,color = Country))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA_3v4 & Country")+
  theme_minimal()

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_3v4_vs_Country.png"),width=10, height=6)


##
# "age" 
p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = age_mean_BP))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA_1v2 & age")+
  theme_minimal()

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_1v2_vs_age.png"),width=10, height=6)

p_pca <- ggplot(data = df_pca ,aes(x=PC3, y=PC4,label = sample,color = age_mean_BP))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA_3v4 & age")+
  theme_minimal()

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_3v4_vs_age.png"),width=10, height=6)

## Type

p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = Type))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA_1v2 & Type")+
  theme_minimal()

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_1v2_vs_Type.png"),width=10, height=6)

p_pca <- ggplot(data = df_pca ,aes(x=PC3, y=PC4,label = sample,color = Type))+
  geom_point()+
  ggrepel::geom_text_repel()+
  labs(subtitle="PCA_3v4 & Type")+
  theme_minimal()

ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_3v4_vs_Type.png"),width=10, height=6)


#########################################

#t-sne 
# Load the required packages
library(Rtsne)
library(ggplot2)
# Load the default dataset
#data(iris)
dim(data_meth)


# Remove Duplicate data present in iris
# data set(Otherwise Error will be generated)
#remove_iris_dup <- unique(iris)
remove_CpG_dup <- unique(data_meth)
dim(remove_CpG_dup)

# Forming the matrix for the first four columns 
meth_matrix <- as.matrix(remove_CpG_dup)


# Calculate tSNE using Rtsne(0 function) 
tsne_out <- Rtsne( X = meth_matrix,
                   perplexity = floor((nrow(meth_matrix) - 1) / 3),
                   dims = 2 )

# Conversion of matrix to dataframe
meta37$tsne_1 <- tsne_out$Y[,1]
meta37$tsne_2 <- tsne_out$Y[,2]

# Plotting the plot using ggplot() function
ggplot2::ggplot(data = meta37,mapping = aes(x = tsne_1, y = tsne_2,label=sample,color = sex)) + # 
  geom_point()+  
  ggrepel::geom_text_repel()+
  labs(subtitle="tsne")+
  theme_minimal()

# Plotting the plot using ggplot() function
ggplot2::ggplot(data = meta37,mapping = aes(x = tsne_1, y = tsne_2,label=sample,color = Country)) + # 
  geom_point()+  
  ggrepel::geom_text_repel()+
  labs(subtitle="tsne")+
  theme_minimal()

# Plotting the plot using ggplot() function
ggplot2::ggplot(data = meta37,mapping = aes(x = tsne_1, y = tsne_2,label=sample,color = Type)) + # 
  geom_point()+  
  ggrepel::geom_text_repel()+
  labs(subtitle="tsne")+
  theme_minimal()
