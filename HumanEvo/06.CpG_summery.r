#################################################################
##
##
##  input:
##           04.methylation_vs_age_CpG
##            
##
##
##  output:    
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
OUTPUT_FOLDER <- "06.CpG_summery"
INPUT_FOLDER <- "04.methylation_vs_age_CpG"
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


meta <- readRDS(file = file.path("03.plots",paste("meta","rds",sep = ".")))
age <- meta$age_mean_BP[order(meta$age_mean_BP)]

df <- readRDS(file = file.path(INPUT_FOLDER,
                                     paste("4b"
                                       ,"H3K27ac",
                                           "with",
                                           "pearson",
                                           N_SAMPLES,
                                           "samples",
                                           ANNOTATION,
                                           "rds",
                                           sep = ".")))

meth <- df[,as.character(levels(meta$sample))]




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
       filename = file.path(OUTPUT_FOLDER,paste("sim vs real pearson_p_val_CpGs","png",sep = ".")),
       width = 6, height = 4)








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


# https://huji.zoom.us/j/89527833780?pwd=anJiMWxYYmRpeDNBRmNDRndENU9NUT09

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