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
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- "04.methylation_vs_Type"
dir.create(OUTPUT_FOLDER,showWarnings = FALSE)
#INPUT_FOLDER <- 
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


meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
sample_order <- as.character(meta$sample[order(meta$age_mean_BP)])
#age <- meta$age_mean_BP[order(meta$age_mean_BP)]


df_CpG <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/07.methylation_p_val/single_CpG.H3K27ac.hg19.39.samples.p_val.delta_meth.rds")

ggplot(data = df_CpG)

TYPES <- unique(meta$Type)
## only big sample groups 
TYPES <- c("Farmer","Steppe","HG")

sample_types <- list()  
for(t in TYPES){
  print(t)
  sample_types[[t]] <- meta$sample[meta$Type %in% t]
}

CpG_per_type <- data.frame()

  
for(i in 1:length(sample_types) ){
  print(i)
  

  
  for(s in sample_types[[i]]){
    print(s)
    CpG_per_type <- rbind(CpG_per_type , data.frame(meth = df_CpG[,s],type = names(sample_types[i])))
  }
}  
  

## kolmogorov-smirnov-test nomrality test 
df_test <- as.numeric(df_CpG$I2978[complete.cases(df_CpG$I2978)])

ks.test(x = df_test, "pnorm", mean=mean(df_test), sd=sd(df_test))

meth_vs_type_test <- function(df_row){
  df_analyis <- data.frame()
  for(s in names(sample_types)){
    #  s <- names(sample_types)[1]
    meth <- as.numeric(df_row[unlist(sample_types[s],use.names = FALSE)])
    # test prerequisite if at least 5 samples 
    if(sum(!is.na(meth)) < 5){
      return(NA)
    } 
    temp <- data.frame(meth = meth ,type = s)
    df_analyis <- rbind(df_analyis,temp)
  }
  kru_test <- kruskal.test(meth ~ type, data = df_analyis)
  return(kru_test$p.value)
}




test <- apply(X = df_CpG, FUN = meth_vs_type_test,MARGIN = 1)

df_CpG$p_val_kruskal <- test

saveRDS(object = df_CpG,file = file.path(OUTPUT_FOLDER,"single_CpG.H3K27ac.hg19.p_val.kruskal.rds"))

## all cpg 

all_CpG.39.samples.merged.hg19 <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/all_CpG.39.samples.merged.hg19.rds")
all_CpG.39.samples.merged.hg19$p_val_kruskal <- NA

for(chr in CHR_NAMES){
  # chr <- CHR_NAMES[22]
  start.time <- Sys.time()
  print(chr)
  df_all_chrom <- all_CpG.39.samples.merged.hg19[all_CpG.39.samples.merged.hg19$chrom == chr,]
  test_chr <- apply(X = df_all_chrom, FUN = meth_vs_type_test,MARGIN = 1)
  all_CpG.39.samples.merged.hg19$p_val_kruskal[all_CpG.39.samples.merged.hg19$chrom == chr] <- test_chr
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}
saveRDS(object = all_CpG.39.samples.merged.hg19,file = file.path(OUTPUT_FOLDER,"single_CpG.all.hg19.p_val.kruskal.rds"))


library("parallel")
clust <- makeCluster(detectCores() - 1)
clusterExport(cl = clust, list("sample_types"))
start.time <- Sys.time()
a <- parApply(cl = clust, X = df_all_chrom,MARGIN = 1,FUN = meth_vs_type_test)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
stopCluster(clust)
'test_chr <- apply(X = df_all_chrom, FUN = meth_vs_type_test,MARGIN = 1)







###########################################################
df <- all_CpG.39.samples.merged.hg19[!is.na(all_CpG.39.samples.merged.hg19$p_val_kruskal) & all_CpG.39.samples.merged.hg19$p_val_kruskal < 0.01,]
df <- df[!is.na(df$p_val_kruskal) & df$p_val_kruskal < 0.001,]
df <- df[!is.na(df$p_val_kruskal) & df$p_val_kruskal < 0.0001,]

sum(df$name != "NO_CHIP")/nrow(df)

sum(all_CpG.39.samples.merged.hg19$name != "NO_CHIP")/nrow(all_CpG.39.samples.merged.hg19)


i <- 200
df_single <- data.frame()
for(s in names(sample_types)){
  df_temp <- data.frame(meth = as.numeric(x = df[i,sample_types[[s]]]), type = s)
  df_single <- rbind(df_single, df_temp)
}

ggplot(data = df_single,mapping = aes(x = type,y = meth, color = type))+
  geom_boxplot() + 
  geom_point(position = "jitter")+
  theme_minimal()+
  ylim(0,1)

