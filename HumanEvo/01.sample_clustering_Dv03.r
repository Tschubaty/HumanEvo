#################################################################
##
##
##  input:
##
##            \HumanEvo\amitai\Matlab_Data\Matlab_data_28.2.19\chr_divided_data
##
##
##  output:   \HumanEvo\HumanEvo\amitai\01.sample_clustering
##
##
##  v_01 14.03.2023
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
library(Rtsne)
library(ggrepel)
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- "clustering"
INPUT_FOLDER <-
  file.path("amitai","Matlab_Data", "Matlab_data_28.2.19", "chr_divided_data")
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
SAMPLE_NAMES = c(
  'I1116',
  'I5319',
  'I4532',
  'I2862',
  'I4529',
  'I3758',
  'I5725',
  'I2861',
  'I3957',
  'I3565',
  'I4792',
  'I4775',
  'I4315',
  'I1053',
  'I7421',
  'I3255',
  'I5835',
  'I2514',
  'I5748',
  'I1633',
  'I5950',
  'I5838',
  'I5742',
  'I5743',
  'I2105',
  'I2520',
  'I2935',
  'I2978',
  'I2980',
  'I3133',
  'I4634',
  'I1965',
  'I1631',
  'I1632',
  'I1961',
  "Stuttgart (LBK)",
  'I1496',
  'I5077',
  'I1962',
  'I4438',
  'I5319',
  'I2134',
  'I1507',
  'I4878',
  'I4432',
  'I4873',
  "Loschbour",
  "La Braña 1",
  'I4596',
  'I5233',
  'I0708',
  'I1960',
  'I4914',
  'I4875',
  'I4877',
  'I2139',
  "SF12",
  'I1734',
  'I5236',
  'I5235',
  "M45 (Houtaomuga)",
  "Ust'-Ishim"
)
SAMPLE_AGE =
  c(
    1,
    1,
    1.2,
    1.6,
    1.8,
    2.2,
    2.5,
    2.8,
    3,
    3.1,
    3.3,
    3.4,
    3.5,
    3.8,
    3.8,
    4,
    4.2,
    4.3,
    4.4,
    4.5,
    4.5,
    4.6,
    4.7,
    4.7,
    4.9,
    5.1,
    5.1,
    5.1,
    5.2,
    5.4,
    5.4,
    5.8,
    6.1,
    6.1,
    6.1,
    7,
    7,
    7,
    7.3,
    7.3,
    7.5,
    7.6,
    7.7,
    7.8,
    7.9,
    7.9,
    8,
    8,
    8,
    8,
    8.1,
    8.1,
    8.1,
    8.4,
    8.5,
    8.9,
    9,
    9.2,
    10,
    10.8,
    12,
    45
  )



################ CODE ###########################


meta_from_XLSX <- as.data.frame(readxl::read_xlsx(path = "Meth samples - human.xlsx"))

meta_from_XLSX <- meta_from_XLSX[,c("Sample","Coverage","tissue","age of sample (kyo)","Culture","Sex","Geographical Location","Deamination %")]

# exclude samples 
include <- meta_from_XLSX$Sample != "Otzi" & meta_from_XLSX$Sample != "Ballito Bay A"
  
meta_from_XLSX <- meta_from_XLSX[include ,]

meta <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/amitai/01.sample_clustering/out/meta.RDS")

meta$Coverage <- NA
meta$tissue <- NA 
meta$"age of sample (kyo)" <- NA
meta$"Culture" <- NA
meta$"Sex" <- NA
meta$"Geographical Location" <- NA
meta$"Deamination %" <- NA


for(r in 1:nrow(meta)){
  meta[r,c("Coverage","tissue","age of sample (kyo)","Culture","Sex","Geographical Location","Deamination %")] <-
    meta_from_XLSX[ SAMPLE_NAMES[r] == meta_from_XLSX$Sample,
                    c("Coverage","tissue","age of sample (kyo)","Culture","Sex","Geographical Location","Deamination %")]
}



# data_type = "sim" #united
# if(data_type == "sim"){
#   filenames <- file.path("Matlab_Data",
#                          "Matlab_data_28.2.19",
#                          list.files(path = file.path("Matlab_Data", "Matlab_data_28.2.19"),
#                                     pattern = "_simu_united.csv"))
# }else{
#   filenames <-
#     file.path(INPUT_FOLDER, list.files(path = INPUT_FOLDER, pattern = ".csv"))
# }
# coordinates<-
#   read.csv(file = file.path("Matlab_Data", "Matlab_data_28.2.19","genomic_locations.csv"),header = F)
# ## initialize meta summery
# meta_colnames <- c("sample", "age" ,paste(CHR_NAMES,"mean_methylation",sep = "_"),paste(CHR_NAMES,"var_methylation",sep = "_"))
# meta <- as.data.frame(matrix(nrow = length(SAMPLE_NAMES),ncol = length(meta_colnames)))
# colnames(meta) <- meta_colnames
# meta$sample <- SAMPLE_NAMES
# meta$age <- SAMPLE_AGE
# # start loop
# for (chr in 1:22) {
#   # debug chr <- 22
#   print(chr)
#   rds_filename <- file.path(OUTPUT_FOLDER,"data",paste(CHR_NAMES[chr],data_type,"rds",sep = "."))
#   if(file.exists(rds_filename)){
#     # load data from rds
#     df_chr <- readRDS(file = rds_filename)
#   }else{ # load data from csv
#     csv_filename <- filenames[grep(pattern = paste("_", chr, "_", sep = ""),
#                                    x = filenames)]
#     df_chr <-
#       read.csv(file = csv_filename,
#                header = F,
#                col.names = SAMPLE_NAMES)
#     rownames(df_chr) <- coordinates[complete.cases(coordinates[,chr]),chr]
#     saveRDS(object = df_chr,file = rds_filename)
#   }
#   meta[,paste(CHR_NAMES[chr],"mean_methylation",sep = "_")] <- colMeans(df_chr,na.rm = T)
#   meta[,paste(CHR_NAMES[chr],"var_methylation",sep = "_")] <- sapply(df_chr,FUN = function(x) var(x,na.rm = T))
# }
# saveRDS(file = file.path(OUTPUT_FOLDER,"out",paste("meta",data_type,"rds",sep = ".")),object = meta)
# ############# plotting ####################
# meta$sample <- factor(x = meta$sample,levels = meta$sample)
# p_age <- ggplot(data = meta,aes(x = sample,y = age)) + geom_point()+
#   ggplot2::theme(
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.background = ggplot2::element_blank(),
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
#   )
# ggsave(filename = file.path(OUTPUT_FOLDER,"plots",paste("age",data_type,"png",sep = ".")),plot = p_age)
# p_meth_chr1 <-ggplot(data = meta,aes(x = sample,y = chr1_mean_methylation,fill = sample)) + geom_bar(stat = "identity")+
#   geom_errorbar(aes(ymin=chr1_mean_methylation - chr1_var_methylation, ymax=chr1_mean_methylation + chr1_var_methylation), width=.2,
#                 position=position_dodge(.9))+
#   ggplot2::theme(
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.background = ggplot2::element_blank(),
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
#   )
# ggsave(filename = file.path(OUTPUT_FOLDER,"plots",paste("mean_methylation_chr1",data_type,"png",sep = ".")),plot = p_meth_chr1)
# ################## correlations ############################
# 
# i <- 1
# 
# cor(x = as.numeric(df_chr[i,]),
#     y = as.vector(SAMPLE_AGE),
#     use = "pairwise.complete.obs",
#     method = c("kendall"))
# 
# 
# cor.test(x = as.numeric(df_chr[i,]),
#          y = SAMPLE_AGE,
#          method = "kendall")
# 
# 
# # p_vals <- apply(
# #   X = df_chr,
# #   MARGIN = 1,
# #   FUN = function(data_row)
# #     if (sum(!is.na(data_row) <= 5)) {
# #       return(NA)
# #     } else{
# #       return(cor.test(
# #         x = as.numeric(data_row),
# #         y = SAMPLE_AGE,
# #         method = "kendall"
# #       ))
# #     }
# # )
# 
# p_vals <- numeric(nrow(df_chr))
# 
# df_p <- data.frame(chr = CHR_NAMES[chr],minus_log_p_vals = -log(p_vals))
# 
# 
# saveRDS(object = df_p,file = file.path(OUTPUT_FOLDER,paste(CHR_NAMES[chr],"p_val","rds",sep = ".")))
# 
# #
# df_p$chr[is.na(df_p$minus_log_p_vals)] <- CHR_NAMES[21]
# df_p$minus_log_p_vals[is.na(df_p$minus_log_p_vals)] <- -log(0.99)
# 
# ggplot(data = df_p[1:10000,],mapping = aes(x = chr,y = minus_log_p_vals))+
#   geom_point(position=position_dodge(width=0.3))+
#   scale_y_continuous(trans='log10')
# 
# position=position_dodge(0.3)
# 
# 
# ggplot(data = df_p[1:578097,],mapping = aes(x =minus_log_p_vals)) +
#   geom_histogram()
# 
# which.max(df_p$minus_log_p_vals)
# 
# 
# 
# df_best <- data.frame(age =SAMPLE_AGE, meth = as.numeric(as.vector(df_chr[which.max(df_p$minus_log_p_vals),])))
# df_best <- data.frame(age =SAMPLE_AGE,
#                       meth = as.numeric(as.vector(
#                         df_chr[
#                           df_p$minus_log_p_vals == head(sort(df_p$minus_log_p_vals,decreasing = T),1000)[400],])))
# 
# 
# #6.164671
# 
# p_best <- ggplot(data = df_best ,aes(x = age,y = meth)) + geom_point()+
#   ggplot2::theme(
#     axis.line = ggplot2::element_line(colour = "black"),
#     panel.background = ggplot2::element_blank(),
#     legend.position = "none",
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
#   )


#   df_chr_unique <- unique(df_chr) # Remove duplicates
#   df_chr_unique_na_omit <-  na.omit(df_chr_unique) # omit NA
#
#   set.seed(42) # Set a seed if you want reproducible results
#
#   perplexity <-  10
#
#
#   #if(file.exists())
#   tsne_out <- Rtsne(X = as.matrix(t(df_chr_unique_na_omit)),
#                     perplexity = perplexity,
#                     num_threads = 7) # Run TSNE
#
#   saveRDS(object = tsne_out, file = file.path(OUTPUT_FOLDER, paste("tsne", CHR_NAMES[1], "RDS", sep = ".")))
#
#   df <-
#     data.frame(
#       x = tsne_out$Y[, 1],
#       y = tsne_out$Y[, 2],
#       age = SAMPLE_AGE,
#       name = SAMPLE_NAMES
#     )
#
#   p <-
#     ggplot(data = df,
#            mapping = aes(
#              x = x,
#              y = y,
#              color = age,
#              label = name
#            )) + geom_point() + geom_label_repel()
#   ggsave(filename = file.path(OUTPUT_FOLDER, paste("tsne", CHR_NAMES[1], "png", sep = ".")),
#          plot = p)
#
# }
