#################################################################
##
##
##  input:
##          bestCpG <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/10.optimal_p_val_and_delta_AGE/min_delta_0.39_max_minus_log_p_9.3/bestCpG.rds")
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
OUTPUT_FOLDER <- "05.annotation_CpG"
dir.create(OUTPUT_FOLDER, showWarnings = FALSE)
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

min_delta <- 0.39
minus_log_p <- 9.3

meta <- as.data.frame(readxl::read_xlsx("AGDP.39.metadata2.xlsx"))
sample_order <- as.character(meta$sample[order(meta$age_mean_BP)])
age <- meta$age_mean_BP[order(meta$age_mean_BP)]

file_name <- paste("min_delta",min_delta,"max_minus_log_p",minus_log_p,sep = "_")
df_CpG <- readRDS(file.path("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/10.optimal_p_val_and_delta_AGE",file_name,"bestCpG.rds"))
nrow(df_CpG)


bed4 <- df_CpG[,c("chrom","start","end","name")]

write.table(x = bed4,
            file = file.path(OUTPUT_FOLDER,paste("delta",min_delta,"results","p",minus_log_p,"bed",sep = ".")),
            quote = FALSE,row.names = FALSE,col.names = FALSE)


### make Genhance query with files https://genome.ucsc.edu/cgi-bin/hgTables

# Database: hg19    Primary Table: geneHancerInteractionsDoubleElite Data last updated: 2019-01-15
# Big Bed File Download: /gbdb/hg19/geneHancer/geneHancerInteractionsDoubleElite.v2.hg19.bb


df_GH <- read.delim(file =  file.path(OUTPUT_FOLDER,paste(min_delta,"delta.GH-Interactions","p",minus_log_p,"bed",sep = ".")),
                  header = FALSE,
                  skip = 1,
                  col.names = colnames(df_CpG)[1:5])

df_CpG$GeneAnnotation <-NA
df_CpG$EnhancerAnnotation <- NA

for(r in 1:nrow(df_CpG)){
  # r <- 1
  annotated_GH <- df_GH[df_CpG$chrom[r] == df_GH$chrom &  df_GH$start <= df_CpG$start[r] & df_CpG$end[r] <= df_GH$end,]
  if(nrow(annotated_GH) > 0){
    

  annotation <- unlist(strsplit(annotated_GH$name[which.max(annotated_GH$score)], split = "/")) #annotated_GH$name[which.max(annotated_GH$score)]
  df_CpG$GeneAnnotation[r] <- annotation[1]
  df_CpG$EnhancerAnnotation[r] <- annotation[2]
  }
}



write.csv2(
  x = unique(df_CpG$GeneAnnotation),
  file = file.path(
    OUTPUT_FOLDER,
    paste(
      "unique.genes.age_vs_meth",
      min_delta,
      "delta.GH_Interactions",
      "p",
      minus_log_p,
      "csv",
      sep = "."
    )
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)




################################### end review 

# plot all CpG
if(FALSE){
for(b in 1:nrow(df_CpG)){
  # b <- 1
  
  dir_name <- file.path(OUTPUT_FOLDER,"correlation_plots", paste("delta",min_delta,"logP",minus_log_p,sep = "_"))
   dir.create(path = dir_name,showWarnings = FALSE)

  meth <- as.numeric(df_CpG[b,sample_order])
  df_plot <- data.frame(methylation = as.numeric(meth), age= as.numeric(age))
  
  plotname <- paste(df_CpG$GeneAnnotation[b], df_CpG$chrom[b],df_CpG$start[b],sep = " ")
  
  p <- ggplot(data = df_plot ,mapping =  aes(x = age, y = methylation)) + 
    geom_point() +
    ggtitle(label = plotname,subtitle = paste("p_val",format.pval(df_CpG$pearson_p_val[b]),sep = " "))+
    geom_smooth(method="lm")+
    theme_minimal()+
    ylim(0,1)+
    scale_x_continuous(name = "years before present (YBP)", breaks = seq(from=1000,by=1000,to=11000), limits = c(500, 11000),guide = guide_axis(n.dodge=2))
  
  
  ggsave(filename =file.path(dir_name,paste(plotname ,"png",sep = ".")),plot = p)

}
}
  
#   p <- ggplot(data = df_CpG,mapping = aes(x = start,y = pearson_cor))+
#     geom_point()+
#     geom_point(data = df_peak[df_peak$p_val == df_sig$p_val[b],],mapping = aes(color = "red"))+
#     xlab(df_sig$name[b])+
#     theme_minimal()+
#     ggtitle(df_sig$annotation_narrow[b])+ theme(legend.position = "none")
#   ggsave(filename = file.path(dir_name,paste(df_CpG$GeneAnnotation[b],"png",sep = ".")),plot = p,width = 4,height = 4)
#   
#   
#   for(r in 1:nrow(df_peak)){
#     meth <- df_peak[r,as.character(levels(meta$sample))]
#     df_plot <- data.frame(methylation = as.numeric(meth), age= age)
#     
#     scatter_plot <- ggplot(data = df_plot ,mapping =  aes(x = age, y = methylation)) +
#       geom_point() +
#       ggtitle(paste(df_peak$chrom[r],df_peak$start[r],"p_val:",format.pval(df_peak$p_val[r])))+
#       geom_smooth(method="lm") #labs(x = "foot length (cm)", y = "height (cm)")
#     ggsave(filename = file.path(OUTPUT_FOLDER,"correlation plots", df_sig$name[b],paste(df_peak$start[r],"cor",round(df_peak$pearson_cor[r],2),"png",sep = ".")),plot = scatter_plot)
#   }
# }
#   
# # best CpG
# dir.create(file.path(OUTPUT_FOLDER,"correlation plots", "best_CpG"),showWarnings = FALSE)
# 
# for(b in 1:nrow(df_sig)){
#   # b <- 1
#   df_peak <- df_CpG[df_CpG$name == df_sig$name[b],]
# 
#   # p <- ggplot(data = df_peak,mapping = aes(x = start,y = pearson_cor))+
#   #   geom_point()+
#   #   geom_point(data = df_peak[df_peak$p_val == df_sig$p_val[b],],mapping = aes(color = "red"))+
#   #   xlab(df_sig$name[b])+
#   #   theme_minimal()+
#   #   ggtitle(df_sig$annotation_narrow[b])+ theme(legend.position = "none")
#   # ggsave(filename = file.path(OUTPUT_FOLDER,"correlation plots", df_sig$name[b],paste(df_sig$annotation_narrow[b],"png",sep = ".")),plot = p,width = 4,height = 4)
# 
#     index <- df_peak$p_val == df_sig$p_val[b]
#     
#     meth <- df_peak[index,as.character(levels(meta$sample))]
#     df_plot <- data.frame(methylation = as.numeric(meth), age= age)
#     
#     scatter_plot <- ggplot(data = df_plot ,mapping =  aes(x = age, y = methylation)) + 
#       geom_point() +
#       ggtitle(paste(df_peak$chrom[index],df_peak$start[index],"p_val:",format.pval(df_peak$p_val[index])))+
#       geom_smooth(method="lm")+
#       theme_minimal()+
#       ylim(0,1)+
#       scale_x_continuous(name = "years before present (YBP)", breaks = seq(from=1000,by=1000,to=11000), limits = c(500, 11000),guide = guide_axis(n.dodge=2))
#     ggsave(filename =file.path(OUTPUT_FOLDER,"correlation plots", "best_CpG",paste(df_sig$annotation_narrow[b] ,"png",sep = ".")),plot = scatter_plot)
# }
# 
# write.table(x = df_sig,
#             file = file.path(OUTPUT_FOLDER,paste("delta_0.25","annotaed_peak_with_best_CpG_results","p",alpha,"bed",sep = ".")),
#             quote = FALSE,row.names = FALSE,col.names = FALSE)
# 
# 
# write.csv2(
#   x = unique(df_sig$annotation_narrow),
#   file = file.path(
#     OUTPUT_FOLDER,
#     paste(
#       "delta_0.25",
#       "annotation_narrow_gene_names",
#       "p",
#       alpha,
#       "txt",
#       sep = "."
#     )
#   ),
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = FALSE
# )
# 
# openxlsx::write.xlsx(x = meta[,c("sample","Coverage","age_mean_BP", "age_std_BP","Group ID", "Locality","Country", "sex")] , file = file.path("AGDP.39.metadata.xlsx"))
# 
# # anaylise trend
# 
# df_all <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/07.methylation_p_val/all.hg19.39.samples.p_val.delta_meth.rds")
# 
# df_test <- df_all[!is.na(df_all$p_val) & df_all$delta_meth >= min_delta & df_all$p_val <= alpha,]
# 
# ggplot(data = df_test,mapping = aes(x = pearson_cor))+geom_histogram(breaks = seq(from=-8,by=0.1,to=8))+xlim(-8,8)+ggtitle("general correlation bias towards negative correlation")
# 
# # openxlsx::write.xlsx(x = df_sig , file = file.path(OUTPUT_FOLDER,paste("delta_0.25","annotaed_peak_with_best_CpG_results","p",alpha,"xlsx",sep = ".")))
# # 
# # 
# # 
# # nrow(df_sig)
# # 
# # length(unique(df_sig$annotation_narrow_GH))
# # 
# # length(unique(df_sig$annotation_narrow))
# # 
# # length(unique(df_sig$annotation_common))
# # 
# # sort(table(df_sig$annotation_narrow_GH),decreasing = T)
# # 
# # sort(table(df_sig$annotation_narrow),decreasing = T)
# # 
# # genes_for_gene_organizer <- sort(table(df_sig$annotation_narrow),decreasing = T)
# # 
# # text_out <- names(genes_for_gene_organizer[genes_for_gene_organizer > 1])
# # 
# # write.csv(x = text_out,
# #           file = file.path(OUTPUT_FOLDER,
# #                            paste("Genes","narrow","GH-Interactions","p",alpha,"txt",sep = ".")),
# #           quote = FALSE,
# #           row.names = FALSE
# #           )
# # 
# # 
# # # histogramm # sample coverage 
# # library(tidyr)
# # 
# # df_p_val<- df[,c("pearson_p_val", "radom_p_val")]
# # 
# # df_p_val<- df_p_val %>% gather(key=Type, value=Value) 
# # 
# # p <- ggplot(data = df_p_val,mapping =  aes(x = -log(Value),fill=Type) )+
# #   geom_histogram(position="dodge")+
# #   theme_minimal()+
# #   ggtitle("- log hist_all_pearson_p_val_CpGs for a CpG")
# # ggsave(plot = p,
# #        filename = file.path(OUTPUT_FOLDER,paste("hist_all_pearson_p_val_CpGs","png",sep = ".")),
# #        width = 6, height = 4)
# # 
# # 
# # # histogramm # sample coverage 
# # p <- ggplot(data = df_p_val,mapping =  aes(x = -log(Value),fill=Type) )+
# #   geom_histogram(position="dodge")+
# #   theme_minimal()+
# #   ggtitle("- log hist_all_pearson_p_val_CpGs for a CpG")
# # ggsave(plot = p,
# #        filename = file.path(OUTPUT_FOLDER,paste("hist_all_pearson_p_val_CpGs","png",sep = ".")),
# #        width = 6, height = 4)
# # 
# # 
# # # small_pearson <- df_p_val[df_p_val$Type == "pearson_p_val",]
# # # small_pearson <- small_pearson[order(small_pearson$Value),]
# # # 
# # # small_random<- df_p_val[df_p_val$Type == "radom_p_val",]
# # # small_random<- small_random[order(small_random$Value),]
# # 
# # # histogramm # sample coverage 
# # p <- ggplot(data = df_p_val[df_p_val$Value <= alpha,],mapping =  aes(x = -log(Value),fill=Type) )+
# #   geom_histogram(position="dodge")+
# #   theme_minimal()+
# #   ggtitle("- log hist_tail_pearson_p_val <= alpha for a CpG")
# # 
# # ggsave(plot = p,
# #        filename = file.path(OUTPUT_FOLDER,paste("hist_tail_pearson_p_val_CpGs","png",sep = ".")),
# #        width = 6, height = 4)
# # 
# # 
# # r <- sample(which(df$pearson_p_val <= alpha),size = 1)
# # example <- data.frame(age = age,meth =as.numeric(df[r,as.character(levels(meta$sample))]))
# # 
# # p <- ggplot(data = example,mapping =  aes(x = age,y = meth ))+
# #   geom_point()+
# #   theme_minimal()+
# #   labs(title= paste(df$chrom[r],".",df$start[r],"-",df$end[r]),
# #        subtitle=paste("corr",df$pearson_cor[r],"p-val",df$pearson_p_val[r]))
# # 
# # p
# # 
# # ggsave(plot = p,
# #        filename = file.path(OUTPUT_FOLDER,paste(paste(df$chrom[r],".",df$start[r],"-",df$end[r]),"png",sep = ".")),
# #        width = 6, height = 4)
# # 
# # 
# # # https://huji.zoom.us/j/89527833780?pwd=anJiMWxYYmRpeDNBRmNDRndENU9NUT09
# # 
# # data_meth <- t(df[,as.character(levels(meta$sample))])
# # 
# # pca_prcomp_res <- prcomp(data_meth)
# # 
# # saveRDS(object = pca_prcomp_res,file =  file.path(OUTPUT_FOLDER,"pca_prcomp_res.rds"))
# # 
# # 
# # df_pca <- cbind(meta,pca_prcomp_res$x)
# # 
# # # sex
# # p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = sex))+
# #   geom_point()+
# #   ggrepel::geom_text_repel()+
# #   labs(subtitle="PCA & sex")
# # 
# # ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_vs_sex.png"),width=8, height=6)
# # 
# # # locality
# # p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = Locality))+
# #   geom_point()+
# #   ggrepel::geom_text_repel()+
# #   labs(subtitle="PCA & locality")
# # 
# # ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_vs_locality.png"),width=16, height=6)
# # 
# # 
# # ##
# # # "Country" 
# # p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = Country))+
# #   geom_point()+
# #   ggrepel::geom_text_repel()+
# #   labs(subtitle="PCA & Country")
# # 
# # ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_vs_Country.png"),width=10, height=6)
# # 
# # 
# # ##
# # # "Country" 
# # p_pca <- ggplot(data = df_pca ,aes(x=PC1, y=PC2,label = sample,color = age))+
# #   geom_point()+
# #   ggrepel::geom_text_repel()+
# #   labs(subtitle="PCA & age")
# # 
# # ggsave(plot = p_pca,filename = file.path(OUTPUT_FOLDER,"PCA_vs_age.png"),width=10, height=6)
# # 
# # #########################################
# # 
# # df_peak <- readRDS("C:/Users/Daniel Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/H3K27ac.peak.mean.methylation.39.samples.merged.hg19.rds")
# # 
# # df_peak <- df_peak[complete.cases(df_peak),]
# # 
# # alpah2 <- 0.05
# # 
# # df_peak$percent_sig
# # 
# # for(r in 1:nrow(df_peak)){
# #   
# #   inside <- df[df_peak$chr[r] == df$chrom & df_peak$start[r] <= df$start & df$end <= df_peak$end[r],]
# #   if(nrow(inside) < 1){
# #     df_peak$n_sig[r] <- 0
# #   }else{
# #     df_peak$n_sig[r] <- sum(inside$pearson_p_val <= alpah2)
# #   }
# #   
# # }
# # 
# # df_peak$percent_sig  <- df_peak$n_sig/df_peak$score 
# # 
# # hist(df_peak$score)
# # 
# # hist(df_peak$percent_sig[df_peak$score > 5] )
# # 
# # saveRDS(object = df_peak,file =  file.path(OUTPUT_FOLDER,"H3K27ac.complete.data.peak.39.samples.hg19.rds"))
# # 
# # ################################################
# # 
# # age <- meta$age_mean_BP[order(meta$age_mean_BP)]
# # 
# # my_peaks <- df_peak[df_peak$percent_sig > 0.90,]
# # 
# # my_peaks <- my_peaks[complete.cases(my_peaks$percent_sig),]
# # 
# # my_peaks <- my_peaks[order(my_peaks$score,decreasing = TRUE),] 
# # 
# # r <- 3
# # 
# # inside <- df[my_peaks$chr[r] == df$chrom & my_peaks$start[r] <= df$start & df$end <= my_peaks$end[r],]
# # 
# # 
# # #chr1  84326322  84326753 H3K27ac    36      . 0.25874001 0.1916754 0.22181286 0.16322068 0.1733193 0.14395641 0.23695774 0.3790677 0.09258414 0.086213526 0.165020625
# # 
# # cpg <- 1
# # 
# # dir.create(file.path(OUTPUT_FOLDER,inside$name[cpg]))
# # 
# # p <- ggplot(data = inside ,mapping =  aes(x = start,y = -log(pearson_p_val) ))+
# #   geom_bar(stat = "identity")+
# #   theme_minimal()+
# #   labs(title= inside$name[cpg])
# # 
# # ggsave(plot = p,
# #        filename = file.path(OUTPUT_FOLDER,inside$name[cpg],paste("pearson_p_val","png",sep = ".")),
# #        width=10, height=6)
# # 
# # p <- ggplot(data = inside ,mapping =  aes(x = start,y = pearson_cor ))+
# #   geom_bar(stat = "identity")+
# #   theme_minimal()+
# #   labs(title= inside$name[cpg])
# # 
# # ggsave(plot = p,
# #        filename = file.path(OUTPUT_FOLDER,inside$name[cpg],paste("pearson_cor","png",sep = ".")),
# #        width=10, height=6)
# # 
# # 
# # for(cpg in 1:nrow(inside)){
# #   # 
# #   example <- data.frame(age = age,meth =as.numeric(inside[cpg,as.character(levels(meta$sample))]))
# #   
# #   
# #   p <- ggplot(data = example,mapping =  aes(x = age,y = meth ))+
# #     geom_point()+
# #     theme_minimal()+
# #     labs(title= paste(inside$chrom[cpg],".",inside$start[cpg],"-",inside$end[cpg]),
# #          subtitle=paste("corr",inside$pearson_cor[cpg],"p-val",inside$pearson_p_val[cpg]))
# #   
# #   ggsave(plot = p,
# #          filename = file.path(OUTPUT_FOLDER,inside$name[cpg],paste(inside$chrom[cpg],inside$start[cpg],"png",sep = ".")),
# #          width=10, height=6)
# #   
# # }
# # 
# # 
# # 
