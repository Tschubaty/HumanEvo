#################################################################
##
##
##  input:
##           "methylation+chip"
##            
##
##
##  output:    03.plots
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
OUTPUT_FOLDER <- "03.plots"
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

meta  <- as.data.frame(readxl::read_xlsx(path = "AGDP.metadata.xlsx"))

colnames(meta)[1] <- "sample"
colnames(meta)[3] <- "age_mean_BP"
colnames(meta)[4] <- "age_std_BP"
colnames(meta)[10] <- "sex"


meta <- meta[,c("sample", "Coverage" ,"age_mean_BP", "age_std_BP","Group ID",                                                                                                                                                                                                                             
                          "Locality" ,                                                                                                                                                                                                                            
                          "Country" ,                                                                                                                                                                                                                             
                          "sex")]


df <- readRDS(file.path(INPUT_FOLDER,
                            paste("all_CpG",N_SAMPLES,"samples","merged",ANNOTATION,"rds",sep = ".")))

# compute in score how many samples 
df$score <- apply(df[,c(7:ncol(df))], 1 ,FUN = function(r) {sum(!is.na(r))})

# compute in score how many samples 
df$name[df$name != "NO_CHIP"] <- "H3K27ac"

df_pie <- table(df$name)
df_pie <- data.frame(df_pie)
df_pie$pr

colnames(df_pie)[1] <- "group"

# Basic piechart
p <- ggplot(df_pie , aes(x="", y=Freq, fill=group)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_void()+
  geom_text(aes(label = Freq),
            position = position_stack(vjust = 0.5))+
   ggtitle("all CpGs in genome")

print(df_pie$Freq[1]/df_pie$Freq[2]*100)

ggsave(plot = p,filename = file.path(OUTPUT_FOLDER,paste("pie_plot_all_CpGs","png",sep = ".")))

# H3K27ac
df_H3K27ac <- df[df$name != "NO_CHIP",]

# histogramm # sample coverage 
p <- ggplot(data = df_H3K27ac,mapping =  aes(score) )+geom_histogram(bins = N_SAMPLES+1)+theme_minimal()+ggtitle("# samples with values for a CpG")
ggsave(plot = p,filename = file.path(OUTPUT_FOLDER,paste("hist_plot_H3K27ac_CpGs","png",sep = ".")))
# NO DATA for I2978
#[1] TRUE TRUE

H3K27ac_means <- apply(df[df$name == "H3K27ac",c(7:ncol(df))], 2 ,FUN = function(s) {mean(s,na.rm = TRUE)})
methylation_means <- as.data.frame(H3K27ac_means)
methylation_means$sample <- rownames(methylation_means)
colnames(methylation_means)[colnames(methylation_means) == "H3K27ac_means"] <- "mean_methylation"
methylation_means$ChIP <- "H3K27ac"


temp <- data.frame(matrix(nrow = N_SAMPLES,ncol = 3))

for(s in 1:N_SAMPLES){
  temp[s,] <- c(mean(x =df[,s+6] ,na.rm = TRUE),colnames(df)[s+6],"All")
}
#NO_CHIP_means <- apply(df[df$name == "NO_CHIP",c(7:ncol(df))], 2 ,FUN = function(s) {mean(s,na.rm = TRUE)})

colnames(temp) <- colnames(methylation_means)

methylation_means <- rbind(methylation_means,temp )

row.names(methylation_means) <- NULL
methylation_means$mean_methylation <- as.numeric(methylation_means$mean_methylation )

p <- ggplot(data = methylation_means ,mapping = aes(x = ChIP,y = mean_methylation,fill = ChIP))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = "jitter")+
  theme_minimal()+
  ggtitle("")


ggsave(plot = p,filename = file.path(OUTPUT_FOLDER,paste("mean meth_all_CpGs","png",sep = ".")))

############# meta plots ###################################################

meta <- meta[meta$sample %in% unique(methylation_means$sample),]
meta$mean_mean_methylation_H3K27ac <- NA
meta$mean_mean_methylation_all <- NA

for(r in 1:nrow(meta)){
  
  meta$mean_mean_methylation_H3K27ac[r] <- methylation_means$mean_methylation[meta$sample[r] == methylation_means$sample & methylation_means$ChIP == "H3K27ac"]
  meta$mean_mean_methylation_all[r] <- methylation_means$mean_methylation[meta$sample[r] == methylation_means$sample & methylation_means$ChIP == "ALL"]
  
}
meta$sample <- factor(x = meta$sample,levels = meta$sample[order(meta$age_mean_BP)])

p <- ggplot(data = meta[order(meta$age_mean_BP),],mapping = aes(x = sample,y = age_mean_BP,color = ))+
  geom_point()+
  geom_errorbar(aes(ymin=age_mean_BP-age_std_BP, ymax=age_mean_BP+age_std_BP), width=.1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))

ggsave(plot = p,filename = file.path(OUTPUT_FOLDER,paste("age_of_samples","png",sep = ".")))


saveRDS(object = meta,file = file.path(OUTPUT_FOLDER,paste("meta","rds",sep = ".")))

######

df_peaks <-  readRDS(file.path(INPUT_FOLDER, 
                               paste("H3K27ac","peak","mean","methylation",N_SAMPLES,"samples","merged",ANNOTATION,"rds",sep = ".")))

# histogramm peak length
p <- ggplot(data = df_peaks,mapping =  aes(x = end-start) )+ 
  geom_histogram(binwidth = 10)+
  theme_minimal()+
  ggtitle("histogramm of H3K27ac ChIP peak length in bp")+
  xlab("length in bp")
# + 
  # geom_vline(xintercept = 300, linetype="dotted", color = "red", size=1)
ggsave(plot = p,
       filename = file.path(OUTPUT_FOLDER,paste("hist_plot_peak_legth_H3K27ac","png",sep = ".")),
        width = 6, height = 4)


# histogramm CpGs in peak 
p <- ggplot(data = df_peaks,mapping =  aes(x = end-start,y = score) )+ 
  geom_point()+
  theme_minimal()+
  ggtitle("peak length vs CpGs contained")+
  xlab("length in bp")+
  ylab("# CpGs in peak bounderies ")


ggsave(plot = p,
       filename = file.path(OUTPUT_FOLDER,paste("peak_legth_H3K27ac_vs__nr_CpGs","png",sep = ".")),
       width = 6, height = 4)