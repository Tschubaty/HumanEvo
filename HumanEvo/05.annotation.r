#################################################################
##
##
##  input:
##           "04.methylation_vs_age_CpG"
##
##
##  output:    "05.annotation"
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
OUTPUT_FOLDER <- "05.annotation"
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

df <- readRDS(file = file.path(INPUT_FOLDER,
                                     paste("H3K27ac",
                                           "with",
                                           "pearson",
                                           N_SAMPLES,
                                           "samples",
                                           ANNOTATION,
                                           "rds",
                                           sep = ".")))


alpha <- 0.001

df_sig <- df[df$pearson_p_val <= 0.001,]
df_sig <- df_sig[!is.na(df_sig$chrom),]

length(unique(df_sig$name))

nrow(df_sig)

print(sum(df$pearson_p_val < 0.001,na.rm = T))

print(sum(df$radom_p_val < 0.001,na.rm = T))


bed3 <- df_sig[!is.na(df_sig$chrom),c(1:3)]


df_sig

write.table(x = bed3[bed3$chrom %in% CHR_NAMES[1:11],],
            file = file.path(OUTPUT_FOLDER,paste("chrom01-11","results","p",alpha,"bed",sep = ".")),
            quote = FALSE,row.names = FALSE,col.names = FALSE)

write.table(x = bed3[bed3$chrom %in% CHR_NAMES[12:22],],
            file = file.path(OUTPUT_FOLDER,paste("chrom12-22","results","p",alpha,"bed",sep = ".")),
            quote = FALSE,row.names = FALSE,col.names = FALSE)

### make Genhance query with files https://genome.ucsc.edu/cgi-bin/hgTables


df1 <- read.delim(file =  file.path(OUTPUT_FOLDER,paste("chrom01-11","GH-Interactions","p",alpha,"bed",sep = ".")),
                  header = FALSE,
                  skip = 1,
                  col.names = colnames(df)[1:5])
                    
df2 <- read.delim(file =  file.path(OUTPUT_FOLDER,paste("chrom12-22","GH-Interactions","p",alpha,"bed",sep = ".")),
                  header = FALSE,
                  skip = 1,
                  col.names = colnames(df)[1:5])

df_GH <- rbind(df1,df2)

df_sig$annotation_narrow <- NA
df_sig$annotation_common <- NA
df_sig$annotation_narrow_GH <-NA

for(r in 1:nrow(df_sig)){
  
  
  
  enhancer_regions <- df_GH[df_GH$chrom == df_sig$chrom[r] & 
                              df_GH$start <=  df_sig$start[r] & 
                              df_sig$end[r] <= df_GH$end,]
  
  if(nrow(enhancer_regions) > 0){
  # most narrow 
  smallest_enhancer_region <- enhancer_regions[which.min(enhancer_regions$end - enhancer_regions$start),]
  most_narrow <- unlist(strsplit(smallest_enhancer_region$name,split='/', fixed=TRUE) ) 

  
  # most common gene mentioning 
  gene_names <- unlist(strsplit(enhancer_regions$name,split='/', fixed=TRUE) )
  most_common <- names(sort(summary(as.factor(gene_names)), decreasing=T)[1])
   
  
  df_sig$annotation_narrow[r] <- most_narrow[1] 
  df_sig$annotation_common[r] <- most_common
  df_sig$annotation_narrow_GH[r] <- most_narrow[2] 
  }
}

nrow(df_sig)

length(unique(df_sig$annotation_narrow_GH))

length(unique(df_sig$annotation_narrow))

length(unique(df_sig$annotation_common))

sort(table(df_sig$annotation_narrow_GH),decreasing = T)

sort(table(df_sig$annotation_narrow),decreasing = T)

genes_for_gene_organizer <- sort(table(df_sig$annotation_narrow),decreasing = T)

text_out <- names(genes_for_gene_organizer[genes_for_gene_organizer > 1])

write.csv(x = text_out,
          file = file.path(OUTPUT_FOLDER,
                           paste("Genes","narrow","GH-Interactions","p",alpha,"txt",sep = ".")),
          quote = FALSE,
          row.names = FALSE
          )


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

#########################################

