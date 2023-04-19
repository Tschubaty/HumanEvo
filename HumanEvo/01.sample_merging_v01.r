#################################################################
##
##
##  input:
##
##            methylation_data
##
##
##  output:   methylation_data_chr_merged
##
##
##  v_01 19.04.2023
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
# library(Rtsne)
# library(ggrepel)
library("dplyr")
##################################### CONSTANTS ########################################
OUTPUT_FOLDER <- "methylation_data_chr_merged"
INPUT_FOLDER <- "methylation_data"
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

meta_from_Reich <- as.data.frame(readxl::read_xlsx(path = "AGDP.metadata.xlsx"))

meta_data <- meta_from_Reich[meta_from_Reich$`I-ID` %in% list.files(INPUT_FOLDER),]

# Rename multiple columns
meta_data <- meta_data %>% 
  dplyr::rename("sample" = "I-ID",
    "age mean (BP)" = "Date mean in BP [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]",
      "age std (BP)" = "Date standard deviation in BP [OxCal sigma for a direct radiocarbon date, and standard deviation of the uniform disribution between the two bounds for a contextual date]",
    "sex"  = "Genetic sex determination")


meta_data <- meta_data[,c("sample", "Coverage" ,"age mean (BP)", "age std (BP)","Group ID",                                                                                                                                                                                                                             
                           "Locality" ,                                                                                                                                                                                                                            
                         "Country" ,                                                                                                                                                                                                                             
                              "sex")]
#debug
# chr <- CHR_NAMES[22]
# file_nr <- 1

for(chr in CHR_NAMES){
  
  print(chr)

for(file_nr in 1:nrow(meta_data)){
  
  start.time <- Sys.time()
  print(file_nr)
  print(meta_data$sample[file_nr])
  
  if(file_nr == 1){
    
    
    bed_header <- read.delim(file = file.path(INPUT_FOLDER,meta_data$sample[file_nr],paste(meta_data$sample[file_nr],chr,"bed",sep = ".")))
    bed <- bed_header
    bed_header$score = NaN
    bed_header$name = NaN
    df <- data.frame(matrix(nrow = nrow(bed_header), ncol = nrow(meta_data)))
 
  }else{
        bed <- read.delim(file = file.path(INPUT_FOLDER,meta_data$sample[file_nr],paste(meta_data$sample[file_nr],chr,"bed",sep = ".")))
  }
  
  colnames(df)[file_nr] <- meta_data$sample[file_nr]
  df[,file_nr] <- bed$score
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}
  
  
  df <- cbind(bed_header,df)
  
  saveRDS(object = df,
          file = file.path(OUTPUT_FOLDER,
                           paste(chr,nrow(meta_data),"samples","rds",sep = ".")))
  print(paste(chr,"saved"))
}


