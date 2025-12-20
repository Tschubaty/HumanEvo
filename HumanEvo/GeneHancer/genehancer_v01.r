# Documentation -----------------------------------------------------------
##
##
##  input: 
##
##
##  output: 
## 
##  
##  
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------

#Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
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
# Load Libraries ----------------------------------------------------------
library(data.table)
library(GenomicRanges)
library(dplyr)

H3K27ac.only.45.samples.histogram.merged2025.hg19 <- readRDS("C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/H3K27ac.only.45.samples.histogram.merged2025.hg19.rds")

## 1) GeneHancer-Dateien einlesen
gh_coords <- fread("GeneHancer_AnnotSV_hg19_v5.25.txt")
gh_gene   <- fread("GeneHancer_AnnotSV_gene_association_scores_v5.25.txt")

## 2) (Double-)Elite-Enhancer-Gene-Paare filtern
##    -> hier nehme ich is_elite == 1 als "double elite" Interaktion
gh_elite <- gh_gene[is_elite == 1]

## 3) Koordinaten an die Elite-Paare anhängen
gh_elite_annot <- gh_elite %>%
  inner_join(gh_coords, by = "GHid") %>%
  rename(chrom = chromosome)

## 4) GRanges-Objekte bauen

# CpGs in H3K27ac-Peaks
cpg_gr <- with(H3K27ac.only.45.samples.histogram.merged2025.hg19,
               GRanges(seqnames = chrom,
                       ranges   = IRanges(start = start, end = end)))

# Elite-GeneHancer-Elemente
gh_gr <- with(gh_elite_annot,
              GRanges(seqnames = chrom,
                      ranges   = IRanges(start = start, end = end),
                      GHid     = GHid,
                      symbol   = symbol,
                      combined_score = combined_score,
                      is_elite = is_elite))

## 5) Overlap: welche (double-)Elite-GeneHancer-Elemente enthalten CpGs?

ov <- findOverlaps(cpg_gr, gh_gr, ignore.strand = TRUE)

hits_idx <- subjectHits(ov)

hits_df <- data.frame(
  GHid     = mcols(gh_gr)$GHid[hits_idx],
  symbol   = mcols(gh_gr)$symbol[hits_idx],
  combined_score = mcols(gh_gr)$combined_score[hits_idx],
  is_elite = mcols(gh_gr)$is_elite[hits_idx],
  chrom    = as.character(seqnames(gh_gr))[hits_idx],
  start    = start(gh_gr)[hits_idx],
  end      = end(gh_gr)[hits_idx],
  stringsAsFactors = FALSE
)

## 6) Einzigartige (GHid, Gen)-Paare mit Koordinaten = Hintergrund-Enrichment-Set

double_elite_enhancers_bg <- hits_df %>%
  distinct(GHid, symbol, combined_score, is_elite, chrom, start, end) %>%
  arrange(symbol, GHid, chrom, start)

## 7) Hintergrund-Gene-Liste
background_genes <- sort(unique(double_elite_enhancers_bg$symbol))

## Optional speichern
saveRDS(double_elite_enhancers_bg,
        "double_elite_enhancers_background_H3K27ac.hg19.rds")

write.table(background_genes,
            "background_genes_double_elite_H3K27ac.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)



## --------------------------------------------------------------------
## APPEND AFTER THE EXISTING CODE (uses background_genes already created)
## --------------------------------------------------------------------

library(dplyr)

## 1) Count all possible double-elite GeneHancer genes (global universe)
all_double_elite_genes <- gh_elite_annot %>%
  distinct(symbol) %>%
  pull(symbol)

length_all_double_elite_genes <- length(all_double_elite_genes)


## 2) Count restricted set = genes hit by double-elite enhancers
##    that also overlap CpGs inside your H3K27ac peaks
length_restricted_genes <- length(background_genes)


## 3) Calculate correction factor
fraction_retained <- length_restricted_genes / length_all_double_elite_genes
fraction_removed  <- 1 - fraction_retained


## 4) Output summary table
correction_summary <- data.frame(
  metric = c("all_double_elite_genes",
             "restricted_genes_in_H3K27ac",
             "fraction_retained",
             "fraction_removed"),
  value  = c(length_all_double_elite_genes,
             length_restricted_genes,
             fraction_retained,
             fraction_removed)
)
options(scipen = 999)   # turn off scientific notation globally
print(correction_summary)

