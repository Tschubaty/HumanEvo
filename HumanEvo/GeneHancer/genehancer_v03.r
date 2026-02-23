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


library(readxl)

# ---- Paths ----
file_S2 <- "C:/Users/Batyrev/Dropbox/R42v/Table S2 significant CpG [1b].xlsx"
file_S5 <- "C:/Users/Batyrev/Dropbox/R42v/Table S5 significant CpG [2b].xlsx"


fg1_tbl <- read_excel(file_S2, sheet = 1) |>
  dplyr::select(
    chrom,
    start,
    end,
    `H3K27ac peak coordinates`,
    GeneAnnotation,
    EnhancerAnnotation
  )


fg2_tbl <- read_excel(file_S5, sheet = 1) |>
  dplyr::select(
    chrom,
    start,
    end,
    name,
    GeneAnnotation,
    EnhancerAnnotation,
    GH_Score
  )




# ---- CRISPR tested genes ----
crispr_symbols <- read_excel(file_S5, sheet = 2) |>
  dplyr::pull("Functional validation of candidate genes using a pooled CRISPR knockout screen in human hepatocytes") |>
  as.character() |>
  trimws() |>
  (\(x) x[!is.na(x) & x != ""])() |>
  unique()


df_significant <- readRDS("C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/12.pipeline/45.pearson.CpG_permutation.min_delta.0.5/sim/df_significant.rds")


df_significant_Q <-df_significant[df_significant$pearson.q_emp < 0.1,]


##  Annotate df_significant_Q with COMPLETE GeneHancer gene database (no elite restriction) ----
##    Output FIRST in log format: annot_long (keeps ALL columns from BOTH tables, with sensible prefixes)

library(data.table)

# --- full GH gene table (no filtering) ---
gh_all <- as.data.table(gh_gene)    # full association table
ghc    <- as.data.table(gh_coords)  # coordinates table

# join gene associations to coordinates
gh_all_annot <- merge(
  gh_all,
  ghc,
  by = "GHid",
  all.x = TRUE,
  allow.cartesian = TRUE
)

# unify coord column naming
if ("chromosome" %in% names(gh_all_annot)) setnames(gh_all_annot, "chromosome", "chrom")
# keep start/end as is (should exist from gh_coords)

# --- prepare CpGs and GH tables with prefixes to keep ALL columns ---
cpg_dt <- as.data.table(df_significant_Q)
cpg_dt[, cpg_row := .I]

gh_dt  <- as.data.table(gh_all_annot)
gh_dt[, gh_row := .I]

setnames(cpg_dt, setdiff(names(cpg_dt), "cpg_row"), paste0("cpg_", setdiff(names(cpg_dt), "cpg_row")))
setnames(gh_dt,  setdiff(names(gh_dt),  "gh_row"),  paste0("gh_",  setdiff(names(gh_dt),  "gh_row")))


library(GenomicRanges)

# --- overlap (GRanges) ---
cpg_gr <- with(cpg_dt,
               GRanges(seqnames = cpg_chrom,
                       ranges   = IRanges(start = cpg_start, end = cpg_end)))

# --- fix: remove GH rows without coordinates (start/end NA) ---
cat("\nGH rows before dropping NA coords:", nrow(gh_dt), "\n")
cat("GH rows with NA start:", sum(is.na(gh_dt$gh_start)), "\n")
cat("GH rows with NA end  :", sum(is.na(gh_dt$gh_end)), "\n")

missing_coord_ghids <- unique(gh_dt[is.na(gh_start) | is.na(gh_end), gh_GHid])
cat("\nUnique GHids missing coordinates:", length(missing_coord_ghids), "\n")
print(head(missing_coord_ghids, 20))

gh_dt <- gh_dt[!is.na(gh_start) & !is.na(gh_end) & !is.na(gh_chrom)]

cat("GH rows after dropping NA coords :", nrow(gh_dt), "\n")

gh_gr <- with(gh_dt,
              GRanges(seqnames = gh_chrom,
                      ranges   = IRanges(start = gh_start, end = gh_end)))

ov <- findOverlaps(cpg_gr, gh_gr, ignore.strand = TRUE)

# --- annot_long: only hits (one row per CpG x GH association) ---
annot_long <- cbind(
  cpg_dt[queryHits(ov)],
  gh_dt[subjectHits(ov)]
)

# --- keep ALL CpGs: add "no hit" placeholder rows for missing CpG rows ---
hit_rows <- unique(queryHits(ov))
missing_rows <- setdiff(cpg_dt$cpg_row, hit_rows)

if (length(missing_rows) > 0) {
  # build a NA-filled gh_dt template with same columns as gh_dt
  gh_na <- gh_dt[0][, lapply(.SD, function(x) NA)]
  
  # one NA gh row per missing CpG row
  add_rows <- rbindlist(lapply(missing_rows, function(r) {
    cbind(cpg_dt[cpg_row == r], gh_na)
  }))
  
  annot_long <- rbindlist(list(annot_long, add_rows), use.names = TRUE, fill = TRUE)
}

# put key columns early (still keeps all columns)
key_first <- c(
  "cpg_row","gh_row",
  intersect(c("cpg_chrom","cpg_start","cpg_end","cpg_name","cpg_pearson.q_emp","cpg_pearson.p_emp","cpg_pearson.statistic"), names(annot_long)),
  intersect(c("gh_GHid","gh_symbol","gh_combined_score","gh_is_elite","gh_chrom","gh_start","gh_end"), names(annot_long))
)
setcolorder(annot_long, c(key_first, setdiff(names(annot_long), key_first)))

# elite-only (BUT keep placeholder CpGs with NA gh_is_elite)
annot_long_elite_2 <- annot_long[(is.na(gh_is_elite) | gh_is_elite == 1)]

# CpGs missing from elite hits
missing_rows_elite <- setdiff(cpg_dt$cpg_row, unique(annot_long_elite_2$cpg_row))

if (length(missing_rows_elite) > 0) {
  gh_na <- gh_dt[0][, lapply(.SD, function(x) NA)]   # NA template with gh_* columns
  add_rows_elite <- rbindlist(lapply(missing_rows_elite, function(r) {
    cbind(cpg_dt[cpg_row == r], gh_na)
  }))
  annot_long_elite_2 <- rbindlist(list(annot_long_elite_2, add_rows_elite),
                                  use.names = TRUE, fill = TRUE)
}

# sanity check: should now be 32
cat("\nElite CpG coverage check:\n")
cat("Unique CpGs in df_significant_Q:", nrow(unique(cpg_dt[, .(cpg_chrom, cpg_start, cpg_end)])), "\n")
cat("Unique CpGs in annot_long_elite_2:", nrow(unique(annot_long_elite_2[, .(cpg_chrom, cpg_start, cpg_end)])), "\n")

cat("\n==================== annot_long_elite_2 (first rows) ====================\n")
print(annot_long_elite_2[1:min(5, .N)])

cat("\n==================== annot_long_elite_2 structure =======================\n")
print(str(annot_long_elite_2))

cat("\n==================== unique GHids per CpG (top) =========================\n")
print(annot_long_elite_2[, .(
  n_pairs = .N,
  n_GHid  = uniqueN(gh_GHid, na.rm = TRUE),
  n_genes = uniqueN(gh_symbol, na.rm = TRUE)
), by = cpg_row][order(-n_GHid)][1:.N])




## ONE-SHOT SUMMARY for annot_long_elite_2:
## unique CpGs, unique enhancers (GHid), unique genes

library(dplyr)

annot_long_elite_2 %>%
  summarise(
    n_unique_CpGs      = n_distinct(cpg_chrom, cpg_start, cpg_end),
    n_unique_peaks     = n_distinct(cpg_name),
    n_unique_enhancers = n_distinct(gh_GHid),
    n_unique_genes     = n_distinct(gh_symbol)
  )


## ONE-SHOT SUMMARY for old fg2_tbl:
## unique CpGs, unique peaks, unique enhancers, unique genes

library(dplyr)

fg2_tbl %>%
  summarise(
    n_unique_CpGs      = n_distinct(chrom, start, end),
    n_unique_peaks     = n_distinct(name),
    n_unique_enhancers = n_distinct(EnhancerAnnotation),
    n_unique_genes     = n_distinct(GeneAnnotation)
  )
######################################################################################################



## CHECK: are NEW enhancers a subset of OLD enhancers?

library(dplyr)

# OLD enhancers (interaction-based)
enh_old <- fg2_tbl %>%
  pull(EnhancerAnnotation) %>%
  unique() %>%
  sort()

# NEW enhancers (strict regulatory-element based)
enh_new <- annot_long_elite_2 %>%
  pull(gh_GHid) %>%
  unique() %>%
  sort()

# enhancers in NEW but NOT in OLD (should be empty or very small)
enh_new_not_in_old <- setdiff(enh_new, enh_old)

list(
  n_enh_old = length(enh_old),
  n_enh_new = length(enh_new),
  n_new_not_in_old = length(enh_new_not_in_old),
  enh_new_not_in_old = enh_new_not_in_old
)


## EXACT gene name differences: annot_long_elite_2 vs old fg2_tbl

library(dplyr)

# genes in new (elite, strict) annotation
genes_new <- annot_long_elite_2 %>%
  pull(gh_symbol) %>%
  unique() %>%
  sort()

# genes in old fg2 table
genes_old <- fg2_tbl %>%
  pull(GeneAnnotation) %>%
  unique() %>%
  sort()

# NEW only
genes_only_in_new <- setdiff(genes_new, genes_old)

# OLD only
genes_only_in_old <- setdiff(genes_old, genes_new)

list(
  n_genes_new = length(genes_new),
  n_genes_old = length(genes_old),
  n_only_in_new = length(genes_only_in_new),
  n_only_in_old = length(genes_only_in_old),
  genes_only_in_new = genes_only_in_new,
  genes_only_in_old = genes_only_in_old
)



library(dplyr)
library(stringr)

# helper: list GHids supporting each gene in NEW (elite-only strict overlaps)
new_gene_to_gh <- annot_long_elite_2 %>%
  distinct(gh_symbol, gh_GHid) %>%
  group_by(gh_symbol) %>%
  summarise(new_GHids = paste(sort(unique(gh_GHid)), collapse=";"), .groups="drop")

# helper: list GHids supporting each gene in OLD (verified elite pairs only)
old_gene_to_gh <- fg2_elite_check %>%
  filter(is_elite_old == 1) %>%
  distinct(GeneAnnotation, EnhancerAnnotation) %>%
  group_by(GeneAnnotation) %>%
  summarise(old_GHids = paste(sort(unique(EnhancerAnnotation)), collapse=";"), .groups="drop") %>%
  rename(gh_symbol = GeneAnnotation)

# gene-level comparison table
gene_cmp <- full_join(new_gene_to_gh, old_gene_to_gh, by="gh_symbol") %>%
  mutate(
    in_new = !is.na(new_GHids),
    in_old_elite = !is.na(old_GHids)
  )

# show ONLY-NEW and ONLY-OLD (elite verified) with supporting GHids
list(
  only_in_new_elite = gene_cmp %>% filter(in_new & !in_old_elite) %>% arrange(gh_symbol),
  only_in_old_elite = gene_cmp %>% filter(!in_new & in_old_elite) %>% arrange(gh_symbol)
)
