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

# ---- CRISPR tested genes ----
crispr_symbols <- read_excel(file_S5, sheet = 2) |>
  dplyr::pull("Functional validation of candidate genes using a pooled CRISPR knockout screen in human hepatocytes") |>
  as.character() |>
  trimws() |>
  (\(x) x[!is.na(x) & x != ""])() |>
  unique()



gh <- toupper(unique(gh_gene$symbol))
bg <- toupper(unique(background_genes))
cr <- toupper(unique(crispr_symbols))
length(cr)
ix<- intersect(toupper(gh ), toupper(cr))
length(ix)

#########################################################################
library(readxl)
library(dplyr)

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



gh_symbols_upper <- toupper(unique(gh_gene$symbol))
fg1_not_in_gh <- fg1_tbl |>
  dplyr::mutate(GeneAnnotation_upper = toupper(GeneAnnotation)) |>
  dplyr::filter(!GeneAnnotation_upper %in% gh_symbols_upper)

fg1_not_in_gh
nrow(fg1_not_in_gh)

fg2_not_in_gh <- fg2_tbl |>
  dplyr::mutate(GeneAnnotation_upper = toupper(GeneAnnotation)) |>
  dplyr::filter(!GeneAnnotation_upper %in% gh_symbols_upper)

fg2_not_in_gh
nrow(fg2_not_in_gh)




#fg2_not_in_gh$

################################################################################

# ---- Analysis 1 foreground ----
fg1_symbols <- unique(fg1_tbl$GeneAnnotation)


length(fg1_symbols)
ix<- intersect(toupper(gh ), toupper(fg1_symbols))
length(ix)

fg2_tbl <- read_excel(file_S5, sheet = 1)

# ---- Analysis 2 foreground ----
fg2_symbols <- unique(fg2_tbl$GeneAnnotation)


length(fg2_symbols)
ix<- intersect(toupper(gh ), toupper(fg2_symbols))
length(ix)


not_in_gh <- function(symbols, gh_symbols) {
  data.frame(
    symbol = symbols[!toupper(symbols) %in% toupper(gh_symbols)],
    stringsAsFactors = FALSE
  ) |> 
    dplyr::distinct(symbol)
}
crispr_not_in_gh <- not_in_gh(crispr_symbols, gh)

crispr_not_in_gh
nrow(crispr_not_in_gh)


fg2_not_in_gh_tbl <- fg2_tbl |>
  dplyr::mutate(GeneAnnotation_upper = toupper(GeneAnnotation)) |>
  dplyr::filter(!GeneAnnotation_upper %in% gh)

fg2_not_in_gh_tbl
nrow(fg2_not_in_gh_tbl)
unique(fg2_not_in_gh_tbl$GeneAnnotation_upper)
# 
# fg1_not_in_gh_tbl <- fg1_not_in_gh_tbl |>
#   dplyr::select(-GeneAnnotation_upper)

fg2_not_in_gh_tbl <- fg2_not_in_gh_tbl |>
  dplyr::select(-GeneAnnotation_upper)



# ghid_fg1_not_in_gh <- unique(fg1_not_in_gh$EnhancerAnnotation)
ghid_fg2_not_in_gh <- unique(fg2_not_in_gh$EnhancerAnnotation)

# gh_gene_fg1 <- gh_gene[GHid %in% ghid_fg1_not_in_gh]
# gh_gene_fg1


gh_gene_fg2 <- gh_gene[GHid %in% ghid_fg2_not_in_gh]
unique(gh_gene_fg2$symbol)
unique(fg2_not_in_gh_tbl$GeneAnnotation)


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


########################## conversion in the Entrez universe ###################

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

to_entrez <- function(x) {
  x <- unique(trimws(x))
  x <- x[x != "" & !is.na(x)]
  
  out <- data.frame(input = x, ENTREZID = NA_character_, keytype = NA_character_, stringsAsFactors = FALSE)
  
  # 1) If ENSG IDs are present, map them directly
  is_ensg <- grepl("^ENSG\\d+", x)
  if (any(is_ensg)) {
    m <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = unique(x[is_ensg]),
                               keytype = "ENSEMBL",
                               columns = "ENTREZID")
    m <- m[!is.na(m$ENTREZID), c("ENSEMBL","ENTREZID")]
    out$ENTREZID[match(m$ENSEMBL, out$input)] <- as.character(m$ENTREZID)
    out$keytype[match(m$ENSEMBL, out$input)] <- "ENSEMBL"
  }
  
  # 2) Try HGNC SYMBOL for everything still unmapped
  unm <- is.na(out$ENTREZID)
  if (any(unm)) {
    m <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = unique(out$input[unm]),
                               keytype = "SYMBOL",
                               columns = "ENTREZID")
    m <- m[!is.na(m$ENTREZID), c("SYMBOL","ENTREZID")]
    out$ENTREZID[match(m$SYMBOL, out$input)] <- as.character(m$ENTREZID)
    out$keytype[match(m$SYMBOL, out$input)] <- "SYMBOL"
  }
  
  # 3) Try ALIAS for everything still unmapped (catches many old names)
  unm <- is.na(out$ENTREZID)
  if (any(unm)) {
    m <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = unique(out$input[unm]),
                               keytype = "ALIAS",
                               columns = "ENTREZID")
    m <- m[!is.na(m$ENTREZID), c("ALIAS","ENTREZID")]
    out$ENTREZID[match(m$ALIAS, out$input)] <- as.character(m$ENTREZID)
    out$keytype[match(m$ALIAS, out$input)] <- "ALIAS"
  }
  
  out
}

# length(unique(gh_gene$symbol))
# gh_gene_symbol_Entrez<-  to_entrez(unique(gh_gene$symbol))
# gh_valid_Entrez <- unique(gh_gene_symbol_Entrez$ENTREZID[!is.na(gh_gene_symbol_Entrez$ENTREZID)])
# gh_valid_Entrez <- sort(as.numeric(gh_valid_Entrez))
# plot(gh_valid_Entrez)
# length(gh_valid_Entrez)
# 
# length(unique(all_double_elite_genes))
# all_double_elite_genes_Entrez <- to_entrez(unique(all_double_elite_genes))
# all_double_elite_genes_valid_Entrez <- unique(all_double_elite_genes_Entrez$ENTREZID[!is.na(all_double_elite_genes_Entrez$ENTREZID)])
# all_double_elite_genes_valid_Entrez <-  sort(as.numeric(all_double_elite_genes_valid_Entrez ))
# plot(all_double_elite_genes_valid_Entrez)
# length(all_double_elite_genes_valid_Entrez)
# 
# length(unique(background_genes)) 
# background_genes_Entrez <- to_entrez(unique(background_genes))
# background_genes_Entrez <- unique(background_genes_Entrez$ENTREZID[!is.na(background_genes_Entrez$ENTREZID)])
# background_genes_valid_Entrez <-  sort(as.numeric(background_genes_Entrez))
# plot(background_genes_valid_Entrez)
# length(background_genes_valid_Entrez)



## =============================================================================
## Summary + plots for Entrez mapping quality across 3 GeneHancer-derived sets
##   - gh_gene$symbol (all)
##   - all_double_elite_genes (all elite GeneHancer genes)
##   - background_genes (elite enhancers that overlap CpGs in H3K27ac peaks)
##
## Output:
##   - informative cat() summaries
##   - stacked barplot: mapped Entrez IDs split by < 1e6 vs >= 1e6
##
## NOTE:
##   Entrez Gene IDs for human genes are typically << 1,000,000. Values >= 1e6
##   often indicate mapping artifacts (e.g., transcript IDs, non-gene entities,
##   or unexpected identifier types) OR non-standard IDs if they slipped in.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- helper: summarize one set ----
summarize_entrez_set <- function(name, symbols_vec) {
  symbols_vec <- unique(symbols_vec)
  symbols_vec <- symbols_vec[!is.na(symbols_vec) & symbols_vec != ""]
  
  mapped <- to_entrez(symbols_vec)
  
  n_input <- length(symbols_vec)
  n_mapped_rows <- sum(!is.na(mapped$ENTREZID))
  
  entrez_all <- mapped$ENTREZID[!is.na(mapped$ENTREZID)]
  entrez_all_num <- suppressWarnings(as.numeric(entrez_all))
  entrez_all_num <- entrez_all_num[!is.na(entrez_all_num)]
  
  entrez_unique <- unique(entrez_all_num)
  n_unique <- length(entrez_unique)
  
  n_lt1m <- sum(entrez_unique < 1e6)
  n_ge1m <- sum(entrez_unique >= 1e6)
  n_non_numeric_dropped <- length(entrez_all) - length(entrez_all_num)
  
  # also count how many inputs remained unmapped
  n_unmapped_inputs <- sum(is.na(mapped$ENTREZID))
  
  # for context: how many mappings came from each keytype
  keytype_tab <- table(mapped$keytype, useNA = "ifany")
  
  # cats
  cat("\n============================================================\n")
  cat("Set:", name, "\n")
  cat("------------------------------------------------------------\n")
  cat("Unique input symbols:", n_input, "\n")
  cat("Unmapped inputs (no Entrez):", n_unmapped_inputs, sprintf("(%.2f%%)\n", 100 * n_unmapped_inputs / n_input))
  cat("Mapped rows (incl. 1-to-many):", n_mapped_rows, "\n")
  cat("Unique numeric Entrez IDs:", n_unique, "\n")
  cat("  - Entrez < 1e6:", n_lt1m, sprintf("(%.2f%% of mapped unique)\n", 100 * n_lt1m / max(1, n_unique)))
  cat("  - Entrez ≥ 1e6:", n_ge1m, sprintf("(%.2f%% of mapped unique)\n", 100 * n_ge1m / max(1, n_unique)))
  cat("Non-numeric Entrez values dropped after coercion:", n_non_numeric_dropped, "\n")
  cat("Keytype usage (where mapping succeeded):\n")
  print(keytype_tab)
  
  # return a small summary row for plotting
  data.table(
    set = name,
    mapped_unique_total = n_unique,
    mapped_lt_1e6 = n_lt1m,
    mapped_ge_1e6 = n_ge1m,
    unmapped_inputs = n_unmapped_inputs,
    input_symbols = n_input
  )
}

# ---- run summaries ----
res1 <- summarize_entrez_set("GeneHancer all symbols", unique(gh_gene$symbol))
res2 <- summarize_entrez_set("GeneHancer double elite", unique(all_double_elite_genes))
res3 <- summarize_entrez_set("Elite and H3K27ac CpG-overlap", unique(background_genes))



summary_dt <- rbindlist(list(res1, res2, res3), use.names = TRUE, fill = TRUE)

cat("\n============================================================\n")
cat("Combined summary table:\n")
print(summary_dt)

# # ---- stacked barplot (mapped unique Entrez split by threshold) ----
# # prepare matrix for base R barplot
# mat <- rbind(
#   "< 1e6" = summary_dt$mapped_lt_1e6,
#   ">= 1e6" = summary_dt$mapped_ge_1e6
# )
# colnames(mat) <- summary_dt$set
# 
# par(mar = c(10, 5, 4, 1))  # more bottom margin for rotated labels
# bp <- barplot(
#   mat,
#   beside = FALSE,
#   las = 2,              # rotate labels
#   ylab = "Unique mapped Entrez IDs",
#   main = "Entrez mapping distribution across GeneHancer-derived sets"
# )
# 
# legend(
#   "topleft",
#   legend = rownames(mat),
#   bty = "n",
#   inset = 0.01
# )
# 
# # add totals on top of bars
# totals <- colSums(mat)
# text(x = bp, y = totals, labels = totals, pos = 3, cex = 0.9)
# 
# # optional: also print unmapped fraction as text under each bar (comment out if you don’t want)
# unmapped_pct <- round(100 * summary_dt$unmapped_inputs / pmax(1, summary_dt$input_symbols), 1)
# mtext(
#   paste0("unmapped: ", unmapped_pct, "%"),
#   side = 1, at = bp, line = 6, cex = 0.8
# )
# 
# par(mar = c(5, 4, 4, 2) + 0.1)


## =============================================================================
## NICE PLOT: per set, show (1) TOTAL input symbols and (2) MAPPED-only,
## both as stacked bars (<1e6 vs >=1e6, plus unmapped for the TOTAL bar)
## =============================================================================
pdf("Entrez_mapping_summary.pdf", width = 18, height = 9)

# summary_dt must exist and contain:
# set, input_symbols, unmapped_inputs, mapped_lt_1e6, mapped_ge_1e6, mapped_unique_total

# ---- derive counts for plotting ----
plot_dt <- summary_dt[, .(
  set,
  input_symbols,
  unmapped_inputs,
  mapped_lt_1e6,
  mapped_ge_1e6
)]

plot_dt[, mapped_total := mapped_lt_1e6 + mapped_ge_1e6]
plot_dt[, unmapped_unique := pmax(0, input_symbols - mapped_total)]  # display-friendly

# ---- build matrices for barplot ----
# left = TOTAL symbols (mapped split + unmapped)
mat_total <- rbind(
  "Mapped Entrez IDs"  = plot_dt$mapped_lt_1e6,
  "provisional Entrez IDs" = plot_dt$mapped_ge_1e6,
  "Unmapped"       = plot_dt$unmapped_unique
)

# right = MAPPED-only (mapped split only)
mat_mapped <- rbind(
  "Mapped Entrez IDs"  = plot_dt$mapped_lt_1e6,
  "provisional Entrez IDs" = plot_dt$mapped_ge_1e6
)

colnames(mat_total)  <- plot_dt$set
colnames(mat_mapped) <- plot_dt$set

# ---- nicer colors ----
cols_total  <- c("Mapped Entrez IDs"="#3B82F6", "provisional Entrez IDs"="#F59E0B", "Unmapped"="#D1D5DB")
cols_mapped <- c("Mapped Entrez IDs"="#3B82F6", "provisional Entrez IDs"="#F59E0B")

# ---- plot layout: two panels side-by-side ----
op <- par(no.readonly = TRUE)
par(mfrow = c(1,2), mar = c(9,5,4,1))

# ========== Panel 1: TOTAL ==========

# totals on top
tot_total <- colSums(mat_total)

bp1 <- barplot(
  mat_total,
  beside = FALSE,
  col = cols_total[rownames(mat_total)],
  border = NA,
  names.arg = rep("", ncol(mat_total)),
  ylab = "Count",
  main = "TOTAL input symbols\n(mapped + unmapped)",
  ylim = c(0, max(tot_total) * 1.08) 
)
axis(
  side = 1,
  at = bp1,
  labels = c(
    "GeneHancer\n(all symbols)",
    "Double-elite\nGeneHancer",
    "H3K27ac\nCpG overlap"
  ),
  tick = FALSE,
  line = 0
)



ymax1 <- max(tot_total)
text(
  x = bp1,
  y = tot_total + 0.03 * ymax1,
  labels = format(tot_total, big.mark = ","),
  cex = 0.9
)


# mapped % under x labels
mapped_pct <- round(100 * plot_dt$mapped_total / pmax(1, plot_dt$input_symbols), 1)
mtext(paste0("mapped: ", mapped_pct, "%"),
      side = 1, at = bp1, line = 6, cex = 0.85)

legend(
  "topright",
  inset = c(-0.22, 0.02),
  legend = rownames(mat_total),
  fill = cols_total[rownames(mat_total)],
  bty = "n",
  cex = 0.9
)

# ========== Panel 2: MAPPED-only ==========
tot_mapped <- colSums(mat_mapped)
bp2 <- barplot(
  mat_mapped,
  beside = FALSE,
  col = cols_mapped[rownames(mat_mapped)],
  border = NA,
  names.arg = rep("", ncol(mat_total)),
  ylab = "Count",
  main = "MAPPED-only symbols\n(Entrez IDs)",
  ylim = c(0, max(tot_mapped) * 1.10)
)


axis(
  side = 1,
  at = bp2,
  labels = c(
    "GeneHancer\n(all symbols)",
    "Double-elite\nGeneHancer",
    "H3K27ac\nCpG overlap"
  ),
  tick = FALSE,
  line = 0
)


# mapped totals on top

ymax2 <- max(tot_mapped)
text(
  x = bp2,
  y = tot_mapped + 0.03 * ymax2,
  labels = format(tot_mapped, big.mark = ","),
  cex = 0.9
)


# # >=1e6 fraction under labels (your “strange IDs” proxy)
# ge1m_pct <- round(100 * plot_dt$mapped_ge_1e6 / pmax(1, plot_dt$mapped_total), 2)
# mtext(paste0(">=1e6: ", ge1m_pct, "%"),
#       side = 1, at = bp2, line = 6, cex = 0.85)

legend(
  "topright",
  inset = c(-0.22, 0.02),
  legend = rownames(mat_mapped),
  fill = cols_mapped[rownames(mat_mapped)],
  bty = "n",
  cex = 0.9
)

par(op)
dev.off()


# =============================================================================
# SAVE Entrez mapping outputs (3 sets) for reuse
# =============================================================================

out_map_dir <- file.path(this.dir, "Entrez_sets_GeneHancer")
dir.create(out_map_dir, showWarnings = FALSE, recursive = TRUE)

# --- helper: save one set cleanly ---
save_entrez_set <- function(prefix, mapping_df) {
  # mapping_df must have columns: input, ENTREZID, keytype
  
  # full mapping (keeps unmapped too)
  saveRDS(mapping_df, file.path(out_map_dir, paste0(prefix, "_mapping_input_to_entrez.rds")))
  fwrite(mapping_df, file.path(out_map_dir, paste0(prefix, "_mapping_input_to_entrez.tsv")), sep = "\t")
  
  # unique numeric Entrez IDs
  entrez_vec <- mapping_df$ENTREZID[!is.na(mapping_df$ENTREZID)]
  entrez_num <- suppressWarnings(as.numeric(entrez_vec))
  entrez_num <- entrez_num[!is.na(entrez_num)]
  entrez_unique <- sort(unique(entrez_num))
  
  saveRDS(entrez_unique, file.path(out_map_dir, paste0(prefix, "_entrez_unique_numeric.rds")))
  write.table(entrez_unique,
              file.path(out_map_dir, paste0(prefix, "_entrez_unique_numeric.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # also save the unmapped inputs (so you can inspect later)
  unmapped_inputs <- mapping_df$input[is.na(mapping_df$ENTREZID)]
  unmapped_inputs <- sort(unique(unmapped_inputs))
  saveRDS(unmapped_inputs, file.path(out_map_dir, paste0(prefix, "_unmapped_inputs.rds")))
  write.table(unmapped_inputs,
              file.path(out_map_dir, paste0(prefix, "_unmapped_inputs.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  invisible(list(entrez_unique = entrez_unique, unmapped_inputs = unmapped_inputs))
}

# =============================================================================
# Build the 3 mapping tables (if not already in memory)
# =============================================================================

gh_gene_symbol_Entrez <- to_entrez(unique(gh_gene$symbol))
all_double_elite_genes_Entrez <- to_entrez(unique(all_double_elite_genes))
background_genes_Entrez <- to_entrez(unique(background_genes))

# =============================================================================
# Save all 3
# =============================================================================

save_entrez_set("01_GeneHancer_all_symbols", gh_gene_symbol_Entrez)
save_entrez_set("02_GeneHancer_double_elite", all_double_elite_genes_Entrez)
save_entrez_set("03_H3K27ac_CpG_overlap_background", background_genes_Entrez)

# =============================================================================
# Write a tiny README summary
# =============================================================================
readme_lines <- c(
  "Entrez mapping outputs (GeneHancer-derived sets)",
  paste("Created:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Files per set:",
  "  *_mapping_input_to_entrez.rds/.tsv  : full mapping table (input symbol -> EntrezID, includes unmapped)",
  "  *_entrez_unique_numeric.rds/.txt    : unique numeric Entrez IDs (sorted)",
  "  *_unmapped_inputs.rds/.txt          : input symbols that could not be mapped",
  "",
  "Set definitions:",
  "  01_GeneHancer_all_symbols                : unique(gh_gene$symbol)",
  "  02_GeneHancer_double_elite               : unique(all_double_elite_genes) (is_elite==1)",
  "  03_H3K27ac_CpG_overlap_background         : unique(background_genes) (elite enhancers overlapping CpGs in H3K27ac peaks)"
)

writeLines(readme_lines, con = file.path(out_map_dir, "README_Entrez_sets.txt"))

cat("\nSaved Entrez sets to:\n", out_map_dir, "\n")
