# Documentation -----------------------------------------------------------
##
##  Purpose:
##    Sanity-check UCSC Table Browser GeneHancer exports against our significant CpGs.
##    We quantify:
##      (1) how many unique CpGs are in df_significant_Q_annot.bed
##      (2) how many of those CpGs overlap at least one *interaction region*
##      (3) how many unique GeneHancer IDs and gene symbols exist in the interaction export
##      (4) whether Regulatory-Element GHids are a subset of Interaction GHids
##      (5) after restricting interactions to GHids present in Regulatory-Elements,
##          how many GHids and gene symbols remain
##
##  Inputs (relative to this.dir):
##    - df_significant_Q_annot.bed
##    - TableBrowserOutPut/S5geneHancerInteractionsDoubleElite.bed
##    - TableBrowserOutPut/S5geneHancerRegElemetsDoubleElite.bed
##
##  Notes:
##    - Interaction export format (BED-like, 5 columns observed):
##        chr, start, end, "GENE/GHID", score
##      where column4 encodes both gene symbol + GHid.
##    - RegElements export format (BED9 observed):
##        chr, start, end, GHID, score, strand, thickStart, thickEnd, itemRgb
##
##  Author: Daniel Batyrev (HUJI 777634015)
##
# ------------------------------------------------------------------------
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
cpg_bed_path <- file.path(this.dir, "df_significant_Q_annot.bed")

tb_dir <- file.path(this.dir, "TableBrowserOutPut")
interaction_bed_path <- file.path(tb_dir, "S5geneHancerInteractionsDoubleElite.bed")
regelem_bed_path     <- file.path(tb_dir, "S5geneHancerRegElemetsDoubleElite.bed")

# ----------------------------- Read CpGs (BED3) --------------------------
# BED is 0-based start, half-open end by convention. We treat them as given and only use overlaps consistently.
cpg_dt <- fread(cpg_bed_path, header = FALSE)
stopifnot(ncol(cpg_dt) >= 3)
setnames(cpg_dt, 1:3, c("chrom", "start", "end"))
cpg_dt <- unique(cpg_dt[, .(chrom = as.character(chrom),
                            start = as.integer(start),
                            end   = as.integer(end))])

n_cpg_unique <- nrow(cpg_dt)

# ----------------------------- Read Interaction export -------------------
int_dt <- fread(interaction_bed_path, header = FALSE)
stopifnot(ncol(int_dt) >= 4)
# expected: V1 chr, V2 start, V3 end, V4 "GENE/GHID", V5 score (optional)
setnames(int_dt, 1:4, c("chrom", "start", "end", "gene_ghid"))
if (ncol(int_dt) >= 5) setnames(int_dt, 5, "score")

int_dt[, chrom := as.character(chrom)]
int_dt[, start := as.integer(start)]
int_dt[, end   := as.integer(end)]
int_dt[, gene_ghid := as.character(gene_ghid)]

# Parse gene symbol + GHid from "GENE/GHID"
int_dt[, c("gene_symbol", "GHid") := tstrsplit(gene_ghid, "/", fixed = TRUE, keep = c(1,2))]
int_dt[, gene_symbol := trimws(gene_symbol)]
int_dt[, GHid := trimws(GHid)]

# Basic counts in interaction file
n_int_unique_GHid  <- uniqueN(int_dt$GHid)
n_int_unique_genes <- uniqueN(int_dt$gene_symbol)

# ----------------------------- CpG overlap with interactions -------------
cpg_gr <- GRanges(seqnames = cpg_dt$chrom,
                  ranges   = IRanges(start = cpg_dt$start, end = cpg_dt$end))

int_gr <- GRanges(seqnames = int_dt$chrom,
                  ranges   = IRanges(start = int_dt$start, end = int_dt$end))

ov_cpg_int <- findOverlaps(cpg_gr, int_gr, ignore.strand = TRUE)
n_cpg_with_at_least_one_interaction <- length(unique(queryHits(ov_cpg_int)))

# ----------------------------- Read Regulatory Elements export -----------
re_dt <- fread(regelem_bed_path, header = FALSE)
stopifnot(ncol(re_dt) >= 4)
# expected BED9: chr, start, end, GHid, score, strand, thickStart, thickEnd, itemRgb
setnames(re_dt, 1:4, c("chrom", "start", "end", "GHid"))
re_dt[, chrom := as.character(chrom)]
re_dt[, start := as.integer(start)]
re_dt[, end   := as.integer(end)]
re_dt[, GHid  := trimws(as.character(GHid))]

re_dt <- unique(re_dt[, .(chrom, start, end, GHid)])  # remove exact duplicates from UCSC export
n_re_unique_GHid <- uniqueN(re_dt$GHid)

# ----------------------------- Subset check: RegElem GHids ⊆ Interaction GHids
re_not_in_int <- setdiff(unique(re_dt$GHid), unique(int_dt$GHid))

# ----------------------------- Restrict interactions to RegElem GHids ----
int_dt_restricted <- int_dt[GHid %in% unique(re_dt$GHid)]
n_int_restricted_unique_GHid  <- uniqueN(int_dt_restricted$GHid)
n_int_restricted_unique_genes <- uniqueN(int_dt_restricted$gene_symbol)

# ----------------------------- Print summary (analysis log) --------------
cat("\n======================= INPUT COUNTS =======================\n")
cat("Unique CpGs in df_significant_Q_annot.bed:", n_cpg_unique, "\n")

cat("\n==================== INTERACTION EXPORT =====================\n")
cat("CpGs with >= 1 overlap to an interaction region:", n_cpg_with_at_least_one_interaction, "of", n_cpg_unique, "\n")
cat("Unique GHids in interaction file:", n_int_unique_GHid, "\n")
cat("Unique gene symbols in interaction file:", n_int_unique_genes, "\n")

cat("\n=================== REG-ELEMENT EXPORT ======================\n")
cat("Unique GHids in RegElements file:", n_re_unique_GHid, "\n")

cat("\n========= SUBSET CHECK: RegElements GHids in Interactions ========\n")
if (length(re_not_in_int) == 0) {
  cat("OK: RegElements GHids are a subset of Interaction GHids (no missing GHids).\n")
} else {
  cat("WARNING: These RegElements GHids are NOT present in the Interaction file:\n")
  print(sort(re_not_in_int))
}

cat("\n====== INTERACTIONS RESTRICTED TO RegElements GHids ======\n")
cat("Unique GHids left after restriction:", n_int_restricted_unique_GHid, "\n")
cat("Unique gene symbols left after restriction:", n_int_restricted_unique_genes, "\n")
cat("Rows in interaction file before restriction:", nrow(int_dt), "\n")
cat("Rows in interaction file after  restriction:", nrow(int_dt_restricted), "\n")


## 1) GeneHancer-Dateien einlesen
gh_coords <- fread("GeneHancer_AnnotSV_hg19_v5.25.txt")
gh_gene   <- fread("GeneHancer_AnnotSV_gene_association_scores_v5.25.txt")


# ==========================================================
# Compare DATABASE vs CORRECTED BROWSER interactions
# corrected = interaction table restricted to GHids in RegElements
# ==========================================================

# --- SETS: GHids ---
ghids_db <- sort(unique(annot_elite$GHid[!is.na(annot_elite$GHid)]))
ghids_corr <- sort(unique(int_dt_restricted$GHid))

ghid_common <- intersect(ghids_db, ghids_corr)
ghid_only_db <- setdiff(ghids_db, ghids_corr)
ghid_only_corr <- setdiff(ghids_corr, ghids_db)

cat("\n==================== GHID COMPARISON (DB vs corrected browser) ====================\n")
cat("Database elite GHids:", length(ghids_db), "\n")
cat("Corrected browser GHids:", length(ghids_corr), "\n")
cat("Common GHids:", length(ghid_common), "\n")
cat("GHids only in DATABASE:", length(ghid_only_db), "\n")
cat("GHids only in corrected browser:", length(ghid_only_corr), "\n")

if (length(ghid_only_db) > 0) {
  cat("\nGHids only in DATABASE:\n")
  print(ghid_only_db)
}
if (length(ghid_only_corr) > 0) {
  cat("\nGHids only in corrected browser:\n")
  print(ghid_only_corr)
}

# --- SETS: genes ---
genes_db <- sort(unique(annot_elite$symbol[!is.na(annot_elite$symbol)]))
genes_corr <- sort(unique(int_dt_restricted$gene_symbol))

gene_common <- intersect(genes_db, genes_corr)
gene_only_db <- setdiff(genes_db, genes_corr)
gene_only_corr <- setdiff(genes_corr, genes_db)

cat("\n==================== GENE COMPARISON (DB vs corrected browser) ====================\n")
cat("Database elite genes:", length(genes_db), "\n")
cat("Corrected browser genes:", length(genes_corr), "\n")
cat("Common genes:", length(gene_common), "\n")
cat("Genes only in DATABASE:", length(gene_only_db), "\n")
cat("Genes only in corrected browser:", length(gene_only_corr), "\n")

if (length(gene_only_db) > 0) {
  cat("\nGenes only in DATABASE (present in gh_gene elite, missing from corrected browser):\n")
  print(gene_only_db)
}
if (length(gene_only_corr) > 0) {
  cat("\nGenes only in corrected browser (present in table browser output, missing from gh_gene elite strict):\n")
  print(gene_only_corr)
}

# --- PER-GHid gene differences ---
cat("\n==================== PER-GHID GENE SET DIFFERENCES ====================\n")

per_gh_diff <- rbindlist(lapply(sort(ghid_common), function(g) {
  g_db <- sort(unique(annot_elite[GHid == g & !is.na(symbol), symbol]))
  g_corr <- sort(unique(int_dt_restricted[GHid == g & !is.na(gene_symbol), gene_symbol]))
  
  data.table(
    GHid = g,
    n_genes_db = length(g_db),
    n_genes_corr = length(g_corr),
    genes_only_db = paste(setdiff(g_db, g_corr), collapse = ";"),
    genes_only_corr = paste(setdiff(g_corr, g_db), collapse = ";")
  )
}))

# show only GHids where something differs
per_gh_diff[
  (nchar(genes_only_db) > 0) | (nchar(genes_only_corr) > 0)
][order(-n_genes_db)][1:.N]

# ==========================================================
# BUILD + EXPORT SUMMARY TABLE (DB strict elite RE vs corrected UCSC)
# Robust: avoids cartesian merge explosions by enforcing unique keys
# Output: CpG_GeneHancer_DB_vs_corrected_UCSC_summary.xlsx
# ==========================================================

library(data.table)
library(GenomicRanges)
library(openxlsx)

# -------------------- helpers --------------------
collapse_unique <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("")
  paste(sort(unique(x)), collapse = ";")
}

split_set <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unique(strsplit(x, ";", fixed = TRUE)[[1]])
}

# ==================================================
# 1) DATABASE MAP: one row per (CpG x GHid)
#    Uses annot_elite produced from database pipeline
# ==================================================
db_map <- as.data.table(annot_elite)

# keep only real GHids; collapse elite symbols per CpG×GHid
db_map <- db_map[!is.na(GHid),
                 .(genes_DATABASE = collapse_unique(symbol)),
                 by = .(cpg_row, chrom, start, end, GHid)]

# ensure uniqueness
setkey(db_map, cpg_row, chrom, start, end, GHid)
db_map <- unique(db_map)

# ==================================================
# 2) CORRECTED UCSC MAP: one row per (CpG x GHid)
#    cpg_dt must exist with columns: cpg_row, chrom, start, end
#    int_dt_restricted must exist with columns: chrom,start,end, GHid, gene_symbol
# ==================================================
cpg_gr <- GRanges(seqnames = cpg_dt$chrom,
                  ranges   = IRanges(start = cpg_dt$start, end = cpg_dt$end))

int_gr <- GRanges(seqnames = int_dt_restricted$chrom,
                  ranges   = IRanges(start = int_dt_restricted$start, end = int_dt_restricted$end),
                  GHid     = int_dt_restricted$GHid,
                  gene     = int_dt_restricted$gene_symbol)

ov <- findOverlaps(cpg_gr, int_gr, ignore.strand = TRUE)

corr_hits <- data.table(
  cpg_row = queryHits(ov),
  GHid    = mcols(int_gr)$GHid[subjectHits(ov)],
  gene    = mcols(int_gr)$gene[subjectHits(ov)]
)

# attach CpG coords + collapse genes per CpG×GHid
corr_map <- merge(
  corr_hits,
  cpg_dt[, .(cpg_row, chrom, start, end)],
  by = "cpg_row",
  all.x = TRUE
)[,
  .(genes_corrected_UCSC = collapse_unique(gene)),
  by = .(cpg_row, chrom, start, end, GHid)
]

# ensure uniqueness
setkey(corr_map, cpg_row, chrom, start, end, GHid)
corr_map <- unique(corr_map)

# ==================================================
# 3) FULL JOIN (safe one-row-per-key on both sides)
# ==================================================
cmp <- merge(
  db_map,
  corr_map,
  by = c("cpg_row", "chrom", "start", "end", "GHid"),
  all = TRUE
)

# ==================================================
# 4) Compute gene differences (DB-only vs UCSC-only)
# ==================================================
cmp[, `:=`(
  genes_DATABASE = fifelse(is.na(genes_DATABASE), "", genes_DATABASE),
  genes_corrected_UCSC = fifelse(is.na(genes_corrected_UCSC), "", genes_corrected_UCSC)
)]

cmp[, `:=`(
  genes_only_DATABASE = {
    db <- split_set(genes_DATABASE); cr <- split_set(genes_corrected_UCSC)
    collapse_unique(setdiff(db, cr))
  },
  genes_only_corrected_UCSC = {
    db <- split_set(genes_DATABASE); cr <- split_set(genes_corrected_UCSC)
    collapse_unique(setdiff(cr, db))
  }
), by = seq_len(nrow(cmp))]

# ==================================================
# 5) Final summary table (sorted by CpG then GHid)
# ==================================================
summary_tbl <- cmp[
  ,
  .(
    chrom,
    start,
    end,
    GHid,
    genes_DATABASE,
    genes_corrected_UCSC,
    genes_only_DATABASE,
    genes_only_corrected_UCSC
  )
][order(chrom, start, GHid)]

print(summary_tbl)

# ==================================================
# 6) Export to Excel (nice formatting)
# ==================================================
out_file <- "CpG_GeneHancer_DB_vs_corrected_UCSC_summary.xlsx"

wb <- createWorkbook()
addWorksheet(wb, "CpG_GeneHancer_Comparison")

writeDataTable(
  wb,
  sheet = 1,
  x = summary_tbl,
  tableStyle = "TableStyleMedium9",
  withFilter = TRUE
)

freezePane(wb, sheet = 1, firstRow = TRUE)

setColWidths(
  wb,
  sheet = 1,
  cols = 1:ncol(summary_tbl),
  widths = "auto"
)

gene_cols <- which(names(summary_tbl) %in% c(
  "genes_DATABASE",
  "genes_corrected_UCSC",
  "genes_only_DATABASE",
  "genes_only_corrected_UCSC"
))

addStyle(
  wb,
  sheet = 1,
  style = createStyle(wrapText = TRUE, valign = "top"),
  rows = 2:(nrow(summary_tbl) + 1),
  cols = gene_cols,
  gridExpand = TRUE,
  stack = TRUE
)

saveWorkbook(wb, out_file, overwrite = TRUE)
cat("Excel file written to:", normalizePath(out_file), "\n")

