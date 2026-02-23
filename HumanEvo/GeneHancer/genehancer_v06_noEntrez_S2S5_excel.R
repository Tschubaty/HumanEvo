# =============================================================================
# GeneHancer annotation (DATABASE-ONLY, AnnotSV v5.25, hg19)
# =============================================================================
#
# GOAL
# ----
# Build a biologically consistent GeneHancer annotation pipeline using ONLY the
# GeneHancer database text files (AnnotSV v5.25).
#
# We define "double-elite" STRICTLY as:
#   (1) enhancer/regulatory element is elite
#       - if an "is_elite" column exists in gh_coords_file, require is_elite==1
#   (2) enhancer-gene association is elite (is_elite==1 in gh_gene_file),
#   (3) CpG must fall INSIDE the GeneHancer regulatory element coordinates (gh_coords).
#
# BACKGROUND (reviewer requirement)
# -------------------------------
# Background universe = all elite enhancer-gene associations where:
#   - the enhancer (GHid) overlaps CpGs inside the *H3K27ac peak CpG universe*
#     (your full CpG-in-peak set)
#   - AND the enhancer-gene association is elite (is_elite==1)
#   - AND (if available) the enhancer/regulatory element is elite (is_elite==1 in coords)
#
# FOREGROUND
# ----------
# BED3 files with CpGs
#   - Table S2 calculated from Analysis1_best_CpGs.bed
#   - Table S5 calculated from df_significant_Q_annot.bed
#
# OUTPUTS
# -------
# (1) Background: GHids + elite genes + enhancer coords + scores
# (2) Foregrounds: for EACH foreground, an Excel file with two tabs:
#     - "wide":  per-CpG collapsed annotation (CpG × all GHids × all elite genes)
#     - "long":  long annotation (one row per CpG × GHid × elite gene)
# (3) Convenience gene lists for enrichment tools (RAW symbols):
#     - background, Table S2, Table S5
#
# Author: Daniel Batyrev (HUJI 777634015)
# =============================================================================

# --------------------------- Set up environment ------------------------------
rm(list = ls())

cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
  picuture_file_extension <- "pdf"
} else {
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "png"
}
setwd(this.dir)

detachAllPackages <- function() {
  basic.packages <- c(
    "package:stats","package:graphics","package:grDevices","package:utils",
    "package:datasets","package:methods","package:base"
  )
  package.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]
  package.list <- setdiff(package.list, basic.packages)
  if (length(package.list) > 0) {
    for (package in package.list) detach(package, character.only = TRUE)
  }
}
detachAllPackages()

# ------------------------------ Libraries ------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(openxlsx)   # for Excel outputs (2 tabs per foreground)
})

# =============================================================================
# 0) INPUT PATHS
# =============================================================================

# --- GeneHancer database files (AnnotSV v5.25, hg19) ---
gh_coords_file <- "GeneHancer_AnnotSV_hg19_v5.25.txt"
gh_gene_file   <- "GeneHancer_AnnotSV_gene_association_scores_v5.25.txt"

# --- CpG-in-H3K27ac-peak universe (background CpG universe) ---
# Must contain columns: chrom, start, end (1bp CpGs are OK)
H3K27ac_cpg_universe_rds <- "C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/H3K27ac.only.45.samples.histogram.merged2025.hg19.rds"

# --- Foreground CpG lists (BED3) ---
# Table S2: from Analysis1_best_CpGs.bed
fg_S2_bed <- file.path(this.dir, "Analysis1_best_CpGs.bed")

# Table S5: from df_significant_Q_annot.bed
fg_S5_bed <- file.path(this.dir, "df_significant_Q_annot.bed")

# --- Output folder ---
out_dir <- file.path(this.dir, "GeneHancer_DB_only_outputs_v06_noEntrez")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1) LOAD DATABASE
# =============================================================================

cat("\n==================== 1) LOAD GeneHancer DATABASE ====================\n")

gh_coords <- fread(gh_coords_file)
gh_gene   <- fread(gh_gene_file)

# coordinate chromosome col
if ("chromosome" %in% names(gh_coords)) setnames(gh_coords, "chromosome", "chrom")
if (!("chrom" %in% names(gh_coords))) stop("gh_coords_file must contain a chromosome column named 'chromosome' or 'chrom'.")

# remove bad coord rows
gh_coords <- gh_coords[!is.na(start) & !is.na(end)]

# optional: require elite enhancer/regulatory element if present in coords
coords_has_elite <- "is_elite" %in% names(gh_coords)
if (coords_has_elite) {
  gh_coords <- gh_coords[is_elite == 1]
  cat("gh_coords: filtered to elite regulatory elements (is_elite==1 in coords)\n")
}

setkey(gh_coords, GHid)

# keep only elite enhancer-gene associations (STRICT)
if (!("is_elite" %in% names(gh_gene))) stop("gh_gene_file must contain an 'is_elite' column (elite enhancer-gene association).")
gh_gene_elite <- gh_gene[is_elite == 1]
setkey(gh_gene_elite, GHid)

cat("gh_coords rows:", nrow(gh_coords), "\n")
cat("gh_gene rows:", nrow(gh_gene), "\n")
cat("gh_gene_elite rows (is_elite==1):", nrow(gh_gene_elite), "\n")

# attach coordinates to elite associations
gh_elite_annot <- merge(
  gh_gene_elite,
  gh_coords,
  by = "GHid",
  all.x = TRUE,
  allow.cartesian = TRUE
)

# some GHids may have missing coordinates
gh_elite_annot <- gh_elite_annot[!is.na(start) & !is.na(end)]
setkey(gh_elite_annot, GHid)

cat("gh_elite_annot rows after coord-attach:", nrow(gh_elite_annot), "\n")
cat("Unique elite GHids in database:", uniqueN(gh_elite_annot$GHid), "\n")
cat("Unique elite gene symbols in database:", uniqueN(gh_elite_annot$symbol), "\n")

# =============================================================================
# 2) BUILD BACKGROUND UNIVERSE
#    (elite enhancers that overlap CpGs in H3K27ac peaks)
# =============================================================================

cat("\n==================== 2) BUILD BACKGROUND UNIVERSE ====================\n")

H3K27ac_cpg_universe <- readRDS(H3K27ac_cpg_universe_rds)
stopifnot(all(c("chrom","start","end") %in% names(H3K27ac_cpg_universe)))

# CpG universe GRanges
cpg_bg_gr <- with(H3K27ac_cpg_universe,
                  GRanges(seqnames = chrom,
                          ranges   = IRanges(start = start, end = end)))

# elite GH regulatory elements GRanges (per association row; OK)
gh_bg_gr <- with(gh_elite_annot,
                 GRanges(seqnames = chrom,
                         ranges   = IRanges(start = start, end = end),
                         GHid = GHid,
                         symbol = symbol,
                         combined_score = combined_score))

# Overlap: CpG ∩ GH enhancer coords
ov_bg <- findOverlaps(cpg_bg_gr, gh_bg_gr, ignore.strand = TRUE)

bg_hits <- data.table(
  GHid           = mcols(gh_bg_gr)$GHid[subjectHits(ov_bg)],
  symbol         = mcols(gh_bg_gr)$symbol[subjectHits(ov_bg)],
  combined_score = mcols(gh_bg_gr)$combined_score[subjectHits(ov_bg)],
  chrom          = as.character(seqnames(gh_bg_gr))[subjectHits(ov_bg)],
  start          = start(gh_bg_gr)[subjectHits(ov_bg)],
  end            = end(gh_bg_gr)[subjectHits(ov_bg)]
)

background_db <- unique(bg_hits)
setorder(background_db, symbol, GHid, chrom, start)

background_genes <- sort(unique(background_db$symbol))
background_ghids <- sort(unique(background_db$GHid))

cat("Unique background elite GHids (overlap CpG universe):", length(background_ghids), "\n")
cat("Unique background elite genes:", length(background_genes), "\n")

# Save background universe
saveRDS(background_db, file.path(out_dir, "background_double_elite_GHid_gene_pairs_overlap_H3K27acCpGs.hg19.rds"))
fwrite(background_db, file.path(out_dir, "background_double_elite_GHid_gene_pairs_overlap_H3K27acCpGs.hg19.tsv"), sep = "\t")
write.table(background_genes,
            file.path(out_dir, "background_double_elite_genes_overlap_H3K27acCpGs.raw_symbols.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# =============================================================================
# 3) FOREGROUND ANNOTATION FUNCTION (DATABASE-ONLY, DOUBLE-ELITE)
# =============================================================================

cat("\n==================== 3) DEFINE FOREGROUND ANNOTATION ====================\n")

read_bed3 <- function(path) {
  dt <- data.table::fread(path, header = FALSE, select = 1:3)
  data.table::setnames(dt, c("chrom","start","end"))
  dt[, `:=`(start = as.integer(start), end = as.integer(end))]
  dt
}


annotate_cpgs_db_elite <- function(cpg_dt, gh_elite_annot_dt, label = "FG") {
  cpg_dt <- as.data.table(cpg_dt)
  stopifnot(all(c("chrom","start","end") %in% names(cpg_dt)))
  
  # stable CpG row id
  cpg_dt <- unique(cpg_dt[, .(chrom, start, end)])
  cpg_dt[, cpg_row := .I]
  setkey(cpg_dt, cpg_row)
  
  # GRanges for CpGs
  cpg_gr <- GRanges(
    seqnames = cpg_dt$chrom,
    ranges   = IRanges(start = cpg_dt$start, end = cpg_dt$end)
  )
  
  # GH elite rows (coords already from BED; may be multiple intervals per GHid)
  gh_dt <- as.data.table(gh_elite_annot_dt)
  gh_dt <- gh_dt[!is.na(chrom) & !is.na(start) & !is.na(end)]
  
  # Identify association score column (same logic as in load step, but safe here too)
  assoc_score_col <- intersect(
    c("combined_score"),
    names(gh_dt)
  )
  if (length(assoc_score_col) < 1) {
    stop(paste0(
      "Could not find association score column in gh_elite_annot_dt. Available columns:\n",
      paste(names(gh_dt), collapse = ", ")
    ))
  }
  assoc_score_col <- assoc_score_col[1]
  
  gh_gr <- GRanges(
    seqnames = gh_dt$chrom,
    ranges   = IRanges(start = gh_dt$start, end = gh_dt$end),
    GHid           = gh_dt$GHid,
    symbol         = gh_dt$symbol,
    combined_score = gh_dt$combined_score,
    assoc_score    = gh_dt[[assoc_score_col]]
  )
  
  ov <- findOverlaps(cpg_gr, gh_gr, ignore.strand = TRUE)
  
  # LONG: use GH coords (from BED) as the main chrom/start/end
  annot_long <- data.table(
    cpg_row = queryHits(ov),
    
    GHid    = mcols(gh_gr)$GHid[subjectHits(ov)],
    symbol  = mcols(gh_gr)$symbol[subjectHits(ov)],
    
    # scores
    combined_score = mcols(gh_gr)$combined_score[subjectHits(ov)],
    association_score = mcols(gh_gr)$assoc_score[subjectHits(ov)],
    
    # GH enhancer coords (requested)
    chrom = as.character(seqnames(gh_gr))[subjectHits(ov)],
    start = start(gh_gr)[subjectHits(ov)],
    end   = end(gh_gr)[subjectHits(ov)]
  )
  
  # attach CpG coords too (as cpg_* so long has both, but GH coords are primary)
  annot_long <- merge(
    annot_long,
    cpg_dt[, .(cpg_row, cpg_chrom = chrom, cpg_start = start, cpg_end = end)],
    by = "cpg_row",
    all.x = TRUE
  )
  
  annot_long <- unique(annot_long)
  setorder(annot_long, chrom, start, end, GHid, symbol, cpg_chrom, cpg_start, cpg_end)
  
  # WIDE: keep CpG coords as the unit; DO NOT include association_score here
  # (combined_score summary remains OK; association_score is per gene↔GHid)
  annot_cpg_collapsed <- annot_long[
    ,
    .(
      elite_GHids_all  = paste(sort(unique(GHid)), collapse = ";"),
      elite_genes_all  = paste(sort(unique(symbol)), collapse = ";"),
      elite_n_GHid     = uniqueN(GHid),
      elite_n_genes    = uniqueN(symbol),
      elite_scores_all = paste(
        sapply(sort(unique(GHid)), function(g) {
          sc <- suppressWarnings(max(combined_score[GHid == g], na.rm = TRUE))
          if (is.infinite(sc)) sc <- NA_real_
          paste0(g, ":", sc)
        }),
        collapse = ";"
      )
    ),
    by = .(cpg_row, cpg_chrom, cpg_start, cpg_end)
  ]
  
  # ensure CpGs with 0 hits are preserved (wide is CpG-centric)
  annot_cpg_all <- merge(
    cpg_dt[, .(cpg_row, cpg_chrom = chrom, cpg_start = start, cpg_end = end)],
    annot_cpg_collapsed,
    by = c("cpg_row","cpg_chrom","cpg_start","cpg_end"),
    all.x = TRUE
  )
  annot_cpg_all[is.na(elite_GHids_all), elite_GHids_all := ""]
  annot_cpg_all[is.na(elite_genes_all), elite_genes_all := ""]
  annot_cpg_all[is.na(elite_scores_all), elite_scores_all := ""]
  annot_cpg_all[is.na(elite_n_GHid), elite_n_GHid := 0L]
  annot_cpg_all[is.na(elite_n_genes), elite_n_genes := 0L]
  
  cat("\n----", label, "DATABASE-ONLY double-elite annotation ----\n")
  cat("Unique CpGs input:", nrow(cpg_dt), "\n")
  cat("CpGs with >=1 elite hit:", uniqueN(annot_long$cpg_row), "of", nrow(cpg_dt), "\n")
  cat("Unique elite GHids:", uniqueN(annot_long$GHid), "\n")
  cat("Unique elite genes:", uniqueN(annot_long$symbol), "\n")
  
  list(
    cpg_dt = cpg_dt,
    annot_long = annot_long,
    annot_cpg_all = annot_cpg_all
  )
}


# =============================================================================
# 4) LOAD FOREGROUNDS (S2, S5) + ANNOTATE
# =============================================================================

cat("\n==================== 4) LOAD + ANNOTATE FOREGROUNDS ====================\n")

if (!file.exists(fg_S2_bed)) stop(paste0("Missing foreground BED for Table S2: ", fg_S2_bed))
if (!file.exists(fg_S5_bed)) stop(paste0("Missing foreground BED for Table S5: ", fg_S5_bed))

fg_S2_cpg <- read_bed3(fg_S2_bed)
fg_S5_cpg <- read_bed3(fg_S5_bed)

res_S2 <- annotate_cpgs_db_elite(fg_S2_cpg, gh_elite_annot, label = "Table S2 (Analysis1_best_CpGs.bed)")
res_S5 <- annotate_cpgs_db_elite(fg_S5_cpg, gh_elite_annot, label = "Table S5 (df_significant_Q_annot.bed)")

# =============================================================================
# 5) SAVE FOREGROUND OUTPUTS (EXCEL with wide + long) + TSV + RAW GENE LISTS
# =============================================================================

cat("\n==================== 5) WRITE OUTPUTS ====================\n")

write_fg_outputs <- function(res, fg_label, out_prefix) {
  # files
  xlsx_path <- file.path(out_dir, paste0(out_prefix, "_GeneHancer_double_elite.xlsx"))
  long_tsv  <- file.path(out_dir, paste0(out_prefix, "_GeneHancer_double_elite_long.tsv"))
  wide_tsv  <- file.path(out_dir, paste0(out_prefix, "_GeneHancer_double_elite_wide.tsv"))
  genes_txt <- file.path(out_dir, paste0(out_prefix, "_elite_genes_raw_symbols.txt"))

  # write TSVs
  fwrite(res$annot_long,    long_tsv, sep = "\t")
  fwrite(res$annot_cpg_all, wide_tsv, sep = "\t")

  # gene list (raw symbols)
  fg_genes <- sort(unique(res$annot_long$symbol))
  fg_genes <- fg_genes[!is.na(fg_genes) & fg_genes != ""]
  writeLines(fg_genes, con = genes_txt, useBytes = TRUE)

  # Excel with 2 tabs
  wb <- createWorkbook()
  addWorksheet(wb, "wide")
  addWorksheet(wb, "long")
  writeDataTable(wb, "wide", res$annot_cpg_all)
  writeDataTable(wb, "long", res$annot_long)
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  cat("\n---", fg_label, "---\n")
  cat("Excel:", normalizePath(xlsx_path, winslash = "/"), "\n")
  cat("Wide TSV:", normalizePath(wide_tsv, winslash = "/"), "\n")
  cat("Long TSV:", normalizePath(long_tsv, winslash = "/"), "\n")
  cat("Raw gene list:", normalizePath(genes_txt, winslash = "/"), "\n")
  cat("Unique genes:", length(fg_genes), "\n")
}

# Foreground outputs
write_fg_outputs(res_S2, "Table S2", "S2")
write_fg_outputs(res_S5, "Table S5", "S5")

# Also save RDS objects for reproducibility
saveRDS(res_S2, file.path(out_dir, "S2_db_only_double_elite_annotation.rds"))
saveRDS(res_S5, file.path(out_dir, "S5_db_only_double_elite_annotation.rds"))

# Convenience gene lists for enrichment tools: background + foregrounds (raw symbols)
writeLines(sort(unique(background_db$symbol)),
           con = file.path(out_dir, "BACKGROUND_elite_genes_raw_symbols.txt"),
           useBytes = TRUE)

cat("\nDONE. Outputs written to:\n", normalizePath(out_dir, winslash = "/"), "\n")
