# =============================================================================
# GeneHancer annotation (DATABASE-ONLY, AnnotSV v5.25, hg19)
# =============================================================================
#
# GOAL
# ----
# Build a biologically consistent GeneHancer annotation pipeline using ONLY the
# GeneHancer database text files (AnnotSV v5.25), without UCSC interaction tables.
#
# We define "double-elite" STRICTLY as:
#   (1) enhancer/regulatory element is elite (is_elite==1 for the enhancer-gene pair),
#   (2) enhancer-gene association is elite (is_elite==1 in gh_gene association table),
#   (3) CpG must fall INSIDE the GeneHancer regulatory element coordinates (gh_coords).
#
# BACKGROUND (reviewer requirement)
# -------------------------------
# Background universe = all elite enhancer-gene associations where:
#   - the enhancer (GHid) overlaps CpGs inside the *H3K27ac peak CpG universe* (your
#     full CpG-in-peak set)
#   - AND the enhancer-gene association is elite (is_elite==1)
#
# FOREGROUND
# ----------
# Two CpG foreground lists are annotated using the SAME rules.
# Provide two inputs (choose ONE of the two options below):
#   Option A: BED3 files with CpGs (recommended)
#   Option B: Excel tables with CpGs (chrom/start/end columns)
#
# OUTPUTS
# -------
# (1) Background: GHids + elite genes + coords + optional scores
# (2) Foreground1, Foreground2: CpG × GHid × elite genes (long + collapsed)
# (3) Convenience gene lists for enrichment tools
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
})

# =============================================================================
# 0) INPUT PATHS
# =============================================================================

# --- GeneHancer database files (AnnotSV v5.25, hg19) ---
gh_coords_file <- "GeneHancer_AnnotSV_hg19_v5.25.txt"
gh_gene_file   <- "GeneHancer_AnnotSV_gene_association_scores_v5.25.txt"

# --- Your CpG-in-H3K27ac-peak universe (background CpG universe) ---
# Must contain columns: chrom, start, end (1bp CpGs are OK)
H3K27ac_cpg_universe_rds <- "C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/methylation+chip/H3K27ac.only.45.samples.histogram.merged2025.hg19.rds"

# --- Foreground CpG lists: choose A or B ---
# Option A: BED3 files (recommended)
fg1_bed <- file.path(this.dir, "df_significant_Q_annot.bed")  # example (edit if needed)
fg2_bed <- file.path(this.dir, "df_significant_Q_annot.bed")  # placeholder; set to your 2nd CpG list bed

# Option B: Excel tables (if you want to reproduce S2/S5)
file_S2 <- "C:/Users/Batyrev/Dropbox/R42v/Table S2 significant CpG [1b].xlsx"
file_S5 <- "C:/Users/Batyrev/Dropbox/R42v/Table S5 significant CpG [2b].xlsx"

# --- Output folder ---
out_dir <- file.path(this.dir, "GeneHancer_DB_only_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1) LOAD DATABASE
# =============================================================================

cat("\n==================== 1) LOAD GeneHancer DATABASE ====================\n")

gh_coords <- fread(gh_coords_file)
gh_gene   <- fread(gh_gene_file)

# sanity: rename coordinate chromosome column
if ("chromosome" %in% names(gh_coords)) setnames(gh_coords, "chromosome", "chrom")

# remove bad coord rows (NA start/end)
gh_coords <- gh_coords[!is.na(start) & !is.na(end)]
setkey(gh_coords, GHid)

# keep only elite enhancer-gene associations (STRICT)
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

# some GHids may have missing coordinates even after gh_coords cleanup (rare)
gh_elite_annot <- gh_elite_annot[!is.na(start) & !is.na(end)]
setkey(gh_elite_annot, GHid)

cat("gh_elite_annot rows after coord-attach:", nrow(gh_elite_annot), "\n")
cat("Unique elite GHids in database:", uniqueN(gh_elite_annot$GHid), "\n")
cat("Unique elite gene symbols in database:", uniqueN(gh_elite_annot$symbol), "\n")

# =============================================================================
# 2) BUILD BACKGROUND UNIVERSE (elite enhancers that overlap CpGs in H3K27ac peaks)
# =============================================================================

cat("\n==================== 2) BUILD BACKGROUND UNIVERSE ====================\n")

# Load CpG universe (all CpGs in H3K27ac peaks)
H3K27ac_cpg_universe <- readRDS(H3K27ac_cpg_universe_rds)

# Convert to GRanges (assumes columns chrom/start/end exist)
cpg_bg_gr <- with(H3K27ac_cpg_universe,
                  GRanges(seqnames = chrom,
                          ranges   = IRanges(start = start, end = end)))

# Build GRanges for elite regulatory elements (one record per association row here,
# but with same GHid coords repeated; that's OK for overlaps, we will deduplicate later)
gh_bg_gr <- with(gh_elite_annot,
                 GRanges(seqnames = chrom,
                         ranges   = IRanges(start = start, end = end),
                         GHid = GHid,
                         symbol = symbol,
                         combined_score = combined_score))

# Overlap: which elite GH regulatory elements contain at least one CpG in H3K27ac peak CpG universe?
ov_bg <- findOverlaps(cpg_bg_gr, gh_bg_gr, ignore.strand = TRUE)

bg_hits <- data.table(
  GHid           = mcols(gh_bg_gr)$GHid[subjectHits(ov_bg)],
  symbol         = mcols(gh_bg_gr)$symbol[subjectHits(ov_bg)],
  combined_score = mcols(gh_bg_gr)$combined_score[subjectHits(ov_bg)],
  chrom          = as.character(seqnames(gh_bg_gr))[subjectHits(ov_bg)],
  start          = start(gh_bg_gr)[subjectHits(ov_bg)],
  end            = end(gh_bg_gr)[subjectHits(ov_bg)]
)

# Background universe = unique elite (GHid, symbol) pairs whose enhancer overlaps CpG universe
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
            file.path(out_dir, "background_double_elite_genes_overlap_H3K27acCpGs.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# =============================================================================
# 3) FOREGROUND ANNOTATION FUNCTION (DATABASE-ONLY)
# =============================================================================

cat("\n==================== 3) DEFINE FOREGROUND ANNOTATION ====================\n")

# Annotate CpGs to elite GeneHancer associations (STRICT):
#   CpG overlaps GH regulatory element coords
#   AND the enhancer-gene association is elite (already in gh_elite_annot)
annotate_cpgs_db_elite <- function(cpg_dt, gh_elite_annot_dt, label = "FG") {
  cpg_dt <- as.data.table(cpg_dt)
  stopifnot(all(c("chrom","start","end") %in% names(cpg_dt)))
  
  # stable CpG row id
  cpg_dt <- unique(cpg_dt[, .(chrom, start, end)])
  cpg_dt[, cpg_row := .I]
  setkey(cpg_dt, cpg_row)
  
  # GRanges for CpGs
  cpg_gr <- GRanges(seqnames = cpg_dt$chrom,
                    ranges   = IRanges(start = cpg_dt$start, end = cpg_dt$end))
  
  # GRanges for elite GH associations
  gh_dt <- as.data.table(gh_elite_annot_dt)
  gh_dt <- gh_dt[!is.na(start) & !is.na(end)]
  
  gh_gr <- GRanges(seqnames = gh_dt$chrom,
                   ranges   = IRanges(start = gh_dt$start, end = gh_dt$end),
                   GHid = gh_dt$GHid,
                   symbol = gh_dt$symbol,
                   combined_score = gh_dt$combined_score)
  
  # overlaps
  ov <- findOverlaps(cpg_gr, gh_gr, ignore.strand = TRUE)
  
  # long table: one row per CpG x (GHid, elite gene association)
  annot_long <- data.table(
    cpg_row = queryHits(ov),
    GHid    = mcols(gh_gr)$GHid[subjectHits(ov)],
    symbol  = mcols(gh_gr)$symbol[subjectHits(ov)],
    combined_score = mcols(gh_gr)$combined_score[subjectHits(ov)]
  )
  
  # attach CpG coords
  annot_long <- merge(
    annot_long,
    cpg_dt,
    by = "cpg_row",
    all.x = TRUE
  )
  
  # keep unique rows
  annot_long <- unique(annot_long)
  setorder(annot_long, chrom, start, end, GHid, symbol)
  
  # collapsed per CpG (all GHids / genes)
  annot_cpg_collapsed <- annot_long[
    ,
    .(
      elite_GHids_all  = paste(sort(unique(GHid)), collapse = ";"),
      elite_genes_all  = paste(sort(unique(symbol)), collapse = ";"),
      elite_n_GHid     = uniqueN(GHid),
      elite_n_genes    = uniqueN(symbol)
    ),
    by = .(cpg_row, chrom, start, end)
  ]
  
  # ensure CpGs with 0 hits are preserved
  annot_cpg_all <- merge(
    cpg_dt[, .(cpg_row, chrom, start, end)],
    annot_cpg_collapsed,
    by = c("cpg_row","chrom","start","end"),
    all.x = TRUE
  )
  annot_cpg_all[is.na(elite_GHids_all), elite_GHids_all := ""]
  annot_cpg_all[is.na(elite_genes_all), elite_genes_all := ""]
  annot_cpg_all[is.na(elite_n_GHid), elite_n_GHid := 0L]
  annot_cpg_all[is.na(elite_n_genes), elite_n_genes := 0L]
  
  # summary cats
  cat("\n----", label, "DATABASE-ONLY elite annotation ----\n")
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
# 4) LOAD TWO FOREGROUND CpG LISTS + ANNOTATE (DATABASE-ONLY)
# =============================================================================

cat("\n==================== 4) LOAD + ANNOTATE FOREGROUNDS ====================\n")

# ---------------- Option A: BED3 foregrounds ----------------
read_bed3 <- function(path) {
  # BED can have >3 columns; keep first 3
  dt <- fread(path, header = FALSE)
  setnames(dt, 1:3, c("chrom","start","end"))
  dt[, .(chrom, as.integer(start), as.integer(end))]
}

use_bed_inputs <- file.exists(fg1_bed) && file.exists(fg2_bed)

if (use_bed_inputs) {
  fg1_cpg <- read_bed3(fg1_bed)
  fg2_cpg <- read_bed3(fg2_bed)
  
} else {
  # ---------------- Option B: Excel foregrounds (S2/S5) ----------------
  suppressPackageStartupMessages(library(readxl))
  
  fg1_tbl <- read_excel(file_S2, sheet = 1) |>
    dplyr::select(chrom, start, end)
  fg2_tbl <- read_excel(file_S5, sheet = 1) |>
    dplyr::select(chrom, start, end)
  
  fg1_cpg <- as.data.table(fg1_tbl)
  fg2_cpg <- as.data.table(fg2_tbl)
}

# Annotate
fg1_res <- annotate_cpgs_db_elite(fg1_cpg, gh_elite_annot, label = "FG1")
fg2_res <- annotate_cpgs_db_elite(fg2_cpg, gh_elite_annot, label = "FG2")

# Save outputs
saveRDS(fg1_res, file.path(out_dir, "FG1_db_only_elite_annotation.rds"))
saveRDS(fg2_res, file.path(out_dir, "FG2_db_only_elite_annotation.rds"))

fwrite(fg1_res$annot_long,    file.path(out_dir, "FG1_db_only_elite_annotation_long.tsv"), sep = "\t")
fwrite(fg2_res$annot_long,    file.path(out_dir, "FG2_db_only_elite_annotation_long.tsv"), sep = "\t")
fwrite(fg1_res$annot_cpg_all, file.path(out_dir, "FG1_db_only_elite_annotation_perCpG.tsv"), sep = "\t")
fwrite(fg2_res$annot_cpg_all, file.path(out_dir, "FG2_db_only_elite_annotation_perCpG.tsv"), sep = "\t")

# Gene lists
fg1_genes <- sort(unique(fg1_res$annot_long$symbol))
fg2_genes <- sort(unique(fg2_res$annot_long$symbol))

write.table(fg1_genes, file.path(out_dir, "FG1_elite_genes_db_only.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(fg2_genes, file.path(out_dir, "FG2_elite_genes_db_only.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# =============================================================================
# 5) MAP Gene symbols -> EntrezID (background + FG1 + FG2) WITH DEBUG PRINTS
# =============================================================================

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(data.table)
})

# --------------------------- helper: safe % formatting ---------------------------
pct <- function(x, denom) sprintf("%.2f%%", 100 * x / max(1, denom))

# --------------------------- helper: mapping function ---------------------------
map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(trimws(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  
  out <- data.table(symbol_input = symbols, ENTREZID = NA_character_, map_keytype = NA_character_)
  
  # (A) ENSEMBL-like IDs
  is_ensg <- grepl("^ENSG\\d+", symbols)
  if (any(is_ensg)) {
    m <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = unique(symbols[is_ensg]),
      keytype  = "ENSEMBL",
      columns  = "ENTREZID"
    )
    m <- as.data.table(m)[!is.na(ENTREZID), .(symbol_input = ENSEMBL, ENTREZID)]
    if (nrow(m) > 0) out[m, on = "symbol_input", `:=`(ENTREZID = i.ENTREZID, map_keytype = "ENSEMBL")]
  }
  
  # (B) SYMBOL
  rem <- is.na(out$ENTREZID)
  if (any(rem)) {
    m <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = unique(out$symbol_input[rem]),
      keytype  = "SYMBOL",
      columns  = "ENTREZID"
    )
    m <- as.data.table(m)[!is.na(ENTREZID), .(symbol_input = SYMBOL, ENTREZID)]
    if (nrow(m) > 0) out[m, on = "symbol_input", `:=`(ENTREZID = i.ENTREZID, map_keytype = "SYMBOL")]
  }
  
  # (C) ALIAS
  rem <- is.na(out$ENTREZID)
  if (any(rem)) {
    m <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = unique(out$symbol_input[rem]),
      keytype  = "ALIAS",
      columns  = "ENTREZID"
    )
    m <- as.data.table(m)[!is.na(ENTREZID), .(symbol_input = ALIAS, ENTREZID)]
    if (nrow(m) > 0) out[m, on = "symbol_input", `:=`(ENTREZID = i.ENTREZID, map_keytype = "ALIAS")]
  }
  
  out[, ENTREZID := as.character(ENTREZID)]
  out[]
}

# --------------------------- helper: debug one set ---------------------------
debug_set_mapping <- function(set_name, symbols_vec, sym2entrez_dt) {
  x <- unique(trimws(as.character(symbols_vec)))
  x <- x[!is.na(x) & x != ""]
  
  # join
  tmp <- data.table(symbol_input = x)
  tmp <- merge(tmp, sym2entrez_dt, by = "symbol_input", all.x = TRUE)
  
  n_in <- length(x)
  n_mapped <- sum(!is.na(tmp$ENTREZID))
  n_unmapped <- n_in - n_mapped
  
  # map type
  kt <- table(tmp$map_keytype, useNA = "ifany")
  
  # unique Entrez (numeric)
  e <- tmp$ENTREZID[!is.na(tmp$ENTREZID)]
  e_num <- suppressWarnings(as.numeric(e))
  e_num <- e_num[!is.na(e_num)]
  e_unique <- sort(unique(e_num))
  
  n_entrez_unique <- length(e_unique)
  n_lt1m <- sum(e_unique < 1e6)
  n_ge1m <- sum(e_unique >= 1e6)
  
  cat("\n=================================================================\n")
  cat("SET:", set_name, "\n")
  cat("-----------------------------------------------------------------\n")
  cat("Unique input symbols:", n_in, "\n")
  cat("Mapped symbols       :", n_mapped, "(", pct(n_mapped, n_in), ")\n")
  cat("Unmapped symbols     :", n_unmapped, "(", pct(n_unmapped, n_in), ")\n")
  cat("Unique numeric Entrez:", n_entrez_unique, "\n")
  cat("  - Entrez < 1e6     :", n_lt1m, "(", pct(n_lt1m, n_entrez_unique), "of mapped unique)\n")
  cat("  - Entrez >= 1e6    :", n_ge1m, "(", pct(n_ge1m, n_entrez_unique), "of mapped unique)\n")
  cat("Keytype usage:\n")
  print(kt)
  
  # return details for saving
  list(
    joined = tmp,
    unmapped = sort(unique(tmp[is.na(ENTREZID), symbol_input])),
    entrez_unique_numeric = e_unique
  )
}

# =============================================================================
# (1) Collect symbols across background + FG1 + FG2, build one master mapping
# =============================================================================
all_symbols <- unique(c(
  background_db$symbol,
  fg1_res$annot_long$symbol,
  fg2_res$annot_long$symbol
))
sym2entrez <- map_symbols_to_entrez(all_symbols)

cat("\n==================== MASTER MAPPING (ALL SYMBOLS TO MAP) ====================\n")
cat("Total unique symbols requested:", length(unique(all_symbols)), "\n")
cat("Mapped:", sum(!is.na(sym2entrez$ENTREZID)), "(", pct(sum(!is.na(sym2entrez$ENTREZID)), nrow(sym2entrez)), ")\n")
cat("Unmapped:", sum(is.na(sym2entrez$ENTREZID)), "(", pct(sum(is.na(sym2entrez$ENTREZID)), nrow(sym2entrez)), ")\n")
cat("Keytype usage overall:\n")
print(table(sym2entrez$map_keytype, useNA = "ifany"))

# Save the master mapping
fwrite(sym2entrez, file.path(out_dir, "symbol_to_entrez_mapping.tsv"), sep = "\t")
saveRDS(sym2entrez, file.path(out_dir, "symbol_to_entrez_mapping.rds"))

# =============================================================================
# (2) Debug per-set mapping quality + save unmapped lists per set
# =============================================================================
dbg_bg  <- debug_set_mapping("BACKGROUND (background_db$symbol)", background_db$symbol, sym2entrez)
dbg_fg1 <- debug_set_mapping("FOREGROUND 1 (fg1_res$annot_long$symbol)", fg1_res$annot_long$symbol, sym2entrez)
dbg_fg2 <- debug_set_mapping("FOREGROUND 2 (fg2_res$annot_long$symbol)", fg2_res$annot_long$symbol, sym2entrez)

write.table(dbg_bg$unmapped,  file.path(out_dir, "unmapped_symbols_BACKGROUND.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(dbg_fg1$unmapped, file.path(out_dir, "unmapped_symbols_FG1.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(dbg_fg2$unmapped, file.path(out_dir, "unmapped_symbols_FG2.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(dbg_bg$entrez_unique_numeric,  file.path(out_dir, "BACKGROUND_entrez_unique_numeric.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(dbg_fg1$entrez_unique_numeric, file.path(out_dir, "FG1_entrez_unique_numeric.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(dbg_fg2$entrez_unique_numeric, file.path(out_dir, "FG2_entrez_unique_numeric.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# =============================================================================
# (3) Add EntrezID columns to long tables (background + FG1 + FG2) + export
# =============================================================================

# ---- BACKGROUND ----
background_db_entrez <- copy(as.data.table(background_db))
background_db_entrez <- merge(
  background_db_entrez,
  sym2entrez,
  by.x = "symbol",
  by.y = "symbol_input",
  all.x = TRUE
)
setnames(background_db_entrez, c("ENTREZID","map_keytype"), c("EntrezID","Entrez_map_keytype"))

cat("\n==================== BACKGROUND TABLE WITH ENTREZ ====================\n")
cat("Rows:", nrow(background_db_entrez), "\n")
cat("Unique symbols:", uniqueN(background_db_entrez$symbol), "\n")
cat("Unique Entrez (non-NA):", uniqueN(background_db_entrez$EntrezID[!is.na(background_db_entrez$EntrezID) & background_db_entrez$EntrezID != ""]), "\n")
cat("Symbols unmapped in background:", uniqueN(background_db_entrez[is.na(EntrezID) | EntrezID == "", symbol]), "\n")

saveRDS(background_db_entrez, file.path(out_dir, "background_double_elite_with_entrez.rds"))
fwrite(background_db_entrez, file.path(out_dir, "background_double_elite_with_entrez.tsv"), sep = "\t")

# ---- FG1 LONG ----
fg1_long_entrez <- copy(as.data.table(fg1_res$annot_long))
fg1_long_entrez <- merge(
  fg1_long_entrez,
  sym2entrez,
  by.x = "symbol",
  by.y = "symbol_input",
  all.x = TRUE
)
setnames(fg1_long_entrez, c("ENTREZID","map_keytype"), c("EntrezID","Entrez_map_keytype"))

cat("\n==================== FG1 LONG TABLE WITH ENTREZ ====================\n")
cat("Rows:", nrow(fg1_long_entrez), "\n")
cat("Unique symbols:", uniqueN(fg1_long_entrez$symbol), "\n")
cat("Unique Entrez (non-NA):", uniqueN(fg1_long_entrez$EntrezID[!is.na(fg1_long_entrez$EntrezID) & fg1_long_entrez$EntrezID != ""]), "\n")
cat("Symbols unmapped in FG1:", uniqueN(fg1_long_entrez[is.na(EntrezID) | EntrezID == "", symbol]), "\n")

saveRDS(fg1_long_entrez, file.path(out_dir, "FG1_db_only_elite_annotation_long_with_entrez.rds"))
fwrite(fg1_long_entrez, file.path(out_dir, "FG1_db_only_elite_annotation_long_with_entrez.tsv"), sep = "\t")

# ---- FG2 LONG ----
fg2_long_entrez <- copy(as.data.table(fg2_res$annot_long))
fg2_long_entrez <- merge(
  fg2_long_entrez,
  sym2entrez,
  by.x = "symbol",
  by.y = "symbol_input",
  all.x = TRUE
)
setnames(fg2_long_entrez, c("ENTREZID","map_keytype"), c("EntrezID","Entrez_map_keytype"))

cat("\n==================== FG2 LONG TABLE WITH ENTREZ ====================\n")
cat("Rows:", nrow(fg2_long_entrez), "\n")
cat("Unique symbols:", uniqueN(fg2_long_entrez$symbol), "\n")
cat("Unique Entrez (non-NA):", uniqueN(fg2_long_entrez$EntrezID[!is.na(fg2_long_entrez$EntrezID) & fg2_long_entrez$EntrezID != ""]), "\n")
cat("Symbols unmapped in FG2:", uniqueN(fg2_long_entrez[is.na(EntrezID) | EntrezID == "", symbol]), "\n")

saveRDS(fg2_long_entrez, file.path(out_dir, "FG2_db_only_elite_annotation_long_with_entrez.rds"))
fwrite(fg2_long_entrez, file.path(out_dir, "FG2_db_only_elite_annotation_long_with_entrez.tsv"), sep = "\t")

cat("\nDONE: Mapping + per-set debug + Entrez columns added. Outputs in:\n", out_dir, "\n")

# =============================================================================
# 6) EnrichR submission files (EntrezID lists): BACKGROUND + FG1 + FG2
#     - one EntrezID per line
#     - writes 3 txt files into a dedicated subfolder under out_dir
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

enrichr_dir <- file.path(out_dir, "EnrichR_submission")
dir.create(enrichr_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helper: clean + write Entrez list ----
write_entrez_list <- function(entrez_vec, out_path, label) {
  entrez_vec <- unique(as.character(entrez_vec))
  entrez_vec <- trimws(entrez_vec)
  entrez_vec <- entrez_vec[!is.na(entrez_vec) & entrez_vec != ""]
  
  # keep numeric-only (EnrichR expects Entrez IDs)
  entrez_num <- suppressWarnings(as.numeric(entrez_vec))
  bad_non_numeric <- entrez_vec[is.na(entrez_num)]
  entrez_num <- entrez_num[!is.na(entrez_num)]
  entrez_num <- sort(unique(entrez_num))
  
  write.table(entrez_num, out_path, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("\n==================== EnrichR file written ====================\n")
  cat("Set:", label, "\n")
  cat("Path:", normalizePath(out_path, winslash = "/"), "\n")
  cat("Unique numeric Entrez IDs:", length(entrez_num), "\n")
  cat("Dropped non-numeric entries:", length(bad_non_numeric), "\n")
  if (length(bad_non_numeric) > 0) {
    cat("First 20 dropped (non-numeric):\n")
    print(head(sort(unique(bad_non_numeric)), 20))
  }
  
  invisible(entrez_num)
}

# ---- BACKGROUND ----
bg_entrez <- background_db_entrez$EntrezID
bg_file <- file.path(enrichr_dir, "BACKGROUND_entrez_ids.txt")
write_entrez_list(bg_entrez, bg_file, "BACKGROUND (elite enhancers overlapping H3K27ac CpGs)")

# ---- FOREGROUND 1 ----
fg1_entrez <- fg1_long_entrez$EntrezID
fg1_file <- file.path(enrichr_dir, "FOREGROUND1_entrez_ids.txt")
write_entrez_list(fg1_entrez, fg1_file, "FOREGROUND 1 (FG1 long annotation EntrezIDs)")

# ---- FOREGROUND 2 ----
fg2_entrez <- fg2_long_entrez$EntrezID
fg2_file <- file.path(enrichr_dir, "FOREGROUND2_entrez_ids.txt")
write_entrez_list(fg2_entrez, fg2_file, "FOREGROUND 2 (FG2 long annotation EntrezIDs)")

cat("\nDONE. EnrichR submission folder:\n", normalizePath(enrichr_dir, winslash = "/"), "\n")


# =============================================================================
# 6b) EnrichR submission files (RAW SYMBOLS): BACKGROUND + FG1 + FG2
#      - one symbol per line
#      - writes 3 txt files into the SAME enrichr_dir
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# reuse existing folder from previous step (create if missing)
enrichr_dir <- file.path(out_dir, "EnrichR_submission")
dir.create(enrichr_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helper: clean + write symbol list ----
write_symbol_list <- function(sym_vec, out_path, label) {
  sym_vec <- unique(trimws(as.character(sym_vec)))
  sym_vec <- sym_vec[!is.na(sym_vec) & sym_vec != ""]
  sym_vec <- sort(sym_vec)
  
  writeLines(sym_vec, con = out_path, useBytes = TRUE)
  
  cat("\n==================== EnrichR RAW SYMBOL file written ====================\n")
  cat("Set:", label, "\n")
  cat("Path:", normalizePath(out_path, winslash = "/"), "\n")
  cat("Unique symbols:", length(sym_vec), "\n")
  cat("First 20 symbols:\n")
  print(head(sym_vec, 20))
  
  invisible(sym_vec)
}

# ---- BACKGROUND (raw symbols) ----
bg_sym <- background_db_entrez$symbol
bg_sym_file <- file.path(enrichr_dir, "BACKGROUND_symbols_raw.txt")
write_symbol_list(bg_sym, bg_sym_file, "BACKGROUND (raw GeneHancer symbols)")

# ---- FOREGROUND 1 (raw symbols) ----
fg1_sym <- fg1_long_entrez$symbol
fg1_sym_file <- file.path(enrichr_dir, "FOREGROUND1_symbols_raw.txt")
write_symbol_list(fg1_sym, fg1_sym_file, "FOREGROUND 1 (raw GeneHancer symbols)")

# ---- FOREGROUND 2 (raw symbols) ----
fg2_sym <- fg2_long_entrez$symbol
fg2_sym_file <- file.path(enrichr_dir, "FOREGROUND2_symbols_raw.txt")
write_symbol_list(fg2_sym, fg2_sym_file, "FOREGROUND 2 (raw GeneHancer symbols)")

cat("\nDONE. RAW symbol EnrichR files created in:\n", normalizePath(enrichr_dir, winslash = "/"), "\n")

