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
#   BED3 files with CpGs
#  table S2 clalcutae from Analysis1_best_CpGs.bed 
#  and S5 from df_significant_Q_annot.bed 

# OUTPUTS
# -------
# (1) Background: GHids + elite genes + Enhancer coords + scores
# (2) Foreground: for each forground: excel with two tabs: 
#    wide:  CpG × GHid × elite genes
#    and long:  elite genes x 
# (3) Convenience gene lists for enrichment tools background forground (raw symbols)
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

# --- Foreground CpG lists:
# BED3 files 
tbl2 <- file.path(this.dir, "Analysis1_best_CpGs.bed")  # example (edit if needed)
tbl5 <- file.path(this.dir, "df_significant_Q_annot.bed")  # placeholder; set to your 2nd CpG list bed


# --- Output folder ---
out_dir <- file.path(this.dir, "GeneHancer_DB_only_outputs2")
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

