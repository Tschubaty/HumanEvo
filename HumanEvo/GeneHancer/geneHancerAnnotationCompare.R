# =============================================================================
# REBUILD Table S2 + Table S5 with CORRECT GeneHancer annotation (DB-only, RAW symbols)
#   - DOES NOT use any Entrez conversion files
#   - Re-annotates CpGs directly from:
#       GeneHancer_AnnotSV_hg19_v5.25.txt                 (regulatory elements coords)
#       GeneHancer_AnnotSV_gene_association_scores_v5.25.txt (GHid -> gene symbol + is_elite)
#   - Logic (the “correct” one we discussed):
#       1) CpG must fall INSIDE a GeneHancer regulatory element (coords table)
#       2) keep ONLY elite gene associations (is_elite==1) for those GHids
#       3) per CpG: collapse ALL GHids + ALL elite symbols (semicolon-separated)
#   - Prints: genes removed vs added (old table vs rebuilt)
#   - Writes: new Excel outputs into a NEW timestamped folder under:
#       C:\Users\Batyrev\Documents\GitHub\HumanEvo\HumanEvo\GeneHancer\
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(openxlsx)
  library(GenomicRanges)
})

# -------------------- PATHS --------------------
this.dir <- "C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/GeneHancer"
setwd(this.dir)

# Your original tables (use the same files you already use in Dropbox)
file_S2 <- "C:/Users/Batyrev/Dropbox/R42v/Table S2 significant CpG [1b].xlsx"
file_S5 <- "C:/Users/Batyrev/Dropbox/R42v/Table S5 significant CpG [2b].xlsx"

stopifnot(file.exists(file_S2), file.exists(file_S5))

# GeneHancer DB files (must be in this.dir; adjust if needed)
gh_coords_file <- file.path(this.dir, "GeneHancer_AnnotSV_hg19_v5.25.txt")
gh_gene_file   <- file.path(this.dir, "GeneHancer_AnnotSV_gene_association_scores_v5.25.txt")

stopifnot(file.exists(gh_coords_file), file.exists(gh_gene_file))

# Output folder (new)
out_dir_tables <- file.path(
  this.dir,
  paste0("Rebuilt_Tables_DBonly_RAW_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# -------------------- HELPERS --------------------
collapse_unique <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  paste(sort(x), collapse = ";")
}

split_semicolon_unique <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  out <- unique(unlist(strsplit(x, ";", fixed = TRUE)))
  out <- trimws(out)
  out <- out[!is.na(out) & out != ""]
  sort(unique(out))
}

extract_unique_genes_from_table <- function(dt, gene_col) {
  split_semicolon_unique(dt[[gene_col]])
}

print_gene_diffs <- function(old_genes, new_genes, label) {
  only_old <- setdiff(old_genes, new_genes)
  only_new <- setdiff(new_genes, old_genes)
  common   <- intersect(old_genes, new_genes)
  
  cat("\n============================================================\n")
  cat("GENE SET COMPARISON:", label, "\n")
  cat("------------------------------------------------------------\n")
  cat("Old unique genes:", length(old_genes), "\n")
  cat("New unique genes:", length(new_genes), "\n")
  cat("Common genes    :", length(common), "\n")
  cat("Removed (old-only):", length(only_old), "\n")
  if (length(only_old) > 0) print(only_old)
  cat("Added (new-only)  :", length(only_new), "\n")
  if (length(only_new) > 0) print(only_new)
}

write_pretty_xlsx <- function(dt, out_xlsx, sheet_name) {
  wb <- createWorkbook()
  addWorksheet(wb, sheet_name)
  
  writeDataTable(
    wb,
    sheet = sheet_name,
    x = as.data.frame(dt),
    tableStyle = "TableStyleMedium9",
    withFilter = TRUE
  )
  
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = 1:ncol(dt), widths = "auto")
  
  # Wrap corrected annotation cols
  wrap_cols <- which(names(dt) %in% c("GeneAnnotation_corrected","EnhancerAnnotation_corrected"))
  if (length(wrap_cols) > 0) {
    addStyle(
      wb, sheet = sheet_name,
      style = createStyle(wrapText = TRUE, valign = "top"),
      rows = 2:(nrow(dt) + 1),
      cols = wrap_cols,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  cat("\nWrote:", normalizePath(out_xlsx, winslash = "/"), "\n")
}

# -------------------- LOAD TABLES --------------------
# Sheet 1 in both, as in your prior code
fg1_old <- read_excel(file_S2, sheet = 1) |> as.data.table()
fg2_old <- read_excel(file_S5, sheet = 1) |> as.data.table()

stopifnot(all(c("chrom","start","end") %in% names(fg1_old)))
stopifnot(all(c("chrom","start","end") %in% names(fg2_old)))

# Identify peak-name columns (keep original column name for merging)
peak_col_fg1 <- if ("H3K27ac peak coordinates" %in% names(fg1_old)) "H3K27ac peak coordinates" else if ("name" %in% names(fg1_old)) "name" else NA_character_
peak_col_fg2 <- if ("name" %in% names(fg2_old)) "name" else if ("H3K27ac peak coordinates" %in% names(fg2_old)) "H3K27ac peak coordinates" else NA_character_

if (is.na(peak_col_fg1) || is.na(peak_col_fg2)) stop("Could not find the peak name column in one of the tables.")

# Add stable row ids (for joins that do not depend on duplicate coords)
fg1_old[, cpg_row := .I]
fg2_old[, cpg_row := .I]

# -------------------- LOAD GeneHancer DB --------------------
cat("\nLoading GeneHancer DB...\n")
gh_coords <- fread(gh_coords_file)  # chromosome/start/end/GHid
gh_gene   <- fread(gh_gene_file)    # GHid/symbol/combined_score/is_elite

# Standardize coords column name to chrom
setnames(gh_coords, "chromosome", "chrom")

# Keep ONLY elite gene associations (as agreed)
gh_gene_elite <- gh_gene[is_elite == 1]

cat("GeneHancer coords rows:", nrow(gh_coords), "\n")
cat("GeneHancer gene assoc rows (all):", nrow(gh_gene), "\n")
cat("GeneHancer elite gene assoc rows:", nrow(gh_gene_elite), "\n")
cat("Unique GHids in coords:", uniqueN(gh_coords$GHid), "\n")
cat("Unique GHids in elite assoc:", uniqueN(gh_gene_elite$GHid), "\n")

# -------------------- FUNCTION: annotate one table (fg1/fg2) --------------------
annotate_table_db_only <- function(fg_dt, peak_col_name) {
  
  # CpG GRanges
  cpg_gr <- GRanges(
    seqnames = fg_dt$chrom,
    ranges   = IRanges(start = as.integer(fg_dt$start), end = as.integer(fg_dt$end))
  )
  
  # GH coords GRanges (regulatory elements)
  gh_gr <- GRanges(
    seqnames = gh_coords$chrom,
    ranges   = IRanges(start = as.integer(gh_coords$start), end = as.integer(gh_coords$end)),
    GHid     = gh_coords$GHid
  )
  
  # Overlaps: CpG inside GH regulatory element
  ov <- findOverlaps(cpg_gr, gh_gr, ignore.strand = TRUE)
  
  # Long map: CpG row -> GHid (from coords)
  map_cpg_gh <- data.table(
    cpg_row = fg_dt$cpg_row[queryHits(ov)],
    GHid    = mcols(gh_gr)$GHid[subjectHits(ov)]
  )
  map_cpg_gh <- unique(map_cpg_gh)
  
  # Join elite genes: (CpG row, GHid) -> symbol
  # allow.cartesian=TRUE because a GHid maps to many symbols
  map_cpg_gh_gene <- merge(
    map_cpg_gh,
    gh_gene_elite[, .(GHid, symbol)],
    by = "GHid",
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  # Aggregate per CpG row
  per_cpg <- map_cpg_gh_gene[
    ,
    .(
      EnhancerAnnotation_corrected = collapse_unique(GHid),
      GeneAnnotation_corrected     = collapse_unique(symbol),
      n_GHids_corrected            = uniqueN(GHid, na.rm = TRUE),
      n_genes_corrected            = uniqueN(symbol, na.rm = TRUE)
    ),
    by = cpg_row
  ]
  
  # Merge back to fg_dt (KEEP ALL CpGs)
  out <- merge(
    fg_dt,
    per_cpg,
    by = "cpg_row",
    all.x = TRUE,
    sort = FALSE
  )
  
  # Fill NAs for non-annotated CpGs
  out[is.na(EnhancerAnnotation_corrected), EnhancerAnnotation_corrected := ""]
  out[is.na(GeneAnnotation_corrected),     GeneAnnotation_corrected := ""]
  out[is.na(n_GHids_corrected),            n_GHids_corrected := 0L]
  out[is.na(n_genes_corrected),            n_genes_corrected := 0L]
  
  # Coverage prints
  cat("\n==================== DB-only annotation coverage ====================\n")
  cat("Rows:", nrow(out), "\n")
  cat("CpGs with >=1 GHid:", sum(out$n_GHids_corrected > 0), "of", nrow(out), "\n")
  cat("Unique corrected GHids (table-level):", length(split_semicolon_unique(out$EnhancerAnnotation_corrected)), "\n")
  cat("Unique corrected genes (table-level):", length(split_semicolon_unique(out$GeneAnnotation_corrected)), "\n")
  
  out
}

# -------------------- RUN ANNOTATION --------------------
cat("\nAnnotating Table S2 (FG1) from DB...\n")
fg1_rebuilt <- annotate_table_db_only(fg1_old, peak_col_fg1)

cat("\nAnnotating Table S5 (FG2) from DB...\n")
fg2_rebuilt <- annotate_table_db_only(fg2_old, peak_col_fg2)

# -------------------- PRINT GENE DIFFS (old vs corrected) --------------------
# Old gene column is GeneAnnotation in your tables
old_fg1_genes <- extract_unique_genes_from_table(fg1_old, "GeneAnnotation")
old_fg2_genes <- extract_unique_genes_from_table(fg2_old, "GeneAnnotation")

new_fg1_genes <- extract_unique_genes_from_table(fg1_rebuilt, "GeneAnnotation_corrected")
new_fg2_genes <- extract_unique_genes_from_table(fg2_rebuilt, "GeneAnnotation_corrected")

print_gene_diffs(old_fg1_genes, new_fg1_genes, "Table S2 (FG1) old vs corrected(DB-only)")
print_gene_diffs(old_fg2_genes, new_fg2_genes, "Table S5 (FG2) old vs corrected(DB-only)")

# -------------------- OPTIONAL: also print enhancer diffs --------------------
old_fg1_enh <- split_semicolon_unique(fg1_old$EnhancerAnnotation)
old_fg2_enh <- split_semicolon_unique(fg2_old$EnhancerAnnotation)

new_fg1_enh <- split_semicolon_unique(fg1_rebuilt$EnhancerAnnotation_corrected)
new_fg2_enh <- split_semicolon_unique(fg2_rebuilt$EnhancerAnnotation_corrected)

cat("\n==================== ENHANCER (GHid) DIFFS ====================\n")
cat("FG1 old enh:", length(old_fg1_enh), " | new enh:", length(new_fg1_enh), "\n")
cat("FG2 old enh:", length(old_fg2_enh), " | new enh:", length(new_fg2_enh), "\n")

# -------------------- WRITE NEW EXCEL FILES --------------------
out_S2_new <- file.path(out_dir_tables, "Table_S2_significant_CpG_1b_CORRECTED_DBonly_RAW.xlsx")
out_S5_new <- file.path(out_dir_tables, "Table_S5_significant_CpG_2b_CORRECTED_DBonly_RAW.xlsx")

write_pretty_xlsx(fg1_rebuilt, out_S2_new, "S2_corrected_DBonly")
write_pretty_xlsx(fg2_rebuilt, out_S5_new, "S5_corrected_DBonly")

cat("\nDONE. Output folder:\n", normalizePath(out_dir_tables, winslash = "/"), "\n")




