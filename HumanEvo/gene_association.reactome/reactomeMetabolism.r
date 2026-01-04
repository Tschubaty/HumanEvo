#####################################################################################
# =============================================================================
# Documentation
# =============================================================================
# Purpose
#   Create a reproducible gene list for "metabolic genes" using curated Reactome
#   pathway annotations, **in Entrez Gene ID space**.
#
#   Definition used here:
#     "Metabolic genes" = all genes associated with the Reactome pathway
#     "Metabolism" (R-HSA-1430728) AND all of its descendant pathways in the
#     Reactome pathway hierarchy.
#
#   This script ONLY creates the Reactome-based metabolic gene set and the
#   corresponding Reactome gene universe (both as Entrez Gene IDs).
#
#   A later step (NOT implemented here) will intersect this gene set with the
#   analysis-matched background/foreground gene universes (e.g. derived from
#   bone H3K27ac peaks + CpG→gene mapping), as requested by the reviewer.
#
# Inputs (expected in the same folder as this script)
#   1) NCBI2Reactome_All_Levels.txt
#        - Reactome mapping: Entrez Gene ID ↔ Reactome pathway membership.
#        - Columns:
#            1 EntrezID, 2 PathwayID, 3 URL, 4 PathwayName, 5 Evidence, 6 Species
#   2) ReactomePathways.txt
#        - Reactome pathway ID ↔ pathway name ↔ species.
#   3) ReactomePathwaysRelation.txt
#        - Reactome pathway hierarchy (parent↔child relations).
#
# Outputs (written to ./reactome_metabolism_gene_set/)
#   1) reactome_metabolism_summary.tsv
#        - Root pathway ID, number of descendant pathways, gene counts, fraction.
#   2) reactome_metabolism_pathway_ids.tsv
#        - All Reactome pathway IDs under "Metabolism" (human only).
#   3) reactome_metabolism_genes_entrez.tsv
#        - Unique Entrez Gene IDs under Metabolism (root + descendants).
#   4) reactome_all_genes_entrez_reactome_universe.tsv
#        - Unique Entrez Gene IDs present anywhere in Reactome (human only).
#
# Author
#   Daniel Batyrev (HUJI 777634015)
# =============================================================================


# =============================================================================
# Set up Work Environment
# =============================================================================

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


# =============================================================================
# Load Libraries
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})


# =============================================================================
# Inputs
# =============================================================================
file_map <- file.path(this.dir, "NCBI2Reactome_All_Levels.txt")
file_pw  <- file.path(this.dir, "ReactomePathways.txt")
file_rel <- file.path(this.dir, "ReactomePathwaysRelation.txt")

stopifnot(file.exists(file_map), file.exists(file_pw), file.exists(file_rel))


# =============================================================================
# Step 1: Read Reactome pathway name table (ID ↔ name ↔ species)
# =============================================================================
pw <- fread(
  file_pw,
  sep = "\t",
  header = FALSE,
  col.names = c("pathway_id", "pathway_name", "species"),
  quote = ""
)

# Keep only Homo sapiens pathways
pw_hsa_ids <- pw[species == "Homo sapiens", unique(pathway_id)]

# Identify the root pathway ID for "Metabolism"
root_id <- pw[species == "Homo sapiens" & pathway_name == "Metabolism", pathway_id]
stopifnot(length(root_id) == 1)


# =============================================================================
# Step 2: Read pathway hierarchy and get ALL descendants of Metabolism root
# =============================================================================
rel <- fread(
  file_rel,
  sep = "\t",
  header = FALSE,
  col.names = c("parent", "child"),
  quote = ""
)

# Restrict relations to human pathway IDs
rel_hsa <- rel[parent %in% pw_hsa_ids & child %in% pw_hsa_ids]

# Breadth-first traversal to collect all descendants of root_id
met_paths <- root_id
frontier <- root_id

repeat {
  kids <- unique(rel_hsa[parent %in% frontier, child])
  kids <- setdiff(kids, met_paths)
  if (length(kids) == 0) break
  met_paths <- c(met_paths, kids)
  frontier <- kids
}
met_paths <- unique(met_paths)


# =============================================================================
# Step 3: Read NCBI(Entrez)→Reactome mapping and extract metabolism genes
# =============================================================================
# Expected columns in NCBI2Reactome_All_Levels.txt:
#   1 EntrezID, 2 PathwayID, 3 URL, 4 PathwayName, 5 Evidence, 6 Species

map <- fread(
  file_map,
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE
)

setnames(
  map,
  c("entrez_id", "pathway_id", "url", "pathway_name", "evidence", "species")
)

# Keep only Homo sapiens
map_hsa <- map[species == "Homo sapiens"]

# Coerce Entrez IDs to integer safely (drop non-numeric if any)
map_hsa[, entrez_id := suppressWarnings(as.integer(entrez_id))]
map_hsa <- map_hsa[!is.na(entrez_id)]

# DENOMINATOR: all human Entrez genes that appear anywhere in Reactome
all_reactome_entrez <- unique(map_hsa[, .(entrez_id)])

# NUMERATOR: metabolic genes (subset to Metabolism root + descendants)
map_met <- map_hsa[pathway_id %in% met_paths]
met_entrez <- unique(map_met[, .(entrez_id)])

# fraction within Reactome universe
metabolic_fraction_in_reactome_universe <- nrow(met_entrez) / nrow(all_reactome_entrez)


# =============================================================================
# Outputs
# =============================================================================
out_dir <- file.path(this.dir, "reactome_metabolism_gene_set")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Summary
fwrite(
  data.table(
    metabolism_root_pathway_id = root_id,
    n_metabolism_descendant_pathways = length(met_paths),
    n_unique_metabolic_genes_entrez = nrow(met_entrez),
    n_unique_reactome_genes_entrez_universe = nrow(all_reactome_entrez),
    metabolic_fraction_in_reactome_universe = metabolic_fraction_in_reactome_universe
  ),
  file.path(out_dir, "reactome_metabolism_summary.tsv"),
  sep = "\t"
)

# Pathway IDs under Metabolism
fwrite(
  data.table(pathway_id = met_paths),
  file.path(out_dir, "reactome_metabolism_pathway_ids.tsv"),
  sep = "\t"
)

# Metabolic genes (Entrez IDs)
fwrite(
  met_entrez,
  file.path(out_dir, "reactome_metabolism_genes_entrez.tsv"),
  sep = "\t"
)

# Reactome universe (Entrez IDs)
fwrite(
  all_reactome_entrez,
  file.path(out_dir, "reactome_all_genes_entrez_reactome_universe.tsv"),
  sep = "\t"
)


# =============================================================================
# Done
# =============================================================================
cat("\nDone.\n")
cat("Reactome Metabolism root pathway:", root_id, "\n")
cat("Number of metabolism pathways (root + descendants):", length(met_paths), "\n")
cat("Number of unique metabolic genes (Entrez IDs):", nrow(met_entrez), "\n")
cat("Number of unique human genes in Reactome universe (Entrez IDs):",
    nrow(all_reactome_entrez), "\n")
cat("Metabolic fraction in Reactome universe:",
    round(metabolic_fraction_in_reactome_universe, 4), "\n")
cat("Output folder:", out_dir, "\n\n")

#####################################################################################
#####################################################################################
# Reactome Entrez composition plot (canonical vs provisional)
#####################################################################################
library(data.table)
# # all_reactome_entrez and met_entrez already exist as data.tables with column: entrez_id
# 
# # ---- Build reactome_summary (THIS WAS MISSING) ----
# get_entrez_split <- function(entrez_dt) {
#   x <- unique(entrez_dt$entrez_id)
#   x <- x[!is.na(x)]
#   data.table(
#     canonical_entrez   = sum(x <  1e6),
#     provisional_entrez = sum(x >= 1e6),
#     total_entrez       = length(x)
#   )
# }
# 
# reactome_summary <- rbindlist(list(
#   cbind(data.table(set = "Reactome\n(all pathways)"),  get_entrez_split(all_reactome_entrez)),
#   cbind(data.table(set = "Reactome\n(metabolism)"),    get_entrez_split(met_entrez))
# ), use.names = TRUE, fill = TRUE)
# 
# # ---- Prepare matrix ----
# mat_reactome <- rbind(
#   "Mapped Entrez IDs"       = reactome_summary$canonical_entrez,
#   "Provisional Entrez IDs"  = reactome_summary$provisional_entrez
# )
# colnames(mat_reactome) <- reactome_summary$set
# 
# cols_reactome <- c(
#   "Mapped Entrez IDs"      = "#3B82F6",
#   "Provisional Entrez IDs" = "#F59E0B"
# )
# 
#
# # ---- Plot ----
# op <- par(no.readonly = TRUE)
# par(mar = c(7, 5, 5, 2), xpd = FALSE)
# 
# totals <- colSums(mat_reactome)
# 
# bp <- barplot(
#   mat_reactome,
#   beside = FALSE,
#   col = cols_reactome[rownames(mat_reactome)],
#   border = NA,
#   ylim = c(0, max(totals) * 1.15),
#   ylab = "Number of Entrez genes",
#   main = "Reactome gene universe\nEntrez ID composition",
#   names.arg = rep("", ncol(mat_reactome))
# )
# 
# # ---- X-axis labels (fixed) ----
# axis(
#   side = 1,
#   at = bp,
#   labels = c(
#     "Reactome\n(all pathways)",
#     "Reactome\n(metabolism)"
#   ),
#   tick = FALSE,
#   las = 1,
#   line = 0
# )
# 
# # ---- Totals on top (FIXED) ----
# text(
#   x = bp,
#   y = totals + 0.04 * max(totals),
#   labels = format(totals, big.mark = ","),
#   cex = 0.9,
#   font = 2
# )
# 
# legend(
#   "topright",
#   inset = 0.02,
#   legend = rownames(mat_reactome),
#   fill = cols_reactome[rownames(mat_reactome)],
#   bty = "n",
#   cex = 0.9
# )
# 
# par(op)


library(data.table)

# all_reactome_entrez and met_entrez already exist as data.tables with column: entrez_id

# =========================
# Build reactome_summary
# =========================
get_entrez_split <- function(entrez_dt) {
  x <- unique(entrez_dt$entrez_id)
  x <- x[!is.na(x)]
  data.table(
    canonical_entrez   = sum(x <  1e6),
    provisional_entrez = sum(x >= 1e6),
    total_entrez       = length(x)
  )
}

reactome_summary <- rbindlist(list(
  cbind(data.table(set = "Reactome\n(all pathways)"), get_entrez_split(all_reactome_entrez)),
  cbind(data.table(set = "Reactome\n(metabolism)"),   get_entrez_split(met_entrez))
), use.names = TRUE, fill = TRUE)

# =========================
# TOTAL panel: unmapped = 0 (Reactome is already Entrez space)
# (no "not in metabolism" category)
# =========================
reactome_summary[, `:=`(
  input_symbols   = total_entrez,
  mapped_lt_1e6   = canonical_entrez,
  mapped_ge_1e6   = provisional_entrez,
  mapped_total    = canonical_entrez + provisional_entrez,
  unmapped_unique = 0L
)]

# =========================
# Build matrices (same style as other plots)
# =========================
mat_total <- rbind(
  "Mapped Entrez IDs"       = reactome_summary$mapped_lt_1e6,
  "Provisional Entrez IDs"  = reactome_summary$mapped_ge_1e6,
  "Unmapped"                = reactome_summary$unmapped_unique
)

mat_mapped <- rbind(
  "Mapped Entrez IDs"       = reactome_summary$mapped_lt_1e6,
  "Provisional Entrez IDs"  = reactome_summary$mapped_ge_1e6
)

colnames(mat_total)  <- reactome_summary$set
colnames(mat_mapped) <- reactome_summary$set

cols_total  <- c("Mapped Entrez IDs"="#3B82F6", "Provisional Entrez IDs"="#F59E0B", "Unmapped"="#D1D5DB")
cols_mapped <- c("Mapped Entrez IDs"="#3B82F6", "Provisional Entrez IDs"="#F59E0B")

# =========================
# Save (PDF + high-res PNG)
# =========================
out_dir <- file.path(getwd(), "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pdf_file <- file.path(out_dir, "Reactome_Entrez_mapping_quality_two_panel.pdf")
png_file <- file.path(out_dir, "Reactome_Entrez_mapping_quality_two_panel.png")

draw_two_panel_reactome <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mfrow = c(1,2), mar = c(7,5,4,1))
  
  # ---------- Panel 1: TOTAL ----------
  tot_total <- colSums(mat_total)
  bp1 <- barplot(
    mat_total,
    beside = FALSE,
    col = cols_total[rownames(mat_total)],
    border = NA,
    names.arg = rep("", ncol(mat_total)),
    ylab = "Number of Entrez genes",
    main = "TOTAL input symbols\n(mapped + unmapped)",
    ylim = c(0, max(tot_total) * 1.10)
  )
  axis(
    side = 1,
    at = bp1,
    labels = colnames(mat_total),
    tick = FALSE,
    las = 1,
    line = 0
  )
  ymax1 <- max(tot_total)
  text(
    x = bp1,
    y = tot_total + 0.03 * ymax1,
    labels = format(tot_total, big.mark = ","),
    cex = 0.9,
    font = 2
  )
  
  mapped_pct <- round(100 * reactome_summary$mapped_total / pmax(1, reactome_summary$input_symbols), 1)
  mtext(paste0("mapped: ", mapped_pct, "%"),
        side = 1, at = bp1, line = 4, cex = 0.85)
  
  legend(
    "topright",
    inset = c(0, 0.02),
    legend = rownames(mat_total),
    fill = cols_total[rownames(mat_total)],
    bty = "n",
    cex = 0.9
  )
  
  # ---------- Panel 2: MAPPED-only ----------
  tot_mapped <- colSums(mat_mapped)
  bp2 <- barplot(
    mat_mapped,
    beside = FALSE,
    col = cols_mapped[rownames(mat_mapped)],
    border = NA,
    names.arg = rep("", ncol(mat_mapped)),
    ylab = "Number of Entrez genes",
    main = "MAPPED-only\n(Entrez IDs)",
    ylim = c(0, max(tot_mapped) * 1.10)
  )
  axis(
    side = 1,
    at = bp2,
    labels = colnames(mat_mapped),
    tick = FALSE,
    las = 1,
    line = 0
  )
  ymax2 <- max(tot_mapped)
  text(
    x = bp2,
    y = tot_mapped + 0.03 * ymax2,
    labels = format(tot_mapped, big.mark = ","),
    cex = 0.9,
    font = 2
  )
  
  legend(
    "topright",
    inset = c(0, 0.02),
    legend = rownames(mat_mapped),
    fill = cols_mapped[rownames(mat_mapped)],
    bty = "n",
    cex = 0.9
  )
}

pdf(pdf_file, width = 14, height = 6)
draw_two_panel_reactome()
dev.off()

png(png_file, width = 4200, height = 1800, res = 300)
draw_two_panel_reactome()
dev.off()

cat("\nSaved:\n", pdf_file, "\n", png_file, "\n")



## =============================================================================
## 
## Purpose: 
## Inputs used here:

## Notes:

## =============================================================================

# ---- Paths (edit if needed) ----

dir_in <- file.path("C:/Users/Batyrev/Documents/GitHub/HumanEvo/HumanEvo/GeneHancer/", "Entrez_sets_GeneHancer")

GH_bg_entrez <- readRDS(file.path(dir_in, "03_H3K27ac_CpG_overlap_background_entrez_unique_numeric.rds"))
GH_elite_entrez <- readRDS(file.path(dir_in, "02_GeneHancer_double_elite_entrez_unique_numeric.rds"))
GH_all_entrez <- readRDS(file.path(dir_in, "01_GeneHancer_all_symbols_entrez_unique_numeric.rds"))



# Reactome
reactome_all <- all_reactome_entrez$entrez_id
reactome_met <- met_entrez$entrez_id

# GeneHancer
GH_all <- GH_all_entrez
GH_elite <- GH_elite_entrez
GH_bg <- GH_bg_entrez


# Intersections with metabolism
bg_metabolic <- intersect(GH_bg, reactome_met)
elite_metabolic <- intersect(GH_elite, reactome_met)
allGH_metabolic <- intersect(GH_all, reactome_met)

fractions <- data.frame(
  set = c(
    "Reactome universe",
    "GeneHancer all",
    "GeneHancer double-elite",
    "H3K27ac CpG background"
  ),
  total_genes = c(
    length(reactome_all),
    length(GH_all),
    length(GH_elite),
    length(GH_bg)
  ),
  metabolic_genes = c(
    length(reactome_met),
    length(allGH_metabolic),
    length(elite_metabolic),
    length(bg_metabolic)
  )
)

fractions$metabolic_fraction <- round(
  fractions$metabolic_genes / fractions$total_genes,
  4
)

print(fractions)


expected_bg_metabolic <- length(GH_bg) * (length(reactome_met) / length(reactome_all))

cat("\nExpected metabolic genes in background (null):",
    round(expected_bg_metabolic, 1), "\n")
cat("Observed metabolic genes in background:",
    length(bg_metabolic), "\n")
cat("Fold enrichment:",
    round(length(bg_metabolic) / expected_bg_metabolic, 2), "\n")


#######################################################################################

########################## conversion in the Entrez universe ###################

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

to_entrez <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[x != "" & !is.na(x)]
  
  out <- data.frame(
    input = x,
    ENTREZID = NA_character_,
    keytype = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # helper: safe select (won't crash on invalid keys)
  safe_select <- function(keys, keytype) {
    tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = keys,
        keytype = keytype,
        columns = "ENTREZID"
      ),
      error = function(e) NULL
    )
  }
  
  # 1) ENSEMBL (ENSG...) — normalize ENSG IDs first (remove version suffix etc.)
  is_ensg <- grepl("^ENSG\\d+", x)
  if (any(is_ensg)) {
    ensg_raw <- x[is_ensg]
    ensg_norm <- sub("^(ENSG\\d+).*", "\\1", ensg_raw)   # keep only ENSG + digits
    ensg_norm <- unique(ensg_norm)
    
    m <- safe_select(ensg_norm, "ENSEMBL")
    if (!is.null(m)) {
      m <- m[!is.na(m$ENTREZID), c("ENSEMBL", "ENTREZID")]
      # match using normalized IDs
      idx <- match(sub("^(ENSG\\d+).*", "\\1", out$input), m$ENSEMBL)
      hit <- which(!is.na(idx))
      out$ENTREZID[hit] <- as.character(m$ENTREZID[idx[hit]])
      out$keytype[hit]  <- "ENSEMBL"
    }
  }
  
  # 2) SYMBOL
  unm <- is.na(out$ENTREZID)
  if (any(unm)) {
    keys <- unique(out$input[unm])
    m <- safe_select(keys, "SYMBOL")
    if (!is.null(m)) {
      m <- m[!is.na(m$ENTREZID), c("SYMBOL", "ENTREZID")]
      idx <- match(out$input, m$SYMBOL)
      hit <- which(unm & !is.na(idx))
      out$ENTREZID[hit] <- as.character(m$ENTREZID[idx[hit]])
      out$keytype[hit]  <- "SYMBOL"
    }
  }
  
  # 3) ALIAS
  unm <- is.na(out$ENTREZID)
  if (any(unm)) {
    keys <- unique(out$input[unm])
    m <- safe_select(keys, "ALIAS")
    if (!is.null(m)) {
      m <- m[!is.na(m$ENTREZID), c("ALIAS", "ENTREZID")]
      idx <- match(out$input, m$ALIAS)
      hit <- which(unm & !is.na(idx))
      out$ENTREZID[hit] <- as.character(m$ENTREZID[idx[hit]])
      out$keytype[hit]  <- "ALIAS"
    }
  }
  
  out
}


library(readxl)

# ---- Paths ----
file_S2 <- "C:/Users/Batyrev/Dropbox/R42v/Table S2 significant CpG [1b].xlsx"
file_S5 <- "C:/Users/Batyrev/Dropbox/R42v/Table S5 significant CpG [2b].xlsx"

# ---- Analysis 1 foreground ----
fg1_symbols <- read_excel(file_S2, sheet = 1) |>
  dplyr::pull(GeneAnnotation) |>
  unique()

# ---- Analysis 2 foreground ----
fg2_symbols <- read_excel(file_S5, sheet = 1) |>
  dplyr::pull(GeneAnnotation) |>
  unique()

# ---- CRISPR tested genes ----
crispr_symbols <- read_excel(file_S5, sheet = 2) |>
  dplyr::pull("Functional validation of candidate genes using a pooled CRISPR knockout screen in human hepatocytes") |>
  as.character() |>
  trimws() |>
  (\(x) x[!is.na(x) & x != ""])() |>
  unique()

#
fg1_entrez_map <- to_entrez(fg1_symbols)
fg2_entrez_map <- to_entrez(fg2_symbols)
crispr_entrez_map <- to_entrez(crispr_symbols)

fg1_entrez <- unique(as.numeric(fg1_entrez_map$ENTREZID))
fg2_entrez <- unique(as.numeric(fg2_entrez_map$ENTREZID))
crispr_entrez <- unique(as.numeric(crispr_entrez_map$ENTREZID))

fg1_entrez <- fg1_entrez[!is.na(fg1_entrez)]
fg2_entrez <- fg2_entrez[!is.na(fg2_entrez)]
crispr_entrez <- crispr_entrez[!is.na(crispr_entrez)]


summarize_fg <- function(name, symbols, entrez_map, entrez_vec) {
  cat("\n====================================================\n")
  cat("Foreground:", name, "\n")
  cat("----------------------------------------------------\n")
  cat("Input gene symbols:", length(unique(symbols)), "\n")
  cat("Mapped to Entrez:", sum(!is.na(entrez_map$ENTREZID)), "\n")
  cat("Unique Entrez IDs:", length(entrez_vec), "\n")
  cat("Mapping rate:",
      round(100 * length(entrez_vec) / length(unique(symbols)), 1), "%\n")
}

summarize_fg("Analysis 1 (CpGs 1b)", fg1_symbols, fg1_entrez_map, fg1_entrez)
summarize_fg("Analysis 2 (CpGs 2b)", fg2_symbols, fg2_entrez_map, fg2_entrez)
summarize_fg("CRISPR screen", crispr_symbols, crispr_entrez_map, crispr_entrez)

library(data.table)


# -------------------------------------------------------------------------
# INPUT: expects these mapping data.frames exist:
#   fg1_entrez_map, fg2_entrez_map, crispr_entrez_map
# each has at least: input, ENTREZID (from your to_entrez())
# -------------------------------------------------------------------------

summarize_map_for_plot <- function(set_name, map_df) {
  # input symbols used (unique)
  input_syms <- unique(trimws(as.character(map_df$input)))
  input_syms <- input_syms[!is.na(input_syms) & input_syms != ""]
  n_input <- length(input_syms)
  
  # mapped ENTREZIDs (may have NA)
  entrez_raw <- map_df$ENTREZID
  entrez_raw <- entrez_raw[!is.na(entrez_raw) & entrez_raw != ""]
  entrez_num <- suppressWarnings(as.numeric(entrez_raw))
  entrez_num <- entrez_num[!is.na(entrez_num)]
  entrez_unique <- unique(entrez_num)
  
  mapped_lt_1e6 <- sum(entrez_unique < 1e6)
  mapped_ge_1e6 <- sum(entrez_unique >= 1e6)
  
  # unmapped inputs = how many input symbols have no ENTREZID at all
  tmp <- data.table(
    input   = trimws(as.character(map_df$input)),
    ENTREZID = map_df$ENTREZID
  )
  tmp <- tmp[!is.na(input) & input != ""]
  
  tmp[, has_map := any(!is.na(ENTREZID) & ENTREZID != ""), by = input]
  
  # FIX: reference column properly inside i
  unmapped_inputs <- tmp[has_map == FALSE, uniqueN(input)]
  
  data.table(
    set = set_name,
    input_symbols = n_input,
    unmapped_inputs = unmapped_inputs,
    mapped_lt_1e6 = mapped_lt_1e6,
    mapped_ge_1e6 = mapped_ge_1e6
  )
}


# ---- build summary_dt for the 3 FOREGROUND sets ----
summary_dt <- rbindlist(list(
  summarize_map_for_plot("CpGs 1b\n(foreground)", fg1_entrez_map),
  summarize_map_for_plot("CpGs 2b\n(foreground)", fg2_entrez_map),
  summarize_map_for_plot("CRISPR\n(test set)",    crispr_entrez_map)
), use.names = TRUE, fill = TRUE)

print(summary_dt)

# =============================================================================
# TWO-PANEL PLOT (TOTAL + MAPPED-only)  — same style as your GeneHancer plot
# =============================================================================

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
mat_total <- rbind(
  "Mapped Entrez IDs"        = plot_dt$mapped_lt_1e6,
  "Provisional Entrez IDs"   = plot_dt$mapped_ge_1e6,
  "Unmapped"                 = plot_dt$unmapped_unique
)

mat_mapped <- rbind(
  "Mapped Entrez IDs"        = plot_dt$mapped_lt_1e6,
  "Provisional Entrez IDs"   = plot_dt$mapped_ge_1e6
)

colnames(mat_total)  <- plot_dt$set
colnames(mat_mapped) <- plot_dt$set

cols_total  <- c("Mapped Entrez IDs"="#3B82F6", "Provisional Entrez IDs"="#F59E0B", "Unmapped"="#D1D5DB")
cols_mapped <- c("Mapped Entrez IDs"="#3B82F6", "Provisional Entrez IDs"="#F59E0B")

# ---- save paths ----
out_dir <- file.path(getwd(), "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf_file <- file.path(out_dir, "Foreground_Entrez_mapping_TOTAL_vs_MAPPED.pdf")
png_file <- file.path(out_dir, "Foreground_Entrez_mapping_TOTAL_vs_MAPPED.png")

# ---- function to draw the plot (so we can save pdf+png) ----
draw_two_panel <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(mfrow = c(1,2), mar = c(9,5,4,1))
  
  # ========== Panel 1: TOTAL ==========
  tot_total <- colSums(mat_total)
  bp1 <- barplot(
    mat_total,
    beside = FALSE,
    col = cols_total[rownames(mat_total)],
    border = NA,
    names.arg = rep("", ncol(mat_total)),
    ylab = "Count",
    main = "TOTAL input symbols\n(mapped + unmapped)",
    ylim = c(0, max(tot_total) * 1.10)
  )
  axis(
    side = 1,
    at = bp1,
    labels = colnames(mat_total),
    tick = FALSE,
    line = 0,
    las = 1
  )
  ymax1 <- max(tot_total)
  text(
    x = bp1,
    y = tot_total + 0.03 * ymax1,
    labels = format(tot_total, big.mark = ","),
    cex = 0.9
  )
  mapped_pct <- round(100 * plot_dt$mapped_total / pmax(1, plot_dt$input_symbols), 1)
  mtext(paste0("mapped: ", mapped_pct, "%"),
        side = 1, at = bp1, line = 6, cex = 0.85)
  
  legend(
    "topright",
    inset = c(0, 0.02),
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
    names.arg = rep("", ncol(mat_mapped)),
    ylab = "Count",
    main = "MAPPED-only symbols\n(Entrez IDs)",
    ylim = c(0, max(tot_mapped) * 1.10)
  )
  axis(
    side = 1,
    at = bp2,
    labels = colnames(mat_mapped),
    tick = FALSE,
    line = 0,
    las = 1
  )
  ymax2 <- max(tot_mapped)
  text(
    x = bp2,
    y = tot_mapped + 0.03 * ymax2,
    labels = format(tot_mapped, big.mark = ","),
    cex = 0.9
  )
  
  legend(
    "topright",
    inset = c(0, 0.02),
    legend = rownames(mat_mapped),
    fill = cols_mapped[rownames(mat_mapped)],
    bty = "n",
    cex = 0.9
  )
}

# ---- save PDF ----
pdf(pdf_file, width = 16, height = 7)
draw_two_panel()
dev.off()

# ---- save PNG (high resolution) ----
png(png_file, width = 4800, height = 2100, res = 300)
draw_two_panel()
dev.off()

cat("\nSaved:\n", pdf_file, "\n", png_file, "\n")


################################# Metabolic enrichment vs analysis-matched background
# Key reviewer-facing statistical test
#
# H0: foreground metabolic fraction ≤ background metabolic fraction
# H1: foreground metabolic fraction > background metabolic fraction
#
# Universe: GH_bg_entrez (H3K27ac-matched background)
# Metabolic definition: Reactome Metabolism (root + descendants)

library(data.table)

# --------------------------------------------------
# 1) Define analysis-matched universe
# --------------------------------------------------
universe <- unique(as.integer(GH_bg_entrez))
universe <- universe[!is.na(universe)]

# Restrict Reactome metabolism to universe
reactome_met_u <- intersect(
  universe,
  unique(as.integer(reactome_met))
)
reactome_met_u <- reactome_met_u[!is.na(reactome_met_u)]

U <- length(universe)
M <- length(reactome_met_u)
bg_rate <- M / U

cat("\nAnalysis-matched background (GH_bg_entrez)\n")
cat("Universe size:", U, "\n")
cat("Metabolic genes in universe:", M,
    sprintf("(%.4f)\n", bg_rate))

# --------------------------------------------------
# 2) One-sided Fisher enrichment test
# --------------------------------------------------
fg_enrichment_test <- function(set_name, fg_vec, universe, reactome_met_u) {
  
  fg <- unique(as.integer(fg_vec))
  fg <- fg[!is.na(fg)]
  
  # IMPORTANT: restrict foreground to universe
  fg <- intersect(fg, universe)
  
  n <- length(fg)
  k <- length(intersect(fg, reactome_met_u))
  
  # 2x2 contingency table
  #                metabolic   non-metabolic
  # foreground        k         n-k
  # background      M-k       (U-M)-(n-k)
  
  mat <- matrix(
    c(k, n - k,
      M - k, (U - M) - (n - k)),
    nrow = 2,
    byrow = TRUE
  )
  
  dimnames(mat) <- list(
    group  = c("foreground", "background_rest"),
    status = c("metabolic", "non_metabolic")
  )
  
  ft <- fisher.test(mat, alternative = "greater")
  
  fg_rate <- if (n == 0) NA_real_ else k / n
  expected <- n * (M / U)
  fold_enrichment <- ifelse(expected == 0, NA_real_, k / expected)
  
  data.table(
    set = set_name,
    fg_genes_in_universe = n,
    fg_metabolic = k,
    fg_fraction = fg_rate,
    bg_fraction = M / U,
    expected_metabolic = expected,
    fold_enrichment = fold_enrichment,
    odds_ratio = unname(ft$estimate),
    p_enrichment = ft$p.value
  )
}

# --------------------------------------------------
# 3) Run ONLY the two CpG foreground sets
# --------------------------------------------------
res <- rbindlist(list(
  fg_enrichment_test("CpGs 1b", fg1_entrez, universe, reactome_met_u),
  fg_enrichment_test("CpGs 2b", fg2_entrez, universe, reactome_met_u)
), use.names = TRUE, fill = TRUE)

# --------------------------------------------------
# 4) Pretty output
# --------------------------------------------------
res_print <- copy(res)
res_print[, `:=`(
  fg_fraction = round(fg_fraction, 4),
  bg_fraction = round(bg_fraction, 4),
  expected_metabolic = round(expected_metabolic, 2),
  fold_enrichment = round(fold_enrichment, 2),
  odds_ratio = round(odds_ratio, 3),
  p_enrichment = signif(p_enrichment, 3)
)]

cat("\nMetabolic enrichment vs H3K27ac-matched background (one-sided Fisher)\n")
print(res_print)

# Optional save
fwrite(
  res,
  "metabolic_enrichment_fisher_vs_GH_bg_entrez_foregrounds_only.tsv",
  sep = "\t"
)
