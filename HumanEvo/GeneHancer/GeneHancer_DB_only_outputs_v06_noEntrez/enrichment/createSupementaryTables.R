# =============================================================================
# GeneHancer enrichment export (S3 + S6) — Excel workbook glue (per folder)
# =============================================================================
#
# GOAL
# ----
# From the project folder:
#   GeneHancer_DB_only_outputs_v06_noEntrez/enrichment/
# create ONE Excel file per enrichment subfolder (e.g., S3, S6), by gluing ALL
# Enrichr output *.txt tables inside each subfolder into a single workbook
# (one sheet per txt file).
#
# INPUTS
# ------
# Root folder (script location):
#   .../GeneHancer_DB_only_outputs_v06_noEntrez/enrichment/
# Subfolders:
#   S3/, S6/  (and optionally any S* folders present)
#
# RULES
# -----
# - For each target subfolder (S3, S6):
#     - include EVERY *.txt file found directly inside it
#     - sheet name = filename without extension, with trailing "_table" removed
#     - sanitize sheet names for Excel and truncate to 31 chars
#     - resolve duplicates by adding a numeric suffix
# - Output file is written INSIDE each subfolder:
#     S3/Table S3 enrichment analysis.xlsx
#     S6/Table S6 enrichment analysis.xlsx
#
# OUTPUTS
# -------
# One workbook per folder: "Table <folder> enrichment analysis.xlsx"
#
# Author: Daniel Batyrev (HUJI 777634015)
# =============================================================================

# --------------------------- Set up environment ------------------------------
rm(list = ls())

cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
} else {
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
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
  library(openxlsx)
})

# ------------------------------ Parameters -----------------------------------

# (Optional alternative) automatically pick all S* folders:
target_folders <- list.dirs(this.dir, full.names = FALSE, recursive = FALSE)
target_folders <- target_folders[grepl("^S\\d+$", target_folders)]

# ------------------------------ Helpers --------------------------------------
make_sheet_name <- function(txt_file) {
  nm <- tools::file_path_sans_ext(basename(txt_file))
  nm <- gsub("_table$", "", nm)
  nm <- gsub("[\\[\\]\\*\\?:/\\\\]", "_", nm)  # invalid Excel sheet chars
  nm <- substr(nm, 1, 31)
  nm
}

safe_add_sheet <- function(wb, sheet_name) {
  if (!(sheet_name %in% names(wb))) return(sheet_name)
  i <- 2
  base <- substr(sheet_name, 1, 28)
  while (paste0(base, "_", i) %in% names(wb)) i <- i + 1
  paste0(base, "_", i)
}

# ------------------------------ Main -----------------------------------------
for (folder in target_folders) {
  folder_path <- file.path(this.dir, folder)
  
  if (!dir.exists(folder_path)) {
    message("⚠️ Skip (folder not found): ", folder_path)
    next
  }
  
  txt_files <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)
  
  if (length(txt_files) == 0) {
    message("⚠️ Skip (no .txt files): ", folder_path)
    next
  }
  
  out_xlsx <- file.path(folder_path, paste0("Table ", folder, " enrichment analysis.xlsx"))
  
  wb <- createWorkbook()
  
  for (txt_file in txt_files) {
    df <- read.delim(
      txt_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    sheet_name <- make_sheet_name(txt_file)
    sheet_name <- safe_add_sheet(wb, sheet_name)
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, x = df)
  }
  
  saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
  message("✅ Saved workbook: ", out_xlsx)
}
