#!/usr/bin/env Rscript
# Writes ref/gene_sets.csv (every gene set the engines use), ref/state_cutpoints.csv and
# www/downloads/state_cutpoints.csv (the CMI cut points of the grades) and ref/SHA256SUMS (checksums of
# the files in ref/). Run after changing ref/ or engine/constants.R:
#   Rscript tools/export_gene_sets.R .
args <- commandArgs(trailingOnly = TRUE)
app <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
for (f in c("constants.R", "io.R", "models.R")) source(file.path(app, "engine", f))
a <- load_assets(file.path(app, "ref"))
# LF line endings on every platform (a text-mode file on Windows would get CRLF and other checksums).
write_lf <- function(path, fun) { con <- file(path, open = "wb"); on.exit(close(con)); fun(con) }
row <- function(set, desc, genes) data.frame(set = set, description = desc, gene = genes, stringsAsFactors = FALSE)
gs <- rbind(
  row("Ig_constant_heavy", "Immunoglobulin heavy-chain constant regions (Ig-score metrics)", IG_G_HEAVY),
  row("Ig_constant_light", "Immunoglobulin light-chain constant regions (Ig-score metrics)", IG_G_LIGHT),
  row("Ig_J_chain", "J chain (secretory metrics)", "JCHAIN"),
  row("Ig_variable_total_burden", "Immunoglobulin variable-region genes counted in total Ig burden only", setdiff(IG_G_PANEL, c(IG_G_HEAVY, IG_G_LIGHT))),
  do.call(rbind, lapply(names(a$mtp$programs), function(p) row(paste0("MTP_", p), paste("Malignant transcriptional program:", gsub("_", " ", p)), a$mtp$programs[[p]]))),
  row("Cell_cycle_S", "S-phase genes (Seurat cc.genes.updated.2019)", a$ann$cell_cycle_genes$s.genes),
  row("Cell_cycle_G2M", "G2/M genes (Seurat cc.genes.updated.2019)", a$ann$cell_cycle_genes$g2m.genes),
  row("Upload_check_B_markers", "B-cell markers used only for the upload check", B_MARKERS),
  row("Upload_check_non_B_markers", "T-cell, myeloid, NK, erythroid and platelet markers used only for the upload check", NON_B_MARKERS))
write_lf(file.path(app, "ref", "gene_sets.csv"), function(con) utils::write.csv(gs, con, row.names = FALSE))
# CMI grades; a value on a boundary belongs to the lower grade.
num <- function(x) ifelse(is.na(x), "", sprintf("%.16g", x))
cut <- data.frame(grade = c("Below threshold", "Mal-1", "Mal-2", "Mal-3"),
                  cmi_from = num(c(NA, CMI_THRESHOLD, CMI_INDEX2, CMI_INDEX3)), cmi_to = num(c(CMI_THRESHOLD, CMI_INDEX2, CMI_INDEX3, NA)),
                  rule = c(sprintf("CMI < %s", num(CMI_THRESHOLD)), sprintf("%s <= CMI <= %s", num(CMI_THRESHOLD), num(CMI_INDEX2)),
                           sprintf("%s < CMI <= %s", num(CMI_INDEX2), num(CMI_INDEX3)), sprintf("CMI > %s", num(CMI_INDEX3))),
                  stringsAsFactors = FALSE)
for (f in c(file.path(app, "ref", "state_cutpoints.csv"), file.path(app, "www", "downloads", "state_cutpoints.csv")))
  write_lf(f, function(con) utils::write.csv(cut, con, row.names = FALSE))
files <- sort(setdiff(list.files(file.path(app, "ref")), "SHA256SUMS"))
sums <- vapply(files, function(f) digest::digest(file = file.path(app, "ref", f), algo = "sha256"), character(1))
write_lf(file.path(app, "ref", "SHA256SUMS"), function(con) writeLines(paste0(sums, "  ", files), con))
cat(nrow(gs), "genes in", length(unique(gs$set)), "sets;", length(files), "files checksummed\n")
