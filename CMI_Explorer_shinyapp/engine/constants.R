# Fixed definitions shared by the engines and the interface. Do not edit without re-running
# tests/test_equivalence.R.

CMI_VERSION <- "2.0.0"
CMI_RELEASE_DATE <- "2026-10-18"

# ---- Immunoglobulin panel -------------------------------------------------------------------------
IG_MIN_UMI <- 3
IG_G_HEAVY <- c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE")
IG_G_IGG <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4")
IG_G_SWIT <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2")
IG_G_LIGHT <- c("IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC7")
IG_G_PANEL <- unique(c(IG_G_HEAVY, IG_G_LIGHT, "IGHV3-30", "IGHV4-34", "IGHV1-18", "IGHV5-78",
                       "IGKV4-1", "IGLV1-40", "IGLV6-57"))
IG_LOCUS <- c(IGHM = 0, IGHD = 0, IGHG3 = 1, IGHG1 = 2, IGHA1 = 3, IGHG2 = 4, IGHG4 = 5, IGHA2 = 7)
IG_CANON <- unique(c(IG_G_PANEL, IG_G_HEAVY, IG_G_IGG, IG_G_SWIT, "IGHE", names(IG_LOCUS), "JCHAIN", "IGJ"))

# ---- CMI grades -----------------------------------------------------------------------------------
# Malignancy threshold and the Mal-1/Mal-2 and Mal-2/Mal-3 boundaries on the CMI scale.
# A value on a boundary belongs to the lower grade; a value equal to the threshold is malignant.
CMI_THRESHOLD <- 0.9755057544904363
CMI_INDEX2 <- 0.98926105                 # Mal-1 / Mal-2 boundary
CMI_INDEX3 <- 0.99619844                 # Mal-2 / Mal-3 boundary
# Lineage markers for the upload check (not part of the CMI).
B_MARKERS <- c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5")
NON_B_MARKERS <- c("CD3E", "CD3D", "LYZ", "CD14", "NKG7", "HBB", "PPBP")
STATE_LEVELS <- c("Below threshold", "Mal-1", "Mal-2", "Mal-3")
assign_state <- function(cmi) {
  out <- ifelse(cmi < CMI_THRESHOLD, "Below threshold",
         ifelse(cmi <= CMI_INDEX2, "Mal-1", ifelse(cmi <= CMI_INDEX3, "Mal-2", "Mal-3")))
  factor(out, levels = STATE_LEVELS)
}
