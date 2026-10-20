.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- gsub("~+~", " ", sub("^--file=", "", a[grep("^--file=", a)]), fixed = TRUE)
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
PANEL_DIR <- normalizePath(.script_dir)
source_release_dependency(file.path(.release_root, "shared", "figure1_config.R"))
DATASETS_TXT <- "Analysis cohorts"
suppressMessages({ library(readr); library(dplyr); library(ggplot2) })
W_MM <- 60; H_MM <- 56
TRAD <- c(Ig     = "Immunoglobulin /\nBCR clonality",
          CNV    = "inferCNV /\nCopyKAT",
          MTP    = "Malignant gene\nprograms",
          Diff   = "Differentiation\nstate",
          Prolif = "Cell-cycle\nstate")
freq <- validate_anonymous_data(readr::read_csv(P("panelD_algorithm_frequency.csv"), show_col_types = FALSE))
UPSTREAM <- c(Ig     = "Immunoglobulin / BCR clonality",
              CNV    = "Copy-number variation (inferCNV / CopyKAT)",
              MTP    = "Malignant transcriptional program",
              Diff   = "Differentiation state",
              Prolif = "Proliferation (cell cycle)")
N_STUDIES <- 10
stopifnot(all(c("method", "full", "count", "pct") %in% names(freq)), !anyDuplicated(freq$method),
          setequal(freq$method, names(TRAD)), all(freq$full == UPSTREAM[freq$method]),
          all(freq$count >= 0 & freq$count <= N_STUDIES),
          isTRUE(all.equal(freq$pct, freq$count / N_STUDIES * 100)))
freq <- freq |>
  dplyr::mutate(trad  = unname(TRAD[method]),
                trad  = factor(trad, levels = trad[order(count, -match(method, names(TRAD)))]),
                label = sprintf("%d (%d%%)", count, pct),
                inside = count >= 5,
                lab_x  = ifelse(inside, count - 0.25, count + 0.25),
                lab_h  = ifelse(inside, 1, 0))
pD <- ggplot(freq, aes(count, trad)) +
  geom_col(fill = FREQ_COL, width = 0.66) +
  geom_text(aes(x = lab_x, label = label, hjust = lab_h), size = SZ_TICK / PT,
            family = FIG_FONT, colour = "#1A1A1A") +
  scale_x_continuous(sprintf("Prior studies (n = %d)", N_STUDIES), limits = c(0, N_STUDIES),
                     breaks = c(0, 5, 10), expand = c(0, 0)) +
  labs(y = NULL) +
  theme_fig1() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = SZ_TICK, colour = "black", lineheight = 0.9),
        plot.margin = ggplot2::margin(SZ_TICK * 1.6, SZ_TICK * 0.9, SZ_TICK * 0.3, SZ_TICK * 0.3))
pD <- PANEL_TAG(pD, "D")
pD <- dataset_line(pD, sprintf("Literature: %d curated studies", N_STUDIES))
save_panel(pD, W_MM, H_MM)
print(data.frame(method = freq$method, label = gsub("\n", " ", as.character(freq$trad)), count = freq$count,
                 pct = freq$pct)[order(-freq$count), ], row.names = FALSE)
