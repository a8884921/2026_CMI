.release_script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
.release_root <- normalizePath(file.path(.release_script_dir, "..", ".."))
source(file.path(.release_root, "privacy.R"), encoding = "UTF-8")
PANEL_DIR <- local({ a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd() })
PANEL_DIR <- normalizePath(PANEL_DIR)
source_release_dependency(file.path(.release_root, "shared", "figure3_config.R"))
suppressMessages({ library(readr); library(pROC); library(ggplot2) })
cells <- validate_anonymous_data(readr::read_csv(plot_csv("Figure3B_train_cell_CMI.csv"),
                         show_col_types = FALSE))
stopifnot(all(c("id", "GROUP", "CMI") %in% names(cells)))
r  <- pROC::roc(factor(cells$GROUP, levels = c("normal", "tumor")), cells$CMI,
                levels = c("normal", "tumor"), direction = "<", quiet = TRUE)
ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
auc <- as.numeric(pROC::auc(r))
n_tumor  <- sum(cells$GROUP == "tumor")
n_normal <- sum(cells$GROUP == "normal")
stopifnot(abs(auc - 0.9835) < 5e-4, abs(ci[1] - 0.9814) < 5e-4, abs(ci[3] - 0.9857) < 5e-4)
d <- data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities)
d <- d[order(d$fpr, d$tpr), ]
d <- d[!duplicated(d$fpr, fromLast = TRUE), ]
lab <- sprintf("AUC = %.3f\n95%% CI: %.3f - %.3f", auc, ci[1], ci[3])
pB <- ggplot(d, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = SEM[["neutral"]], linewidth = LW_AXIS) +
  geom_area(fill = "#FDD9B5", alpha = 0.55) +
  geom_path(colour = SEM[["case"]], linewidth = LW_DATA) +
  annotate("text", x = 0.05, y = 0.30, label = lab, hjust = 0, vjust = 0,
           size = PT(SZ_TICK), family = FIG_FONT, colour = "black") +
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = c(0, 0.4, 0.8),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = c(0, 0.4, 0.8),
                     expand = c(0.01, 0.01)) +
  coord_equal() +
  theme_fig3() +
  theme(plot.margin = ggplot2::margin(SZ_TICK * 0.6, SZ_TICK * 0.6,
                                      SZ_TICK * 0.4, SZ_TICK * 0.4))
pB <- tag_plot(pB, "B")
pB <- pB + ggplot2::theme(plot.tag.position = c(0, 0.96))
DS_B <- "B-cell lymphoma scRNA-seq \u00b7 internal training set \u00b7 n = 38,374 tumor"
pB <- dataset_line(pB, DS_B)
if (requireNamespace("systemfonts", quietly = TRUE)) {
  w_cap <- systemfonts::string_width(DS_B, family = FIG_FONT, size = SZ_DATA) * 0.3528
  w_lim <- 64
  invisible(NULL)
  stopifnot(w_cap <= w_lim)
}
save_panel(pB, "Figure3B_internal_train_CMI_ROC", 84, 76,
           dataset = DS_B)
