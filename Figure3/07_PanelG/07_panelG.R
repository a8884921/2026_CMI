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
suppressMessages({ library(readr); library(ggplot2) })
wf <- validate_anonymous_data(readr::read_csv(plot_csv("Figure3G_waterfall_values.csv"), show_col_types = FALSE))
stopifnot(all(c("feature", "feature_value", "shap_value",
                "baseline_efx", "prediction_fx") %in% names(wf)))
col_shap <- if ("shap_value_as_published" %in% names(wf)) "shap_value_as_published" else "shap_value"
col_feat <- if ("feature_as_published" %in% names(wf)) "feature_as_published" else "feature"
wf$shap_value <- wf[[col_shap]]
wf$feature    <- wf[[col_feat]]
stopifnot(all(wf$feature %in% SCORE_LEVELS))
wf <- wf[order(-abs(wf$shap_value)), ]
stopifnot(length(unique(wf$baseline_efx)) == 1, length(unique(wf$prediction_fx)) == 1)
E_raw  <- unique(wf$baseline_efx); fx_raw <- unique(wf$prediction_fx)
E <- 0.0624; fx <- 2.61
stopifnot(abs(E - E_raw) < 1e-3, abs(fx - fx_raw) < 5e-3)
wf$cum_end   <- E + cumsum(wf$shap_value)
wf$cum_start <- c(E, head(wf$cum_end, -1))
wf$idx       <- seq_len(nrow(wf))
wf$xmin <- wf$idx - 0.42
wf$xmax <- wf$idx + 0.42
wf$ymin <- pmin(wf$cum_start, wf$cum_end)
wf$ymax <- pmax(wf$cum_start, wf$cum_end)
wf$mid  <- (wf$cum_start + wf$cum_end) / 2
wf$feature_label <- sprintf("%s=%s", SCORE_LABELS[wf$feature],
                            formatC(wf$feature_value, format = "g", digits = 4))
ylim    <- range(c(E, fx, wf$cum_start, wf$cum_end))
pad     <- diff(ylim) * 0.30
label_y <- max(c(E, fx, wf$cum_start, wf$cum_end)) + pad * 0.85
wf$thick <- (wf$ymax - wf$ymin) > pad * 0.22
wf$lab_x <- (wf$xmin + wf$xmax) / 2
MM_W   <- 170; MM_H <- 42
MM_PLOT <- MM_H - (SZ_TICK * 0.5 + SZ_TICK * 0.4) - (SZ_TICK * 0.5 + SZ_TICK * 0.4)
mm2y   <- diff(ylim) / MM_PLOT
h_lab  <- (SZ_TICK * 1.15) * 0.3528 * mm2y
wf$lab_y <- ifelse(wf$thick, wf$mid, wf$ymax + h_lab * 0.45)
wf$lab_v <- ifelse(wf$thick, 0.5, 1)
wf$top_y <- ifelse(wf$thick, wf$lab_y - h_lab * 0.45, wf$ymax + h_lab * 0.15)
pG <- ggplot(wf) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = shap_value > 0), colour = "grey20", linewidth = LW_AXIS) +
  geom_hline(yintercept = E, linetype = "dashed", colour = SEM[["neutral"]],
             linewidth = LW_ERR) +
  geom_label(aes(x = lab_x, y = lab_y, label = sprintf("%+.3g", shap_value), vjust = lab_v),
             size = PT(SZ_TICK), family = FIG_FONT, colour = "black",
             fill = "white", label.size = 0, label.padding = ggplot2::unit(0.03, "lines")) +
  geom_segment(aes(x = lab_x, xend = lab_x, y = label_y, yend = top_y),
               linewidth = LW_AXIS, colour = SEM[["neutral"]]) +
  geom_text(aes(x = lab_x, y = label_y, label = feature_label),
            size = PT(SZ_TICK), family = FIG_FONT, colour = "grey15", vjust = 0) +
  scale_fill_manual(values = c(`TRUE` = SEM[["case"]], `FALSE` = SEM[["control"]]),
                    guide = "none") +
  scale_x_continuous("Prediction", limits = c(0.30, nrow(wf) + 1.90),
                     breaks = seq_len(nrow(wf)), labels = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(ylim) - pad * 1.35, label_y + pad * 0.45),
                     breaks = NULL) +
  annotate("text", x = nrow(wf) + 0.45, y = min(ylim) - pad * 1.1, hjust = 1, vjust = 1,
           label = sprintf("f(x)=%.3g", fx),
           size = PT(SZ_TICK), family = FIG_FONT, fontface = "bold") +
  annotate("text", x = 0.55, y = E - pad * 0.10, hjust = 0, vjust = 1,
           label = sprintf("E[f(x)]=%.4g", E),
           size = PT(SZ_TICK), family = FIG_FONT, colour = "black") +
  theme_fig3() +
  labs(y = NULL) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = ggplot2::margin(SZ_TICK * 0.5, SZ_TICK * 0.5,
                                      SZ_TICK * 0.4, SZ_TICK * 0.9))
pG <- tag_plot(pG, "G")
pG <- pG + ggplot2::theme(plot.tag.position = c(0, 0.96))
DS_G <- "GBM g1-19 (H2O AutoML) \u00b7 single test cell \u00b7 n = 9,750 test cells \u00b7 SHAP waterfall"
pG <- dataset_line(pG, DS_G)
if (requireNamespace("systemfonts", quietly = TRUE)) {
  w_cap <- systemfonts::string_width(DS_G, family = FIG_FONT, size = SZ_DATA) * 0.3528
  w_lim <- 170 * 0.85
  invisible(NULL)
  stopifnot(w_cap <= w_lim)
}
save_panel(pG, "Figure3G_SHAP_waterfall_test", 170, 42,
           dataset = DS_G)
stopifnot(abs(E + sum(wf$shap_value) - fx) < 5e-3)
