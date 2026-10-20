# =============================================================================
# CMI Explorer: Analysis, Help and About.
# =============================================================================
# The interface process stays small (shiny, bslib, DT); every analysis runs in a separate R
# worker (engine/worker.R) that loads only the model parameters it needs. The workspaces are in
# R/workspace.R, the Help and About pages in R/pages.R.

library(shiny)
library(bslib)
library(DT)

APP_ROOT <- normalizePath(".", mustWork = TRUE)
# A UTF-8 locale before anything is parsed: labels contain Delta, mu, minus and dash characters.
if (!isTRUE(l10n_info()[["UTF-8"]])) for (loc in c("en_US.UTF-8", "C.UTF-8", "C.utf8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
for (f in c("constants.R", "io.R", "bulk_engine.R", "spatial_engine.R")) source(file.path(APP_ROOT, "engine", f), local = TRUE)
for (f in c("figures.R", "jobs.R", "precheck.R", "workspace.R", "pages.R")) source(file.path(APP_ROOT, "R", f), local = TRUE)

if (!isTRUE(ensure_utf8())) warning("CMI Explorer: no UTF-8 locale is available; some characters may not display correctly.")

LIM <- limits_for()
options(shiny.maxRequestSize = LIM$max_upload_mb * 1024^2, shiny.sanitize.errors = TRUE)
addResourcePath("examples", file.path(APP_ROOT, "examples"))   # example inputs and their precomputed result bundles
try(sweep_jobs(idle_min = 30, max_h = 6), silent = TRUE)

# ---- Page ----------------------------------------------------------------------------------------
theme <- bs_theme(version = 5, primary = "#1C2340", secondary = "#4B536E", "body-bg" = "#F5F6F9", "body-color" = "#1C2340",
                  base_font = font_collection("IBM Plex Sans", "system-ui", "-apple-system", "Segoe UI", "Roboto", "Arial", "sans-serif"),
                  "border-radius" = "6px", "link-color" = "#24418F", "focus-ring-color" = "rgba(28, 35, 64, .6)")
brand <- tags$a(class = "brand-link", href = "#analysis", `data-go` = "analysis",
  `aria-label` = "Peking University Third Hospital, CMI Explorer, analysis",
  tags$img(class = "brand-logo", src = "pku3h-logo.jpg", alt = "Peking University Third Hospital",
           width = "44", height = "44", loading = "eager", decoding = "async"),
  tags$span(class = "brand-divider", `aria-hidden` = "true"),
  tags$span(class = "brand-name", "CMI Explorer"))

ui <- page_navbar(
  title = brand, id = "nav", selected = "analysis", theme = theme, lang = "en", fillable = FALSE,
  window_title = "CMI Explorer",
  header = tags$head(tags$link(rel = "stylesheet", href = "styles.css?v=20261021"), tags$script(src = "app.js"),
                     tags$link(rel = "icon", type = "image/svg+xml", href = "favicon.svg"),
                     tags$meta(name = "description", content = "Score single-cell, spatial and bulk B-cell lymphoma data with the clonal malignancy index (CMI).")),
  nav_panel("Analysis", value = "analysis", analysis_page(LIM)),
  nav_panel("Help", value = "help", help_page(LIM)),
  nav_panel("About", value = "about", about_page()),
  footer = tags$footer(class = "app-footer", div(class = "inner",
    span(sprintf("CMI Explorer %s \u00B7 research use only", CMI_VERSION)),
    span("Uploaded files are deleted when you close the page."))))

server <- function(input, output, session) {
  session$allowReconnect(TRUE)
  # Remove job folders of closed or idle sessions now and every 30 minutes while anyone is connected
  # (an instance without sessions is stopped by the host, which clears its temporary storage).
  observe({ invalidateLater(30 * 60 * 1000, session); try(sweep_jobs(idle_min = 30, max_h = 6), silent = TRUE) })
  ws <- list(workspace_server("sc", "single_cell", APP_ROOT, LIM),
             workspace_server("sp", "spatial", APP_ROOT, LIM),
             workspace_server("bulk", "bulk", APP_ROOT, LIM))
  # The page sends a keepalive every few minutes; it marks this session's job folders as in use.
  observeEvent(input$keepalive, for (w in ws) touch_alive(w$dirs()), ignoreInit = TRUE)
}

shinyApp(ui, server)
