set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
"$RSCRIPT_BIN" 04_panelD.R
"$RSCRIPT_BIN" 05_sensitivity_kruskal_wallis.R
"$RSCRIPT_BIN" 06_sensitivity_ordinal_trend.R
