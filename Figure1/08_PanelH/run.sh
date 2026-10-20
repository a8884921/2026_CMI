set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
"$RSCRIPT_BIN" 08_panelH.R
