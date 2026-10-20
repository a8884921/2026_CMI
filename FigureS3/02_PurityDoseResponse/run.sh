set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
PYTHON_BIN="${PYTHON_BIN:-python3}"
"$PYTHON_BIN" 00_extract_plot_data.py
"$PYTHON_BIN" 02_purity_dose_response.py
