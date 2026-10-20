from __future__ import annotations
import os
import argparse, datetime, json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import figure_s4 as M
MM = 1 / 25.4
SPECS = {'A': ('KM curves, risk tables, HR and log-rank P in 11 cohorts', M.page1, (0, 0, 180, 152)), 'B': ('Time-dependent AUC at 1, 3 and 5 years in 11 cohorts', M.page1, (0, 152, 180, 240)), 'C': ('Integrated AUC over 12-60 months in 9 cohorts', M.page2, (0, 0, 180, 70)), 'D': ('Time-dependent AUC curves in 9 cohorts', M.page2, (0, 70.5, 180, 148)), 'E': ('Multivariable Cox models and 15 sensitivity analyses', M.page2, (0, 148, 180, 242)), 'F': ('Distributions of 18 clinical and molecular features', M.page3, (0, 0, 180, 225))}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', required=True, choices=sorted(SPECS))
    a = ap.parse_args()
    L = a.panel
    title, fn, (x0, y0, x1, y1) = SPECS[L]
    M.setup_fonts()
    D = M.load()
    pad = 2.0
    W, H = (x1 - x0 + 2 * pad, y1 - y0 + 2 * pad)
    fig = plt.figure(figsize=(W * MM, H * MM))
    fn(fig, W, H, pad - x0, pad - y0, D)
    out = Path(os.environ.get('SCI_OUTPUT_DIR', str(HERE / 'results')))
    out.mkdir(parents=True, exist_ok=True)
    base = out / f'FigureS4_Panel{L}'
    fig.savefig(str(base) + '.pdf')
    fig.savefig(str(base) + '.svg')
    fig.savefig(str(base) + '.png', dpi=600)
    plt.close(fig)
    rec = {'panel': L, 'canvas_mm': [round(W, 2), round(H, 2)], 'font_size_pt': M.SZ}
    (out / f'FigureS4_Panel{L}_validation.json').write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[OK] FigureS4 Panel {L} -> {os.path.basename(str(base))}.pdf ({W:.1f} x {H:.1f} mm)')
if __name__ == '__main__':
    main()
