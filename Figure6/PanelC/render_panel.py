from __future__ import annotations
import os
import argparse, datetime, json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import figure6 as M
MM = 1 / 25.4
PW, PH, POX, POY = (260.0, 330.0, 20.0, 20.0)
ROWS_Y, ROWS_H = (97.2, 13 * M.LH * 1.05)
YF = ROWS_Y + ROWS_H + 9.0
YG = YF + 41.5
SPECS = {'A': ('Study design', lambda cv, D: M.panel_A(cv, D, 0, 0, 47.0, 41.0)), 'B': ('C-index of 101 algorithms in 13 cohorts', lambda cv, D: M.panel_B(cv, D, 0, 0, 117.0, 15.4)), 'C': ('Training and test KM curves and time-dependent ROC', lambda cv, D: M.panel_C(cv, D, 0, 0, 15.0, 20.0)), 'D': ('CMI-high versus CMI-low HR and log-rank P in 13 cohorts', lambda cv, D: M.panel_D(cv, D, 0, 0, 8.2, ROWS_H)), 'E': ('C-index comparison with published models in 9 cohorts', lambda cv, D: M.panel_E(cv, D, 0, 0, 8.2, ROWS_H, 110.0)), 'F': ('Risk scores by clinical and molecular features', lambda cv, D: M.panel_F(cv, D, 0, 0, 40.0, 15.0, 23.4, 4.6)), 'G': ('Multivariable Cox models in four disease groups', lambda cv, D: M.panel_G(cv, D, 0, 0, M.H - YG, (M.W - 4.0) / 2))}

def bbox(draw, D):
    fig = plt.figure(figsize=(PW * MM, PH * MM))
    cv = M.Canvas(fig, PW, PH, POX, POY)
    cv.tables = []
    cv.n_points = {}
    draw(cv, D)
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    plt.close(fig)
    return (bb.x0 * 25.4 - POX, PH - bb.y1 * 25.4 - POY, bb.x1 * 25.4 - POX, PH - bb.y0 * 25.4 - POY)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', required=True, choices=sorted(SPECS))
    a = ap.parse_args()
    L = a.panel
    title, draw = SPECS[L]
    M.setup_fonts()
    D = M.load()
    x0, yt, x1, yb = bbox(draw, D)
    cw, ch = (max(x1 - x0, 18.0), max(yb - yt, 10.0))
    pad = 2.0
    W, H = (cw + 2 * pad, ch + 2 * pad)
    fig = plt.figure(figsize=(W * MM, H * MM))
    cv = M.Canvas(fig, W, H, pad - x0, pad - yt)
    cv.tables = []
    cv.n_points = {}
    draw(cv, D)
    out = Path(os.environ.get('SCI_OUTPUT_DIR', str(HERE / 'results')))
    out.mkdir(parents=True, exist_ok=True)
    base = out / f'Figure6_Panel{L}'
    fig.savefig(str(base) + '.pdf')
    fig.savefig(str(base) + '.svg')
    fig.savefig(str(base) + '.png', dpi=600)
    plt.close(fig)
    rec = {'panel': L, 'canvas_mm': [round(W, 2), round(H, 2)], 'font_size_pt': M.SZ}
    (out / f'Figure6_Panel{L}_validation.json').write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[OK] Figure6 Panel {L} -> {os.path.basename(str(base))}.pdf ({W:.1f} x {H:.1f} mm)')
if __name__ == '__main__':
    main()
