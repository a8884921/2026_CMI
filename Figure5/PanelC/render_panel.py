from __future__ import annotations
from publication_data import read_csv, read_parquet, numeric_id, numeric_ids, configured_id, validate_cohorts
import os
import argparse, datetime, hashlib, json, os, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import style2 as S
import fig5_v2 as F
import data2 as D
MM = 1 / 25.4
PROBE_W, PROBE_H, PROBE_OFF = (240.0, 310.0, (20.0, 20.0))
SPECS = {'A': ('Subject-level malignant CMI state composition (42 Visium sections)', 'Visium | 42 sections from 9 public datasets', lambda pg: F.panel_A(pg, 12.5, 14.0, 103.0, S.H_PLOT)), 'B': ('Self-enrichment z of three malignant states', 'Visium | all entities pooled', lambda pg: F.panel_B(pg, 136.5, 14.0, 42.5, S.H_PLOT)), 'C': ('Mal_3 largest connected territory share', 'Visium | QC-retained cohort', lambda pg: F.panel_C(pg, 14.5, 76.0, 30.0, S.H_PLOT)), 'D': ('Visium HD Mal_3 transformation boundary', 'Visium HD', lambda pg: F.panel_D(pg, 50.0, 76.0, 23.0)), 'E': ('GeoMx paired tumor burden, CMI severity and effects', 'GeoMx | Cohort_B | 9 paired / 38 patients', lambda pg: F.panel_E(pg, 96.0, 76.0, 11.0, S.H_PLOT, 14.0, 137.0, 42.0)), 'F': ('Mal_3 territory and CD8 T spatial relationship', 'Visium', lambda pg: F.panel_F(pg, 0.5, 129.5, 38.0, 38.0, 3.0)), 'G': ('Paired CD8 T contact z by malignant state', 'Visium | QC-retained cohort', lambda pg: F.panel_G(pg, 97.5, 136.0, 31.0, S.H_PLOT)), 'H': ('CD8 T abundance (%) by malignant state', 'Visium | QC-retained cohort', lambda pg: F.panel_H(pg, 147.5, 136.0, 31.5, S.H_PLOT)), 'I': ('Discrimination: position and number', 'Visium | QC-retained cohort | 23 patients', lambda pg: F.panel_I(pg, 12.5, 187.0, S.H_PLOT, 49.0)), 'J': ('CosMx immune neighborhoods by CMI decile', 'CosMx | Cohort_P | 77 patients (curves 76)', lambda pg: F.panel_J(pg, 110.0, 187.0, 44.0, S.H_PLOT))}

def ink_bbox(draw):
    pg = S.Page(PROBE_W, PROBE_H, PROBE_OFF)
    draw(pg)
    pg.fig.canvas.draw()
    bb = pg.fig.get_tightbbox(pg.fig.canvas.get_renderer())
    plt.close(pg.fig)
    x0, x1 = (bb.x0 * 25.4 - PROBE_OFF[0], bb.x1 * 25.4 - PROBE_OFF[0])
    y_top = PROBE_H - bb.y1 * 25.4 - PROBE_OFF[1]
    y_bot = PROBE_H - bb.y0 * 25.4 - PROBE_OFF[1]
    return (x0, y_top, x1, y_bot)

def wrap_tag(pg, x, y, text, maxw, size):
    parts = text.split(' | ')
    lines, cur = ([], '')
    for p in parts:
        cand = p if not cur else cur + ' | ' + p
        if S.text_w(cand, size) <= maxw or not cur:
            cur = cand
        else:
            lines.append(cur)
            cur = p
    lines.append(cur)
    for i, ln in enumerate(lines):
        pg.tag(x, y + i * size * 0.3528 * 1.35, ln)
    return len(lines)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', required=True, choices=sorted(SPECS))
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    letter = a.panel
    title, tag, draw = SPECS[letter]
    if letter == 'D':
        tag = f"Visium HD | {configured_id('FIG5_HD_SAMPLE_ID')} | 8 um"
    elif letter == 'F':
        tag = f"Visium | {configured_id('FIG5_SAMPLE_ID')}"
    outdir = Path(a.out) if a.out else Path(os.environ.get('SCI_OUTPUT_DIR', str(HERE / 'results')))
    outdir.mkdir(parents=True, exist_ok=True)
    x0, y_top, x1, y_bot = ink_bbox(draw)
    cw, ch = (max(x1 - x0, 20.0), max(y_bot - y_top, 12.0))
    pad, top = (2.5, 4.0)
    W = cw + 2 * pad
    nlines = max(1, len(tag.split(' | ')) if S.text_w(tag, S.SZ_SMALL) > W - 2 * pad else 1)
    top = 4.0 + nlines * S.SZ_SMALL * 0.3528 * 1.35 + 1.5
    H = ch + top + pad
    pg = S.Page(W, H, (pad - x0, top - y_top))
    F.LOG.clear() if hasattr(F, 'LOG') else None
    draw(pg)
    pg.letter(pad, 1.0, letter)
    wrap_tag(pg, pad, 4.6, tag, W - 2 * pad, S.SZ_SMALL)
    base = outdir / f'Figure5_Panel{letter}'
    pg.fig.savefig(str(base) + '.pdf')
    pg.fig.savefig(str(base) + '.svg')
    pg.fig.savefig(str(base) + '.png', dpi=600)
    plt.close(pg.fig)
    import pandas as pd
    log = getattr(F, 'LOG', [])
    if log:
        pd.DataFrame(log).to_csv(outdir / f'Panel{letter}_values.csv', index=False, encoding='utf-8-sig')
    rec = {'panel': letter, 'canvas_mm': [round(W, 2), round(H, 2)], 'font_size_pt': {'ticks': S.SZ_TICK, 'small_text': S.SZ_SMALL}, 'n_values': len(log)}
    (outdir / f'Figure5_Panel{letter}_validation.json').write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[OK] Panel {letter} -> {os.path.basename(str(base))}.pdf  ({W:.1f} x {H:.1f} mm)')
if __name__ == '__main__':
    main()
