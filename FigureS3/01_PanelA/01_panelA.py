from publication_data import read_csv, read_parquet, numeric_id, numeric_ids, configured_id, validate_cohorts
import json
import os
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
STEP = HERE
PROJECT_ROOT = os.path.dirname(os.path.dirname(HERE))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
input_dir = os.environ.get('SCI_INPUT_DIR', os.path.join(HERE, 'inputs'))
input_dir = os.path.normpath(input_dir)
plot_data_dir = os.path.join(input_dir, 'plot_data')
results_dir = os.environ.get('SCI_OUTPUT_DIR', os.path.join(HERE, 'results'))
os.makedirs(results_dir, exist_ok=True)
figure_name = 'panelA_spatial_atlas'
PT = 1.0
canvas_width = 595.276
page_left_margin = 20.0
page_right_margin = 20.0
grid_height = 59.5
grid_gap = 12.0
column_gap = 8.0
disease_header_height = 22.0
card_gap = 16.0
legend_row_height = 8.0
legend_start_y = 52.0
legend_row_y = [52.0, 60.0, 68.3]
font_sizes = dict(main_title=24.0, panel_letter=18.0, disease_title=7.3, disease_title_small=5.8, legend=5.3, sample_title=3.4, subtitle=6.0)
disease_order = ['DLBCL', 'PCNSL', 'MZL', 'BL']
disease_columns = {'DLBCL': 6, 'PCNSL': 4, 'MZL': 3, 'BL': 1}
cell_types = ['CD4T', 'CD8T', 'DC', 'Endothelial', 'Macrophage', 'Mast', 'NK', 'Normal_B', 'Stromal', 'Mal_1', 'Mal_2', 'Mal_3']
malignant_types = {'Mal_1', 'Mal_2', 'Mal_3'}
colors = {'Mal_1': '#0072B2', 'Mal_2': '#D55E00', 'Mal_3': '#009E73', 'CD4T': '#F5B7B1', 'CD8T': '#B3CDE3', 'NK': '#DECBE4', 'Macrophage': '#FED9A6', 'DC': '#E5D8BD', 'Mast': '#FADBD8', 'Normal_B': '#A9DFBF', 'Stromal': '#D5D8DC', 'Endothelial': '#AED6F1'}
fp_spot = os.path.join(plot_data_dir, 'panelA_r01_spatial_atlas.csv')
fp_smp = os.path.join(plot_data_dir, 'panelA_r02_sample_statistics.csv')
fp_ct = os.path.join(plot_data_dir, 'panelA_r03_cell_type_counts.csv')
fp_ck = os.path.join(plot_data_dir, 'panelA_r05_checks.json')
for fp in (fp_spot, fp_smp, fp_ct, fp_ck):
    if not os.path.exists(fp):
        raise SystemExit('Missing plot data. Run the extraction script first.')
spot = read_csv(fp_spot)
samples = read_csv(fp_smp)
counts = read_csv(fp_ct)
checks = json.load(open(fp_ck, encoding='utf-8'))
assert samples['sample_id'].notna().all() and samples['sample_id'].is_unique
assert spot['sample_id'].notna().all() and set(spot['sample_id']) == set(samples['sample_id'])
assert len(samples) == 42, f'Expected 42 samples; found {len(samples)}'
assert int(checks['numbers']['total_spots']) == len(spot), 'Spot count differs from the validation data.'
class_counts = {r.cell_type: int(r.N) for r in counts.itertuples()}
assert class_counts == checks['numbers']['class_counts'], 'Cell type counts differ from the validation data.'
total_spots = len(spot)
n_malignant = int(sum((class_counts.get(c, 0) for c in malignant_types)))
n_unassigned = class_counts.get('Unassigned', 0)
print(f'42 samples; {total_spots:,} spots; top-1 malignant {n_malignant:,} ({100 * n_malignant / total_spots:.1f}%) | Unassigned {n_unassigned}')

def card_inner_width(n_columns):
    return n_columns * grid_height + (n_columns - 1) * column_gap

def card_rows(disease):
    return int(np.ceil(int((samples.disease == disease).sum()) / disease_columns[disease]))
w_DLBCL = card_inner_width(6)
w_PCNSL = card_inner_width(4)
w_MZL = card_inner_width(3)
w_BL = card_inner_width(1)
x_DLBCL = (page_left_margin + (canvas_width - page_right_margin)) / 2 - w_DLBCL / 2
x_PCNSL = page_left_margin
x_MZL = canvas_width - page_right_margin - w_MZL
x_BL = x_DLBCL
if x_MZL + w_MZL > canvas_width - page_right_margin - 8:
    x_MZL = canvas_width - page_right_margin - w_MZL - 8
if x_MZL < x_PCNSL + w_PCNSL + 30:
    x_MZL = (page_left_margin + canvas_width - page_right_margin) / 2 - w_MZL / 2
y_DLBCL_header = 85.0
y_DLBCL_grid = y_DLBCL_header + disease_header_height
y_second_header = y_DLBCL_grid + card_rows('DLBCL') * (grid_height + grid_gap) - grid_gap + card_gap
y_second_grid = y_second_header + disease_header_height
y_BL_header = y_second_grid + card_rows('PCNSL') * (grid_height + grid_gap) - grid_gap + card_gap
y_BL_grid = y_BL_header + disease_header_height
canvas_height = y_BL_grid + card_rows('BL') * (grid_height + grid_gap) - grid_gap + 12.0
print(f'Canvas: {canvas_width:.1f} x {canvas_height:.1f} pt = {canvas_width / 72 * 25.4:.1f} x {canvas_height / 72 * 25.4:.1f} mm')
fig = plt.figure(figsize=(canvas_width / 72, canvas_height / 72))
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'Arial'
draw_text = lambda x, y, s, size, **kw: fig.text(x / canvas_width, 1 - y / canvas_height, s, fontsize=size, family='Arial', ha=kw.pop('ha', 'left'), va=kw.pop('va', 'top'), **kw)
title_top = 26.0
fig.text(11.0 / canvas_width, 1 - title_top / canvas_height, 'Supplementary Figure 3', fontsize=font_sizes['main_title'], family='Arial', fontweight='bold', ha='left', va='top')
fig.text(20.5 / canvas_width, 1 - 52.0 / canvas_height, 'A', fontsize=font_sizes['panel_letter'], family='Arial', ha='left', va='top')
legend_columns = {}
_row = {'Mal_1': 0, 'Mal_2': 1, 'Mal_3': 2, 'CD4T': 0, 'CD8T': 1, 'DC': 2, 'Endothelial': 0, 'Macrophage': 1, 'Mast': 2, 'NK': 0, 'Normal_B': 1, 'Stromal': 2}
_column_x = {'Mal_1': 81.46, 'Mal_2': 81.46, 'Mal_3': 81.46, 'CD4T': 175.59, 'CD8T': 175.59, 'DC': 175.59, 'Endothelial': 268.8, 'Macrophage': 268.8, 'Mast': 268.8, 'NK': 380.23, 'Normal_B': 380.23, 'Stromal': 380.23}
point_diameter = 5.45
for ct in cell_types:
    dx = _column_x[ct]
    dy = legend_row_y[_row[ct]]
    fig.patches.append(plt.Circle((dx / canvas_width, 1 - (dy + point_diameter / 2) / canvas_height), point_diameter / 2 / canvas_width, transform=fig.transFigure, facecolor=colors[ct], edgecolor='none', zorder=5))
    draw_text(dx + point_diameter + 4.9, dy - 0.3, f'{ct}  (N = {class_counts.get(ct, 0):,} | {100 * class_counts.get(ct, 0) / total_spots:.1f}%)', font_sizes['legend'], va='top')
if n_unassigned:
    draw_text(_column_x['CD4T'], legend_row_y[2] + 8.3, f'Unassigned  (N = {n_unassigned:,} | {100 * n_unassigned / total_spots:.2f}%)', font_sizes['legend'], va='top', color='#666666')
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

def text_width_pt(s, size):
    if not s:
        return 0.0
    p = TextPath((0, 0), s, size=size, prop=FontProperties(family='Arial'))
    return p.get_extents().width

def draw_card(disease, x0, y_header, y_grid):
    subset = samples[samples.disease == disease].sort_values('panel_order')
    n_columns = disease_columns[disease]
    n_rows = card_rows(disease)
    inner_width = card_inner_width(n_columns)
    height = disease_header_height + n_rows * grid_height + (n_rows - 1) * grid_gap
    fig.add_artist(Rectangle((x0 / canvas_width, 1 - (y_header + height) / canvas_height), inner_width / canvas_width, height / canvas_height, facecolor='white', edgecolor='#d8d8d8', lw=0.7, zorder=0))
    n_spot = int(subset.n_spot.sum())
    n_mal = int(subset.n_malignant.sum())
    s = f'{disease}   |   n = {len(subset)} patients   |   N = {n_spot:,} spots   |   top-1 malignant {n_mal:,} ({100 * n_mal / n_spot:.1f}%)'
    available = min(inner_width, 2 * (canvas_width / 2 - max(page_left_margin, canvas_width - page_right_margin - inner_width - x0)))
    sz = font_sizes['disease_title']
    if text_width_pt(s, sz) > available:
        sz = max(font_sizes['disease_title_small'], sz * available / text_width_pt(s, sz))
    draw_text(x0 + inner_width / 2, y_header + disease_header_height * 0.55, s, sz, ha='center', va='center', color='#111111')
    for i, r in enumerate(subset.itertuples()):
        row, column_i = divmod(i, n_columns)
        px = x0 + column_i * (grid_height + column_gap)
        py = y_grid + row * (grid_height + grid_gap)
        ax = fig.add_axes([px / canvas_width, 1 - (py + grid_height) / canvas_height, grid_height / canvas_width, grid_height / canvas_height])
        g = spot[spot.sample_id == r.sample_id]
        for ct in ['Stromal', 'Endothelial', 'Mast', 'DC', 'Normal_B', 'NK', 'CD4T', 'CD8T', 'Macrophage'] + ['Mal_1', 'Mal_2', 'Mal_3']:
            m = g[g.top1_celltype == ct]
            if len(m):
                ax.scatter(m.x, m.y, s=1.6, c=colors[ct], linewidths=0, rasterized=True, zorder=2 if ct not in malignant_types else 3)
        x0s, x1s = (float(g.x.min()), float(g.x.max()))
        y0s, y1s = (float(g.y.min()), float(g.y.max()))
        pad = max(x1s - x0s, y1s - y0s) * 0.03
        ax.set_xlim(x0s - pad, x1s + pad)
        ax.set_ylim(y1s + pad, y0s - pad)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_facecolor('none')
        draw_text(px + grid_height / 2, py - 2.0, f'{r.sample_id} | {r.n_spot:,} spots | Mal {round(100 * r.mal_frac):.0f}%', font_sizes['sample_title'], ha='center', va='bottom', color='#111111')
draw_card('DLBCL', x_DLBCL, y_DLBCL_header, y_DLBCL_grid)
draw_card('PCNSL', x_PCNSL, y_second_header, y_second_grid)
draw_card('MZL', x_MZL, y_second_header, y_second_grid)
draw_card('BL', x_BL, y_BL_header, y_BL_grid)
dataset_label = f'Datasets: Visium spatial transcriptomics, 42 mature B-cell lymphoma slides (DLBCL n=30 | PCNSL n=8 | MZL n=3 | BL n=1) | 1 patient = 1 slide | n = 42 slides | N = {total_spots:,} spots | top-1 malignant {n_malignant:,} ({100 * n_malignant / total_spots:.1f}%)'
draw_text(canvas_width / 2, canvas_height - 6.0, dataset_label, font_sizes['sample_title'] * 1.55, ha='center', va='bottom', color='#8A9197')
base = os.path.join(results_dir, figure_name)
for extension in ('pdf', 'svg', 'png'):
    fig.savefig(base + '.' + extension, dpi=int(os.environ.get('SCI_DPI', '600')))
plt.close(fig)
print('Spatial atlas saved.')
