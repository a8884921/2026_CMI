from publication_data import read_csv, read_parquet, numeric_id, numeric_ids, configured_id, validate_cohorts
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PolyCollection, LineCollection
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('SCI_INPUT_DIR', os.path.join(HERE, 'inputs'))
sys.path.insert(0, HERE)
import style2 as S
import data2 as D
S.apply()
C = S.C
LOG = []
ST = read_csv(os.path.join(DATA, 'stats', 'stats_supplement.csv')).set_index('key')

def stat(key):
    return float(ST.loc[key, 'value'])

def note(panel, item, value, source):
    LOG.append(dict(panel=panel, item=item, value=value))

def blend(ax):
    return mtrans.blended_transform_factory(ax.transData, ax.transAxes)

def box(ax, x, v, col, w=0.52, fill=None, lw_box=S.LW_ERR, lw_med=S.LW_DATA * 1.5, z=3):
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    iqr = q3 - q1
    lo = v[v >= q1 - 1.5 * iqr].min()
    hi = v[v <= q3 + 1.5 * iqr].max()
    ax.add_patch(Rectangle((x - w / 2, q1), w, q3 - q1, fc=fill or 'none', ec=col, lw=lw_box, zorder=z))
    ax.plot([x - w / 2, x + w / 2], [med, med], color=col, lw=lw_med, zorder=z + 1, solid_capstyle='butt')
    ax.plot([x, x], [q3, hi], color=col, lw=lw_box, zorder=z)
    ax.plot([x, x], [lo, q1], color=col, lw=lw_box, zorder=z)
    return med

def jitter_for(keys, width=0.16, seed=11):
    rng = np.random.default_rng(seed)
    return {k: rng.uniform(-width, width) for k in sorted(keys)}

def panel_A(pg, x0, y0, w, h):
    d = D.composition().copy()
    gap = {'PCNSL': 1.3, 'MZL': 1.9, 'BL': 2.8}
    pos, cur, prev = ([], 0.0, None)
    for dz in d.disease:
        if prev is not None and dz != prev:
            cur += gap[dz]
        pos.append(cur)
        cur += 1.0
        prev = dz
    d['xpos'] = pos
    ax = pg.ax(x0, y0, w, h)
    bw = 0.82
    b3, b2, b1 = (d.frac_mal3 * 100, d.frac_mal2 * 100, d.frac_mal1 * 100)
    ax.bar(d.xpos, b3, bw, color=C['Mal_3'], lw=0)
    ax.bar(d.xpos, b2, bw, bottom=b3, color=C['Mal_2'], lw=0)
    ax.bar(d.xpos, b1, bw, bottom=b3 + b2, color=C['Mal_1'], lw=0)
    ax.set_xlim(-0.7, d.xpos.max() + 0.7)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_xticks([])
    ax.set_ylabel('Malignant signal (%)')
    tr = blend(ax)
    centres = {}
    for dz in S.DISEASES:
        q = d[d.disease == dz]
        a, b = (q.xpos.min(), q.xpos.max())
        centres[dz] = (a + b) / 2
        ax.plot([a - 0.42, b + 0.42], [-0.035, -0.035], color=C['neutral'], lw=S.LW_AXIS, transform=tr, clip_on=False)
        ax.text((a + b) / 2, -0.075, dz, ha='center', va='top', fontsize=S.SZ_TICK, transform=tr, clip_on=False)
        ax.text((a + b) / 2, -0.075 - S.SZ_TICK * 1.18 * S.PT / h, f'n = {len(q)}', ha='center', va='top', fontsize=S.SZ_SMALL, color=C['dark'], transform=tr, clip_on=False)
    p = stat('landscape_kw')
    S.omnibus(ax, d.xpos.min() - 0.42, d.xpos.max() + 0.42, 1.05, S.Ptxt(p))
    note('A', 'Kruskal-Wallis P, Mal_3 share across DLBCL/PCNSL/MZL/BL', S.fmtP(p), 'a02 landscape_kw (panelA_values frac_mal3)')
    for dz in S.DISEASES:
        note('A', f'Mal_3-dominant sections {dz}', ST.loc[f'landscape_dominant_{dz}', 'value'], 'a02 (panelA_values)')
    return ax

def panel_B(pg, x0, y0, w, h):
    t = D.self_z()
    pts = t[(t['row_type'] == 'point') & (t.facet != 'ALL (pooled)')].copy()
    pts['label'] = numeric_ids(pts['label']).array
    W = pts.pivot_table(index='label', columns='x_group', values='z')[S.STATES]
    dz = pts.drop_duplicates('label').set_index('label').facet
    ax = pg.ax(x0, y0, w, h)
    ax.set_yscale('asinh', linear_width=5)
    jit = jitter_for(W.index)
    ns = []
    for j, s in enumerate(S.STATES):
        v = W[s].dropna()
        ns.append(len(v))
        med = box(ax, j, v.values, C[s], fill=C['Mal_3_light'] if s == 'Mal_3' else None, lw_box=S.LW_DATA if s == 'Mal_3' else S.LW_ERR)
        for sid, val in v.items():
            S.scatter_marker(ax, j + jit[sid], val, dz[sid], C[s], 11, lw=0.25, edgecolor='white', zorder=5)
        note('B', f'median self-enrichment z {s}', f'{med:.2f} (n = {len(v)})', 'c03_r07 point rows (frozen)')
    S.zero_line(ax)
    ax.set_ylim(-6, 230)
    ax.set_yticks([0, 10, 30, 100])
    ax.set_yticklabels(['0', '10', '30', '100'])
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_xlim(-0.55, 2.55)
    S.state_ticks(ax, ns)
    ax.set_ylabel('Self-enrichment (z)')
    ps = [stat(f'self_vs0_{s}') for s in S.STATES]
    S.group_p(ax, ps, xs=[0, 1, 2.45], levels=[0, 1, 0], ha=['center', 'center', 'right'])
    for s, p, n in zip(S.STATES, ps, ns):
        note('B', f'one-sample Wilcoxon vs 0, {s} (n = {n})', S.fmtP(p), 'a02 self_vs0_* (c03_r07 point rows)')
    note('B', 'Friedman over the 23 complete sections (legend only)', S.fmtP(stat('self_friedman')), 'a02 self_friedman = frozen')
    return ax

def panel_C(pg, x0, y0, w, h):
    t = D.territory_sizes_all()
    t = t[(t.disease == 'DLBCL') & (t.n_state >= 20)]
    W = t.pivot(index='sample_id', columns='state', values='share').dropna()[S.STATES]
    ax = pg.ax(x0, y0, w, h)
    jit = jitter_for(W.index, seed=7)
    for sid, row in W.iterrows():
        ax.plot(np.arange(3) + jit[sid], row.values * 100, color='#D3D6D8', lw=0.45, zorder=1)
    for j, s in enumerate(S.STATES):
        v = W[s].values * 100
        med = box(ax, j, v, C[s], fill=C['Mal_3_light'] if s == 'Mal_3' else None, lw_box=S.LW_DATA if s == 'Mal_3' else S.LW_ERR)
        ax.scatter(j + np.array([jit[k] for k in W.index]), v, s=11, color=C[s], lw=0.25, edgecolor='white', zorder=5)
        note('C', f'median territory share {s}', f'{med:.1f}% (n = {len(v)})', 'a02 territory_share_per_section.csv')
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(-4, 118)
    ax.set_yticks([0, 50, 100])
    ax.set_yticklabels(['0', '50', '100'])
    S.state_ticks(ax, [len(W)] * 3)
    ax.set_ylabel('Largest territory (% of state)')
    p13 = stat('territory_pair_Mal_3_Mal_1')
    S.bracket(ax, 0, 2, 1.05, 0.03, S.Ptxt(p13))
    note('C', 'paired Wilcoxon Mal_3 vs Mal_1', S.fmtP(p13), 'a02 territory_pair_Mal_3_Mal_1')
    note('C', 'paired Wilcoxon Mal_3 vs Mal_2 (legend)', S.fmtP(stat('territory_pair_Mal_3_Mal_2')), 'a02 territory_pair_Mal_3_Mal_2')
    note('C', 'Friedman (legend)', S.fmtP(stat('territory_friedman')), 'a02 territory_friedman')
    return ax

def panel_D(pg, x0, y0, w, case=None):
    case = configured_id('FIG5_HD_SAMPLE_ID', case)
    maps = read_parquet(os.path.join(DATA, 'assets', 'visiumhd_maps.parquet'))
    meta = json.load(open(os.path.join(DATA, 'assets', 'visiumhd_meta.json'), encoding='utf-8'))
    for map_key in maps['map'].dropna():
        numeric_id(str(map_key).split('_', 1)[0])
    for map_key in meta['maps']:
        numeric_id(str(map_key).split('_', 1)[0])
    vmax = 0.684428
    ax = pg.ax(x0, y0, w, w)
    S.hd_map(ax, case, 'mal3', meta, maps)
    for name, (u, v) in S.zone_anchor(maps, case).items():
        ax.text(u, v, name, ha='center', va='center', fontsize=S.SZ_SMALL, color='white', bbox=dict(boxstyle='round,pad=0.2', fc='#2B2B2B', ec='none', alpha=0.6))
    cax = pg.ax(x0 + w * 0.12, y0 + w + 2.0, w * 0.76, 1.6)
    S.vector_cbar(cax, 'inferno', 0, vmax, [0, vmax], ['0', f'{vmax:.2f}'])
    pg.text(x0 + w / 2, y0 + w + 7.4, 'Mal_3 weight', size=S.SZ_SMALL, ha='center', color=C['dark'])
    bar = 1000.0 / 7200.0 * w
    pg.line([x0 + w - bar, x0 + w], [y0 - 2.0] * 2, lw=S.LW_DATA)
    pg.text(x0 + w - bar - 1.2, y0 - 2.0, '1 mm', size=S.SZ_SMALL, ha='right', va='center')
    note('D', f'{case}: Mal_3 weight colour scale', f'0-{vmax:.2f} (8 um bins; white line = zone boundary)', 'panel M SVG (c03.6 zones, c05_r01 bins)')
    return ax

def panel_E(pg, x0, y0, w1, h, gap, xf, wf):
    pr = D.geomx_pairs()
    ts = D.geomx_tests()
    md = D.geomx_models()
    specs = [('malignant_burden', 'Tumour burden', (0.3, 1.0), [0.4, 0.7, 1.0]), ('CMI_severity', 'CMI severity', (-0.45, 0.65), [-0.4, 0, 0.4])]
    for i, (ind, ylab, yl, yt) in enumerate(specs):
        ax = pg.ax(x0 + i * (w1 + gap), y0, w1, h)
        q = pr[pr.indicator == ind]
        for _, rr in q.iterrows():
            up = rr.value_tFL > rr.value_FL
            c, lw, z = (C['Mal_3'], S.LW_DATA, 3) if ind == 'CMI_severity' and up else (C['neutral'], S.LW_ERR, 2)
            ax.plot([0, 1], [rr.value_FL, rr.value_tFL], color=c, lw=lw, zorder=z)
            ax.scatter([0, 1], [rr.value_FL, rr.value_tFL], s=9, color=c, zorder=z + 1, lw=0)
        t = ts.loc[ind]
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['FL', 'tFL'])
        ax.set_xlim(-0.35, 1.35)
        ax.set_ylim(*yl)
        ax.set_yticks(yt)
        ax.set_yticklabels([S.minus(f'{v:g}') for v in yt])
        ax.set_ylabel(ylab)
        ax.text(0.5, 1.05, S.Ptxt(t.p_exact), transform=ax.transAxes, ha='center', va='bottom', fontsize=S.SZ_SMALL, color=C['Mal_3'] if t.p_exact < 0.05 else C['dark'], clip_on=False)
        note('E', f'{ind}: n up/down, exact paired Wilcoxon P', f'{int(t.n_up)}/{int(t.n_down)}, P = {S.fmtP(t.p_exact)}', 'c04_r02 paired_test_summary')
    f = md[md.term == 'disease_statetFL'].reset_index(drop=True)
    labs = ['Paired (9)', 'Paired (9), burden-adj.', 'Cohort (38), burden-adj.']
    axf = pg.ax(xf, y0, wf, h)
    yy = np.arange(len(f))[::-1].astype(float)
    for yv, (_, rr), l in zip(yy, f.iterrows(), labs):
        axf.plot([rr.lo, rr.hi], [yv, yv], color=C['Mal_3'], lw=S.LW_ERR)
        axf.scatter([rr.estimate], [yv], s=16, marker='s', color=C['Mal_3'], zorder=3, lw=0)
        axf.text(0.018, yv + 0.38, l, va='center', ha='left', fontsize=S.SZ_SMALL, color=C['ink'], clip_on=False)
        note('E', f'forest {l}: beta [95% CI] (P in legend)', f'{rr.estimate:.3f} [{rr.lo:.3f}, {rr.hi:.3f}], P = {S.fmtP(rr.p)}', 'c05_r01 term disease_statetFL')
    axf.axvline(0, color=C['zero'], lw=S.LW_AXIS, ls=(0, (2.5, 2)), zorder=0)
    axf.set_ylim(-0.5, len(f) - 0.2)
    axf.set_yticks([])
    axf.spines['left'].set_visible(False)
    axf.set_xlim(-0.05, 0.55)
    axf.set_xticks([0, 0.25, 0.5])
    axf.set_xticklabels(['0', '0.25', '0.5'])
    axf.set_xlabel('beta (95% CI)')
    return axf

def spot_s(pitch_px, mm_per_px, fill=1.0):
    return (pitch_px * mm_per_px * fill / S.PT) ** 2

def map_axes(pg, x0, y0, w, h, xs, ys, pad):
    ax = pg.ax(x0, y0, w, h)
    xmin, xmax = (xs.min() - pad, xs.max() + pad)
    ymin, ymax = (ys.min() - pad, ys.max() + pad)
    s = max((xmax - xmin) / w, (ymax - ymin) / h)
    cx, cy = ((xmin + xmax) / 2, (ymin + ymax) / 2)
    ax.set_xlim(cx - s * w / 2, cx + s * w / 2)
    ax.set_ylim(cy + s * h / 2, cy - s * h / 2)
    S.blank(ax)
    return (ax, 1.0 / s)

def panel_F(pg, x0, y0, wm, hm, gap, sample=None):
    sample = configured_id('FIG5_SAMPLE_ID', sample)
    sp = D.spots()
    d = sp[sp.sample_id == sample]
    summ = D.territory_summary()
    pitch = float(summ[summ.sample_id == sample].sx_px.iloc[0])
    tx, ty = D.territory_outline('Mal_3', sample)
    lab = d.CMI_state_top1.values
    for i, mode in enumerate(['state', 'cd8']):
        ax, mmpp = map_axes(pg, x0 + i * (wm + gap), y0, wm, hm, d.x.values, d.y.values, pitch * 0.6)
        col = np.full(len(d), C['spot_bg'], dtype=object)
        z = np.zeros(len(d))
        if mode == 'state':
            for st in S.STATES:
                m = lab == st
                col[m] = C[st]
                z[m] = 1
        else:
            m3 = lab == 'Mal_3'
            col[m3] = C['Mal_3_light']
            z[m3] = 1
            c8 = lab == 'CD8T'
            col[c8] = C['CD8T']
            z[c8] = 2
        o = np.argsort(z, kind='stable')
        ax.scatter(d.x.values[o], d.y.values[o], s=spot_s(pitch, mmpp, 1.12), c=list(col[o]), lw=0, marker='o')
        ax.plot(np.r_[tx, tx[:1]], np.r_[ty, ty[:1]], color=C['ink'], lw=S.LW_ERR, ls=(0, (3, 1.8)))
        if i == 0:
            bar_mm = 10 * pitch * mmpp
            pg.line([x0 + 0.5, x0 + 0.5 + bar_mm], [y0 + hm + 1.6] * 2, lw=S.LW_DATA)
            pg.text(x0 + 0.5 + bar_mm + 1.2, y0 + hm + 1.6, '1 mm', size=S.SZ_SMALL, va='center')
    note('F', 'maps', f"{sample}: Mal_3 territory {int(summ[(summ.sample_id == sample) & (summ.state == 'Mal_3')].n_territory.iloc[0]):,}/{int(summ[(summ.sample_id == sample) & (summ.state == 'Mal_3')].n_state_spots.iloc[0]):,} Mal_3 spots", 'c02.4_r03 (frozen)')
    return None

def panel_G(pg, x0, y0, w, h):
    W = read_csv(os.path.join(DATA, 'stats', 'cd8t_contact_per_section.csv'), index_col=0)[S.STATES]
    acc = D.tme_access().loc['CD8T']
    ax = pg.ax(x0, y0, w, h)
    W.index = numeric_ids(W.index).array
    jit = jitter_for(W.index, seed=5)
    for sid, row in W.iterrows():
        ax.plot(np.arange(3) + jit[sid], row.values, color='#D3D6D8', lw=0.45, zorder=1)
    for j, s in enumerate(S.STATES):
        v = W[s].values
        box(ax, j, v, C[s], fill=C['Mal_3_light'] if s == 'Mal_3' else None, lw_box=S.LW_DATA if s == 'Mal_3' else S.LW_ERR)
        ax.scatter(j + np.array([jit[k] for k in W.index]), v, s=11, color=C[s], lw=0.25, edgecolor='white', zorder=5)
    S.zero_line(ax)
    ax.set_yscale('asinh', linear_width=3)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(-95, 16)
    ax.set_yticks([-50, -10, 0, 10])
    ax.set_yticklabels([S.minus('-50'), S.minus('-10'), '0', '10'])
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    S.state_ticks(ax, [len(W)] * 3)
    ax.set_ylabel('CD8 T contact (z)')
    p13 = float(acc.p_Mal1vMal3)
    S.bracket(ax, 0, 2, 1.05, 0.03, S.Ptxt(p13))
    note('G', 'paired Wilcoxon Mal_1 vs Mal_3 (frozen c05_r03)', S.fmtP(p13), 'c05_r03 p_Mal1vMal3')
    note('G', 'Friedman (legend; frozen = recomputed)', S.fmtP(float(acc.friedman_p)), 'c05_r03 friedman_p')
    note('G', 'medians Mal_1/2/3 (frozen)', f'{acc.median_z_Mal_1:.2f}/{acc.median_z_Mal_2:.2f}/{acc.median_z_Mal_3:.2f}', 'c05_r03; points = a02 recomputed per section')
    return ax

def panel_H(pg, x0, y0, w, h):
    W = read_csv(os.path.join(DATA, 'stats', 'cd8t_number_per_section.csv'), index_col=0)[S.STATES]
    qlab = D.tme_composition().dropna(subset=['label']).set_index('y').label['CD8T']
    ax = pg.ax(x0, y0, w, h)
    W.index = numeric_ids(W.index).array
    jit = jitter_for(W.index, seed=9)
    for sid, row in W.iterrows():
        ax.plot(np.arange(3) + jit[sid], row.values, color='#D3D6D8', lw=0.45, zorder=1)
    for j, s_ in enumerate(S.STATES):
        v = W[s_].values
        med = box(ax, j, v, C[s_], fill=C['Mal_3_light'] if s_ == 'Mal_3' else None, lw_box=S.LW_DATA if s_ == 'Mal_3' else S.LW_ERR)
        ax.scatter(j + np.array([jit[k] for k in W.index]), v, s=11, color=C[s_], lw=0.25, edgecolor='white', zorder=5)
        note('H', f'median CD8 T abundance {s_}', f'{med:.2f}% (n = {len(v)})', 'a02 cd8t_number_per_section.csv')
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(-3, 52)
    ax.set_yticks([0, 25, 50])
    ax.set_yticklabels(['0', '25', '50'])
    S.state_ticks(ax, [len(W)] * 3)
    ax.set_ylabel('CD8 T abundance (%)')
    S.omnibus(ax, 0, 2, 1.05, qlab, color=C['dark'])
    note('H', 'Friedman across the three states, BH (frozen c03_r06)', qlab, 'c03_r06 label column')
    note('H', 'Mal_3 higher than Mal_1 (legend)', f'{int((W.Mal_3 > W.Mal_1).sum())}/{len(W)} sections', 'a02 cd8t_number_per_section.csv')
    return ax

def panel_I(pg, x0, y0, side, xt):
    r = D.roc()
    A = r[r['panel'] == 'A'].copy()
    A['x'] = A.x.astype(float)
    perf = D.roc_perf()
    perf = perf[perf['algorithm'] == 'Gradient boosting GBM'].set_index('feature_set')
    sets = [('Position-only (infiltration depth)', 'Position', C['ink'], S.LW_DATA * 1.5, '-'), ('Number + position (combined)', 'Both', C['neutral'], S.LW_DATA, (0, (2.6, 1.6))), ('Number-only (composition)', 'Number', C['light'], S.LW_DATA * 1.3, '-')]
    ax = pg.ax(x0, y0, side, side)
    ax.plot([0, 1], [0, 1], color=C['faint'], lw=S.LW_AXIS, ls=(0, (1, 1.6)))
    for key, lab, c, lw, ls in sets[::-1]:
        q = A[A.group == f'Gradient boosting GBM | {key}'].sort_values(['x', 'y'])
        ax.plot(q.x, q.y, color=c, lw=lw, ls=ls, solid_capstyle='round')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.set_xticklabels(['0', '0.5', '1'])
    ax.set_yticklabels(['0', '0.5', '1'])
    ax.set_xlabel('1 - specificity')
    ax.set_ylabel('Sensitivity')
    yy = y0 + 0.6
    rows_y = {}
    for key, lab, c, lw, ls in sets:
        pr = perf.loc[key]
        pg.line([xt, xt + 4.5], [yy + 1.6, yy + 1.6], color=c, lw=lw, ls=ls)
        pg.text(xt + 5.6, yy, lab, size=S.SZ_TICK, va='top')
        pg.text(xt + 5.6 + 14.6, yy, f'{pr.Mal3vs1_AUC:.3f}', size=S.SZ_TICK, va='top')
        pg.text(xt + 5.6, yy + 3.75, f'95% CI {pr.Mal3vs1_CI_low:.3f}-{pr.Mal3vs1_CI_high:.3f}', size=S.SZ_SMALL, va='top', color=C['dark'])
        rows_y[lab] = yy + 1.6
        note('I', f'AUC Mal_3 vs Mal_1, {lab}', f'{pr.Mal3vs1_AUC:.3f} ({pr.Mal3vs1_CI_low:.3f}-{pr.Mal3vs1_CI_high:.3f})', 'c07.3_r01 performance (GBM)')
        yy += 8.3
    dd = r[r['panel'] == 'D'].set_index('x')
    yy += 1.4
    pg.line([xt, xt + 20.1], [yy, yy], color=C['faint'], lw=S.LW_AXIS)
    yy += 2.2
    for key, lab in [('Number + position (combined) - Number-only | Mal3vs1', 'DeltaAUC, position added to number')]:
        q = dd.loc[key]
        pg.text(xt, yy, lab, size=S.SZ_SMALL, va='top', color=C['dark'])
        pg.text(xt, yy + 3.55, f'+{q.AUC:.3f}', size=S.SZ_TICK, va='top')
        pg.text(xt + 12.6, yy + 3.75, f"95% CI {S.minus(f'{q.CI_low:.3f}')}-{q.CI_high:.3f}", size=S.SZ_SMALL, va='top', color=C['dark'])
        note('I', f'pre-specified delta AUC, {lab}', f'+{q.AUC:.3f} ({q.CI_low:.3f}-{q.CI_high:.3f})', 'c07.3_r02 ROC data, panel D')
        yy += 8.6
    return ax

def panel_J(pg, x0, y0, w, h):
    q = D.cosmx_deciles()
    m = D.cosmx_models().set_index('outcome')
    med = q.groupby('dec')[['c_cd8T', 'c_cd4T', 'c_mac']].median() * 100
    x = med.index.values + 1
    ax = pg.ax(x0, y0, w, h)
    S.zero_line(ax)
    ax.plot(x, med.c_mac, color=C['neutral'], lw=S.LW_ERR, ls=(0, (2.6, 1.6)), zorder=2)
    ax.plot(x, med.c_cd4T, color=C['dark'], lw=S.LW_ERR, zorder=2)
    ax.plot(x, med.c_cd8T, color=C['CD8T'], lw=S.LW_DATA * 1.5, zorder=3)
    ax.scatter(x, med.c_cd8T, s=10, color=C['CD8T'], lw=0, zorder=4)
    ax.set_xlim(0.6, 10.4)
    ax.set_xticks([1, 5, 10])
    ax.set_ylim(-1.0, 1.4)
    ax.set_yticks([-1, 0, 1])
    ax.set_xlabel('CMI decile')
    ax.set_ylabel('Immune neighbour (Deltapp)')
    tr = mtrans.blended_transform_factory(ax.transAxes, ax.transData)
    specs = [('hit_cd8T', 'CD8 T', C['CD8T'], med.c_cd8T.iloc[-1]), ('hit_cd4T', 'CD4 T', C['dark'], med.c_cd4T.iloc[-1]), ('hit_mac', 'Macrophage', C['neutral'], med.c_mac.iloc[-1])]
    ylab = {'CD8 T': -0.62, 'CD4 T': 0.62, 'Macrophage': 1.18}
    for hit, lab, col, yend in specs:
        b = m.loc[hit]
        b = b[b['item'] == 'CMI_net'].iloc[0]
        ax.text(1.04, ylab[lab], lab, transform=tr, ha='left', va='center', fontsize=S.SZ_TICK, color=col, clip_on=False)
        ax.text(1.04, ylab[lab] - 0.21, S.Ptxt(b.p), transform=tr, ha='left', va='center', fontsize=S.SZ_SMALL, color=col, clip_on=False)
        note('J', f'{lab}: beta per +1 SD CMI_net (pp), P', f'{b.beta * 100:+.2f} pp, P = {S.fmtP(b.p)}', 'c05.1_r01 models 7-9')
    note('J', 'patients with decile curves / in models', f"{q.patient_id.nunique()} / {int(m.loc['hit_cd8T'].n_patient.iloc[0])}", 'c06_r04 / c05.1_r01')
    return ax

def legend_row(pg, x, y):
    items = [('Mal_1', C['Mal_1'], 'sq'), ('Mal_2', C['Mal_2'], 'sq'), ('Mal_3', C['Mal_3'], 'sq'), ('CD8 T', C['CD8T'], 'sq'), ('Other spots', C['spot_bg'], 'sq'), ('Mal_3 territory', C['ink'], 'dash'), ('DLBCL', C['dark'], 'o'), ('PCNSL', C['dark'], 's'), ('MZL', C['dark'], '^'), ('BL', C['dark'], 'D')]
    for t, c, kind in items:
        if kind == 'sq':
            pg.rect(x, y - 1.25, 2.5, 2.5, color=c)
            x0 = x + 3.4
        elif kind == 'dash':
            pg.line([x, x + 5.0], [y, y], color=c, lw=S.LW_ERR, ls=(0, (3, 1.8)))
            x0 = x + 6.0
        else:
            axm = pg.ax(x, y - 1.4, 2.8, 2.8)
            S.blank(axm)
            axm.set_xlim(-1, 1)
            axm.set_ylim(-1, 1)
            axm.scatter([0], [0], s=13 * S.MARK_SCALE[kind] ** 2, marker=kind, facecolor='none', edgecolor=c, lw=0.8)
            x0 = x + 3.6
        pg.text(x0, y, t, size=S.SZ_SMALL, va='center')
        x = x0 + S.text_w(t, S.SZ_SMALL) + (3.0 if t not in ('Mal_3 territory',) else 6.0)
BLOCK_W, BLOCK_H = (180.0, 238.0)

def build(pg):
    L = 0.0
    y = 0.0
    pg.letter(L, y, 'A')
    pg.tag(L + 5.2, y + 1.3, 'Visium | 42 sections from 9 public datasets')
    panel_A(pg, L + 12.5, y + 14.0, 103.0, S.H_PLOT)
    pg.letter(122.0, y, 'B')
    pg.tag(127.2, y + 1.3, 'Visium | all entities pooled')
    panel_B(pg, 136.5, y + 14.0, 42.5, S.H_PLOT)
    y = 62.0
    pg.letter(L, y, 'C')
    pg.tag(L + 5.2, y + 1.3, 'Visium | QC-retained cohort')
    panel_C(pg, L + 14.5, y + 14.0, 30.0, S.H_PLOT)
    pg.letter(43.5, y, 'D')
    pg.tag(48.7, y + 1.3, f"Visium HD | {configured_id('FIG5_HD_SAMPLE_ID')} | 8 um")
    panel_D(pg, 50.0, y + 14.0, 23.0)
    pg.letter(81.5, y, 'E')
    pg.tag(86.7, y + 1.3, 'GeoMx | Cohort_B | 9 paired / 38 patients')
    panel_E(pg, 96.0, y + 14.0, 11.0, S.H_PLOT, 14.0, 137.0, 42.0)
    y = 122.0
    pg.letter(L, y, 'F')
    pg.tag(L + 5.2, y + 1.3, f"Visium | {configured_id('FIG5_SAMPLE_ID')}")
    panel_F(pg, L + 0.5, y + 7.5, 38.0, 38.0, 3.0)
    pg.letter(83.0, y, 'G')
    pg.tag(88.2, y + 1.3, 'Visium | QC-retained cohort')
    panel_G(pg, 97.5, y + 14.0, 31.0, S.H_PLOT)
    pg.letter(133.0, y, 'H')
    pg.tag(138.2, y + 1.3, 'Visium | QC-retained cohort')
    panel_H(pg, 147.5, y + 14.0, 31.5, S.H_PLOT)
    y = 178.0
    pg.letter(L, y, 'I')
    pg.tag(L + 5.2, y + 1.3, 'Visium | QC-retained cohort | 23 patients')
    panel_I(pg, L + 12.5, y + 9.0, S.H_PLOT, L + 49.0)
    pg.letter(96.0, y, 'J')
    pg.tag(101.2, y + 1.3, 'CosMx | Cohort_P | 77 patients (curves 76)')
    panel_J(pg, 110.0, y + 9.0, 44.0, S.H_PLOT)
    legend_row(pg, L + 0.5, 233.5)

def render(out_dir, version):
    os.makedirs(out_dir, exist_ok=True)
    outs = {}
    kinds = {'block': (BLOCK_W, BLOCK_H + 1.0, (0.0, 1.0)), 'A4': (210.0, 297.0, (15.0, 20.0))}
    quick = os.environ.get('QUICK') == '1'
    if quick:
        kinds = {'block': kinds['block']}
    for kind, (w, h, off) in kinds.items():
        LOG.clear()
        pg = S.Page(w, h, off)
        build(pg)
        base = os.path.join(out_dir, f'Figure5_v{version}_{kind}')
        pg.fig.savefig(base + '.pdf')
        if not quick:
            pg.fig.savefig(base + '.svg')
            pg.fig.savefig(base + '.png', dpi=600)
        outs[kind] = (base, pg)
        plt.close(pg.fig)
    pd.DataFrame(LOG).to_csv(os.path.join(out_dir, f'Figure5_v{version}_values.csv'), index=False, encoding='utf-8-sig')
    return outs
if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('SCI_OUTPUT_DIR', os.path.join(HERE, 'results'))
    v = sys.argv[2] if len(sys.argv) > 2 else '2'
    render(out, v)
    print('ok', len(LOG), 'numbers logged')
