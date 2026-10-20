from publication_data import read_csv, read_parquet, numeric_id, numeric_ids, configured_id, validate_cohorts
import os
import argparse, datetime, hashlib, json, os, re, shutil, sys
from pathlib import Path
os.environ.setdefault('FIG6_SZ', '7.0')
os.environ.setdefault('FIG6_H', '247.0')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
HERE = Path(__file__).resolve().parent
STEP = HERE
INDIR = Path(os.environ.get('SCI_INPUT_DIR', str(HERE / 'inputs'))) / 'plot_data'
OUTDIR = Path(os.environ.get('SCI_OUTPUT_DIR', str(HERE / 'results')))
figure_name = 'Figure6'
sys.path.append(str(HERE))
from fig6v2_style import *
LH = SZ * 0.3528 * 1.32

def load():
    D = {f.stem: validate_cohorts(read_csv(f)) for f in INDIR.glob('*.csv')}
    need = ['c20_r00c_design_parameters', 'c20_r01_101_algorithms_Cindex_14_cohorts', 'c20_r01b_RSF_cohort_ranks', 'c20_r02_KM_subject_data', 'c20_r03_KM_cohort_statistics', 'c20_r04_risk_table_counts', 'c20_r05_ROC_coordinates', 'c20_r05c_time_dependent_AUC_95CI', 'c20_r06_published_model_Cindex', 'c20_r06b_published_model_Cindex_cohort_summary', 'c20_r10_clinical_associations_group_statistics', 'c20_r11_clinical_associations_subject_data', 'c20_r12_clinical_analysis_sample_counts', 'c20.1_r01_multivariable_Cox_coefficients', 'c20.1_r02_multivariable_Cox_models']
    miss = [k for k in need if k not in D]
    if miss:
        sys.exit(f'Missing plot data files: {miss} Prepare the required anonymized plot data first.')
    return D
COHORT_ORDER = ['Training', 'Test', 'Cohort_H_DLBCL', 'Cohort_I_DLBCL', 'Cohort_J_DLBCL', 'Cohort_R_DLBCL', 'Cohort_M_DLBCL', 'Cohort_A_DLBCL', 'Cohort_F_DLBCL', 'Cohort_G_FL', 'Cohort_M_BL', 'Cohort_M_MCL', 'Cohort_Q_MCL']
E_ORDER = [g for g in COHORT_ORDER if g not in ('Training', 'Test', 'Cohort_F_DLBCL', 'Cohort_M_BL')]
ENT = {'Mixed': 'Mix'}
MIX_GROUP, MIX_ENT = ('Mix', '(DLBCL, FL, BL, MCL)')
SC_CELLS = 902731

def panel_A(cv, D, x0, y0, w, h):
    cv.letter(x0, y0, 'A')
    k = D['c20_r00c_design_parameters'].set_index('item').value
    steps = [('Single-cell CMI', f'{SC_CELLS:,} B cells'), ('Malignant B-cell genes', f"{int(k['Risk gene (gene_A)'])} risk x {int(k['Protective gene (gene_B)'])} protective"), ('REO gene pairs', f"rank(A) - rank(B) | top {int(k['Selected gene pairs (top 100 by P)'])}"), ('101 ML algorithms', f"{int(k['bulk_cohorts'])} cohorts | {int(k['bulk_tumor_samples']):,} samples"), ('CMI-riskscore (RSF)', '')]
    bx, bw = (x0 + 1.5, w - 1.5)
    top = y0 + 4.6
    gap = 1.35
    bh = (h - 4.6 - gap * (len(steps) - 1)) / len(steps)
    ax = cv.ax(bx, top, bw, len(steps) * (bh + gap) - gap)
    ax.set_xlim(0, bw)
    ax.set_ylim(len(steps) * (bh + gap) - gap, 0)
    ax.axis('off')
    for i, (t1, t2) in enumerate(steps):
        yy = i * (bh + gap)
        last = i == len(steps) - 1
        ax.add_patch(FancyBboxPatch((0.3, yy + 0.15), bw - 0.6, bh - 0.3, boxstyle='round,pad=0,rounding_size=1.0', fc='#FBEAE6' if last else 'white', ec=HIGH if last else MID, lw=0.9 if last else 0.55))
        ax.text(bw / 2, yy + (bh * 0.32 if t2 else bh * 0.5), t1, ha='center', va='center', fontsize=SZ, color=HIGH if last else 'black')
        if t2:
            ax.text(bw / 2, yy + bh * 0.72, t2, ha='center', va='center', fontsize=SZ, color=TXT)
        if not last:
            ax.annotate('', xy=(bw / 2, yy + bh + gap - 0.05), xytext=(bw / 2, yy + bh - 0.12), arrowprops=dict(arrowstyle='-|>', lw=0.6, color=MID, mutation_scale=5))

def polyline(cv, pts, **kw):
    from matplotlib.lines import Line2D
    X = [(cv.ox + x) / cv.fw for x, _ in pts]
    Y = [1 - (cv.oy + y) / cv.fh for _, y in pts]
    return cv.fig.add_artist(Line2D(X, Y, transform=cv.fig.transFigure, solid_capstyle='butt', **kw))

def panel_B(cv, D, x0, y0, pw, ph):
    cv.letter(x0, y0, 'B')
    cl = D['c20_r01_101_algorithms_Cindex_14_cohorts']
    km = D['c20_r03_KM_cohort_statistics'].set_index('group')
    rs = D['c20_r01b_RSF_cohort_ranks']
    cl = cl[cl.group.isin(COHORT_ORDER)]
    n_alg = cl.groupby('group').model_id.nunique()
    if not (n_alg == 101).all():
        sys.exit(f'Each cohort must have C-index values for 101 algorithms; found {sorted(set(n_alg))}')
    px, py = (x0 + 11.0, y0 + 10.6)
    cv.label(x0 + 5.0, y0 + 0.5, f"{int(n_alg.iloc[0])} algorithms | {len(COHORT_ORDER)} cohorts | n = {int(km.loc[COHORT_ORDER, 'n'].sum()):,} | {int(km.loc[COHORT_ORDER, 'events'].sum()):,} deaths")
    assert int(rs.n_display_cohorts.iloc[0]) == len(COHORT_ORDER)
    ax = cv.ax(px, py, pw, ph, row='B')
    xlo, xhi = (-0.6, len(COHORT_ORDER) - 0.4)
    X = lambda d: px + (d - xlo) * pw / (xhi - xlo)
    xs = {g: i for i, g in enumerate(COHORT_ORDER)}
    ent = {'DLBCL': (2, 8), 'FL': (9, 9), 'BL': (10, 10), 'MCL': (11, 12)}
    if COHORT_ORDER[:2] != ['Training', 'Test'] or km.loc[COHORT_ORDER[2:], 'entity'].eq('Mixed').any():
        sys.exit('The first two columns must be training and test sets; the remaining 11 must be disease-specific cohorts.')
    y_ent = py - 0.4
    for k, (a, b) in ent.items():
        if k in ('DLBCL', 'BL'):
            ax.add_patch(Rectangle((a - 0.5, 0.35), b - a + 1, 0.7, fc=BAND, ec='none', zorder=0))
        cv.text(X((a + b) / 2), y_ent, k, ha='center', va='bottom', color=TXT)
    cv.text(X(1.38), y_ent, MIX_ENT, ha='right', va='bottom', color=TXT)
    y_br = y_ent - LH - 0.5
    for (a, b), lab in (((-0.4, 1.4), MIX_GROUP), ((1.5, 12.5), 'Validate')):
        polyline(cv, [(X(a), y_br + 0.9), (X(a), y_br), (X(b), y_br), (X(b), y_br + 0.9)], color=TXT, lw=0.5)
        cv.text(X((a + b) / 2), y_br - 0.35, lab, ha='center', va='bottom', color=TXT)
    rng = np.random.default_rng(20261011)
    oth = cl[~cl.is_RSF]
    ax.scatter(oth.group.map(xs) + rng.uniform(-0.3, 0.3, len(oth)), oth.cindex, s=1.8, c=DOT, lw=0, alpha=0.85, zorder=2)
    r = cl[cl.is_RSF].set_index('group').loc[COHORT_ORDER]
    xr = np.arange(len(COHORT_ORDER))
    ax.plot(xr, r.cindex, color=HIGH, lw=LW_DATA, zorder=3)
    ax.scatter(xr, r.cindex, s=11, c=HIGH, edgecolors='black', linewidths=0.4, zorder=4)
    ax.text(0, r.cindex.iloc[0] + 0.045, 'RSF', ha='center', va='bottom', fontsize=SZ, color=HIGH)
    mean = float(rs.rsf_mean.iloc[0])
    ax.axhline(mean, color=HIGH, lw=0.55, ls=(0, (3, 2)), zorder=1)
    ax.text(len(COHORT_ORDER) - 0.45, mean + 0.012, f'mean {mean:.3f}', ha='right', va='bottom', fontsize=SZ, color=HIGH)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(0.35, 1.05)
    ax.set_yticks([0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.4', '0.6', '0.8', '1.0'])
    ax.set_ylabel('C-index')
    labs = [f"{km.loc[g, 'label']} ({int(km.loc[g, 'n'])})" for g in COHORT_ORDER]
    ax.set_xticks(xr)
    ax.set_xticklabels(labs, rotation=40, ha='right', rotation_mode='anchor')
    ax.tick_params(axis='x', length=1.4)
    return ax
RT_GAP = 6.6

def km_panel(cv, D, g, x, y, w, h, show_ylabel=True):
    pts = D['c20_r02_KM_subject_data']
    st = D['c20_r03_KM_cohort_statistics'].set_index('group').loc[g]
    rt = D['c20_r04_risk_table_counts']
    rt = rt[rt.group == g]
    ax = cv.ax(x, y, w, h, row='C')
    curves = []
    for grp, col in (('Low', LOW), ('High', HIGH)):
        s = pts[(pts.group == g) & (pts.risk_group == grp)]
        X, Y, cx, cy = km_steps(s.os_months, s.os_status, min_risk=MIN_RISK)
        curves.append((X, Y))
        ax.plot(X, Y, color=col, lw=LW_DATA, zorder=3)
        ax.plot(cx, cy, ls='none', marker='|', ms=2.0, mew=0.4, color=col, zorder=3)
    tmax = 160
    ax.set_xlim(0, tmax)
    ax.set_ylim(0, 1.03)
    ax.set_xticks([0, 50, 100, 150])
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(['0', '0.5', '1'])
    ax.set_xlabel('Months', labelpad=0.8)
    if show_ylabel:
        ax.set_ylabel('Overall survival')
    l1 = f'HR {fmt_ci(st.HR_show)} ({fmt_ci(st.HR_lower_show)}-{fmt_ci(st.HR_upper_show)})'
    l2 = f'P = {fmtP(st.logrank_p)}'
    if not place_lines(ax, [[l1, l2]], curves, x_right=tmax * 0.98, x_left=4, y_cands=[0.3, 0.22, 0.38, 0.46, 0.14], line_mm=LH, pad_mm=0.4, fontsize=SZ):
        sys.exit(f'{g}: no non-overlapping HR and P label position.')
    rax = cv.ax(x, y + h + RT_GAP, w, 2 * LH)
    rax.set_xlim(0, tmax)
    rax.set_ylim(0, 2)
    rax.axis('off')
    for row, (grp, col) in enumerate((('High', HIGH), ('Low', LOW))):
        r = rt[rt.risk_group == grp].sort_values('time')
        for t, n in zip(r.time, r.n_at_risk):
            rax.text(t, 1.5 - row, f'{int(n)}', ha='center', va='center', fontsize=SZ, color=col)
    return ax

def roc_panel(cv, D, g, x, y, w, h, show_ylabel=True):
    roc = D['c20_r05_ROC_coordinates']
    roc = roc[roc.group == g]
    ci = D['c20_r05c_time_dependent_AUC_95CI']
    ci = ci[ci.group == g].set_index('eval_time')
    ax = cv.ax(x, y, w, h, row='C')
    ax.plot([0, 1], [0, 1], color=MID, lw=0.45, ls=(0, (2, 2)), zorder=1)
    for i, t in enumerate((12, 36, 60)):
        r = roc[roc.eval_time == t]
        ax.plot(r.fpr, r.tpr, color=ROC3[t], lw=LW_DATA * 0.9, zorder=3 + i)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels(['0', '0.5', '1'])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0', '0.5', '1'])
    ax.set_xlabel('1 - specificity', labelpad=0.8)
    if show_ylabel:
        ax.set_ylabel('Sensitivity')
    y_t = y + h + RT_GAP + LH * 0.5
    for i, t in enumerate((12, 36, 60)):
        s = ci.loc[t]
        yy = y_t + i * LH
        cv.hline(x - 0.5, x + 1.6, yy, color=ROC3[t], lw=LW_DATA * 1.4, solid_capstyle='butt')
        cv.text(x + 2.2, yy, f'{t // 12}-y AUC {s.auc:.2f} ({s.lower:.2f}-{s.upper:.2f})', ha='left', va='center')
    return ax

def panel_C(cv, D, x0, y0, plot_left, kh):
    cv.letter(x0, y0, 'C')
    st = D['c20_r03_KM_cohort_statistics'].set_index('group')
    kw, rw = (35.0, 24.0)
    xs = [plot_left, plot_left + kw + 9.0]
    for g, xx in zip(('Training', 'Test'), xs):
        s = st.loc[g]
        cv.label(x0 + 5.0 if xx == xs[0] else xx - 1.0, y0 + 0.5, f'{s.label} | n = {int(s.n)} | {int(s.events)} deaths')
        km_panel(cv, D, g, xx, y0 + 4.6, kw, kh, show_ylabel=xx == xs[0])
    for row, (nm, col) in enumerate((('CMI-high', HIGH), ('CMI-low', LOW))):
        cv.text(xs[0] - 3.2, y0 + 4.6 + kh + RT_GAP + LH * (0.5 + row), nm, ha='right', va='center', color=col)
    rx = [xs[1] + kw + 12.5, xs[1] + kw + 12.5 + rw + 10.5]
    for g, xx in zip(('Training', 'Test'), rx):
        cv.label(xx - (8.6 if xx == rx[0] else 1.2), y0 + 0.5, f"{st.loc[g, 'label']}")
        roc_panel(cv, D, g, xx, y0 + 4.6, rw, kh, show_ylabel=xx == rx[0])
    return rx[1] + rw

def fmt_ci(v):
    if v >= 1000:
        return f"{float(f'{v:.3g}'):,.0f}"
    if v >= 100:
        return f'{v:.0f}'
    return f'{v:#.3g}'.rstrip('.')

def cohort_labels(lab, st, order, w):
    prev = None
    for i, g in enumerate(order):
        e = ENT.get(st.loc[g, 'entity'], st.loc[g, 'entity'])
        lab.text(0.2, i, '' if e == prev else e, ha='left', va='center', fontsize=SZ, color=TXT)
        prev = e
        lab.text(w - 0.2, i, f"{st.loc[g, 'label']} ({int(st.loc[g, 'n'])})", ha='right', va='center', fontsize=SZ)

def entity_rules(ax, st, order):
    for i in range(len(order) - 1):
        if st.loc[order[i], 'entity'] != st.loc[order[i + 1], 'entity']:
            ax.axhline(i + 0.5, color=MID, lw=0.3, zorder=0)

def panel_D(cv, D, x0, y0, rows_y, rows_h):
    cv.letter(x0, y0, 'D')
    st = D['c20_r03_KM_cohort_statistics'].set_index('group')
    order = COHORT_ORDER
    n = len(order)
    cv.label(x0 + 5.0, y0 + 0.5, f"{n} cohorts | n = {int(st.loc[order, 'n'].sum()):,} | {int(st.loc[order, 'events'].sum()):,} deaths")
    wl, wf, wt = (33.0, 22.0, 35.0)
    lab = cv.ax(x0, rows_y, wl, rows_h)
    lab.axis('off')
    lab.set_xlim(0, wl)
    lab.set_ylim(n - 0.5, -0.5)
    cohort_labels(lab, st, order, wl)
    ax = cv.ax(x0 + wl + 1.0, rows_y, wf, rows_h, row='DE')
    ax.patch.set_alpha(0)
    lo_x, hi_x = (0.8, 1000.0)
    for i, g in enumerate(order):
        s = st.loc[g]
        a, b = (max(s.HR_lower_show, lo_x), min(s.HR_upper_show, hi_x))
        ax.plot([a, b], [i, i], color=HIGH, lw=LW_ERR, zorder=2, solid_capstyle='butt')
        if s.HR_upper_show > hi_x:
            ax.annotate('', xy=(hi_x, i), xytext=(hi_x * 0.62, i), arrowprops=dict(arrowstyle='-|>', lw=LW_ERR, color=HIGH, mutation_scale=4.5), zorder=2)
        ax.scatter([s.HR_show], [i], s=13, marker='s', c=HIGH, lw=0, zorder=3)
    ax.axvline(1, color=MID, lw=0.45, ls=(0, (2, 2)), zorder=0)
    entity_rules(ax, st, order)
    ax.set_xscale('log')
    ax.set_xlim(lo_x, hi_x)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(['1', '10', '100', '1,000'])
    ax.minorticks_off()
    ax.tick_params(left=False, labelleft=False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Hazard ratio')
    num = cv.ax(x0 + wl + wf + 2.0, rows_y, wt, rows_h)
    num.axis('off')
    num.set_xlim(0, wt)
    num.set_ylim(n - 0.5, -0.5)
    num.text(20.5, -1.0, 'HR (95% CI)', ha='right', va='center', fontsize=SZ, color=TXT)
    num.text(wt - 0.2, -1.0, 'P', ha='right', va='center', fontsize=SZ, color=TXT)
    for i, g in enumerate(order):
        s = st.loc[g]
        t1 = f'{fmt_ci(s.HR_show)} ({fmt_ci(s.HR_lower_show)}-{fmt_ci(s.HR_upper_show)})'
        t2 = fmtP(s.logrank_p)
        num.text(20.5, i, t1, ha='right', va='center', fontsize=SZ)
        num.text(wt - 0.2, i, t2, ha='right', va='center', fontsize=SZ)
        cv.tables.append({'panel': 'D', 'texts': [f'{s.label} ({int(s.n)})', t1, t2]})
    return x0 + wl + wf + 2.0 + wt

def panel_E(cv, D, x0, y0, rows_y, rows_h, w):
    cv.letter(x0, y0, 'E')
    c = D['c20_r06_published_model_Cindex']
    summ = D['c20_r06b_published_model_Cindex_cohort_summary'].set_index('group')
    if set(E_ORDER) != set(summ.index[summ.show_in_figure]):
        sys.exit(f'Displayed cohorts do not match the input display flags: {sorted(set(E_ORDER) ^ set(summ.index[summ.show_in_figure]))}')
    c = c[c.group.isin(E_ORDER)]
    n_pub = int(c[~c.is_RSF].model.nunique())
    cv.label(x0 + 5.0, y0 + 0.5, f"{len(E_ORDER)} cohorts | n = {int(summ.loc[E_ORDER, 'n'].sum()):,} | {int(summ.loc[E_ORDER, 'events'].sum()):,} deaths")
    n = len(E_ORDER)
    wl, wq = (30.0, 13.5)
    wp = w - wl - wq - 2.0
    lab = cv.ax(x0, rows_y, wl, rows_h)
    lab.axis('off')
    lab.set_xlim(0, wl)
    lab.set_ylim(n - 0.5, -0.5)
    cohort_labels(lab, summ.assign(entity=summ.entity), E_ORDER, wl)
    ax = cv.ax(x0 + wl + 1.0, rows_y, wp, rows_h, row='DE')
    yi = {g: i for i, g in enumerate(E_ORDER)}
    lo, hi = (c.cindex.min(), c.cindex.max())
    if lo < 0.3 or hi > 1.0:
        sys.exit(f'C-index {lo:.3f}-{hi:.3f} exceeds the axis range and would be clipped.')
    pub = c[~c.is_RSF]
    rsf = c[c.is_RSF]
    ax.scatter(pub.cindex, pub.group.map(yi), s=4, c=DOT, lw=0, alpha=0.95, zorder=2)
    ax.scatter(rsf.cindex, rsf.group.map(yi), s=13, marker='D', c=HIGH, edgecolors='black', linewidths=0.35, zorder=4)
    for i in range(n):
        ax.axhline(i, color='#E6E7E9', lw=0.35, zorder=0)
    entity_rules(ax, summ, E_ORDER)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0.3, 1.0)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(['0.4', '0.6', '0.8', '1.0'])
    ax.tick_params(left=False, labelleft=False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('C-index')
    k_pub = pub.groupby(pub.group.map(summ.entity)).model.nunique()
    ents = [e for e in ('DLBCL', 'FL', 'MCL') if e in k_pub.index]
    if len(ents) != len(k_pub):
        sys.exit(f'Unexpected disease groups in model comparisons: {list(k_pub.index)}')
    key(lab, 0.9, -1.0, 'D', HIGH, 'CMI-riskscore', size=11, edge='black', dx=1.3)
    key(lab, 20.5, -1.0, 'o', DOT, 'published models: ' + ' | '.join((f'{e} {int(k_pub[e])}' for e in ents)), size=8, dx=1.3)
    per = pub.groupby('group').model.nunique()
    if any((per[g] != k_pub[summ.loc[g, 'entity']] for g in E_ORDER)):
        sys.exit('Published model counts must agree across cohorts of the same disease.')
    qa = cv.ax(x0 + wl + wp + 2.0, rows_y, wq, rows_h)
    qa.axis('off')
    qa.set_xlim(0, wq)
    qa.set_ylim(n - 0.5, -0.5)
    qa.text(wq - 0.2, -1.0, 'q vs best', ha='right', va='center', fontsize=SZ, color=TXT)
    for i, g in enumerate(E_ORDER):
        t = fmtP(summ.loc[g, 'q_best'])
        qa.text(wq - 0.2, i, t, ha='right', va='center', fontsize=SZ)
        cv.tables.append({'panel': 'E', 'texts': [f"{summ.loc[g, 'label']} ({int(summ.loc[g, 'n'])})", t]})
    cv.n_points['E'] = len(c)

def panel_F(cv, D, x0, y0, h, plot_left, w_each, gap):
    cv.letter(x0, y0, 'F')
    rows = D['c20_r10_clinical_associations_group_statistics']
    pts = D['c20_r11_clinical_associations_subject_data']
    scale = D['c20_r12_clinical_analysis_sample_counts'].set_index('disease')
    keys = rows[rows.panel == 'F'].drop_duplicates('row_order').sort_values('row_order')
    ents = list(dict.fromkeys(keys.disease))
    n_eff = lambda d: int(scale.loc[d, 'n_scored'] - scale.loc[d, 'n_QC_excluded'])
    cv.label(x0 + 5.0, y0 + 0.5, f'{len(ents)} lymphoma entities | ' + ' | '.join((f'{d} {n_eff(d):,}' for d in ents)) + ' scored tumours')
    y0 = y0 + LH + 0.4
    ylo, yhi = (-3.3, 8.0)
    head = 2 * LH + 0.6
    ph = h - (LH + 0.4) - head - 3 * LH - 1.2
    axes = []
    for k, (_, r) in enumerate(keys.iterrows()):
        x = plot_left + k * (w_each + gap)
        g = rows[rows.row_order == r.row_order].sort_values('order')
        cv.text(x + w_each / 2, y0 + 0.5, f'{r.disease}', ha='center', color=TXT)
        t = cv.text(x + w_each / 2, y0 + 0.5 + LH, r.row_label, ha='center')
        bb = t.get_window_extent(cv.fig.canvas.get_renderer())
        over = bb.x1 * 25.4 / cv.fig.dpi - cv.ox - (W - 0.8)
        if over > 0:
            t.set_x(t.get_position()[0] - over / cv.fw)
        ax = cv.ax(x, y0 + head, w_each, ph, row='F')
        lv = g.level_label.tolist()
        cols = [LOW, HIGH] if len(lv) == 2 else [LOW] + [DOT] * (len(lv) - 2) + [HIGH]
        data = [pts[(pts.row_order == r.row_order) & (pts.level_label == L)].score.values for L in lv]
        pos = violin_box(ax, data, cols, seed=20261011 + k)
        top = max((np.max(v) for v in data))
        p_txt = fmtP(float(g.p_show.iloc[0]))
        p_bracket(ax, pos[0], pos[-1], top + 0.45, f'P = {p_txt}')
        ax.set_xlim(0.4, len(lv) + 0.6)
        ax.set_ylim(ylo, yhi)
        ax.set_yticks([-3, 0, 3, 6])
        if k == 0:
            ax.set_ylabel('CMI-riskscore (z)')
        else:
            ax.tick_params(labelleft=False)
        ax.set_xticks(pos)
        ax.set_xticklabels([f'{L}\n{int(n):,}' for L, n in zip(lv, g.n)], linespacing=1.15)
        ax.tick_params(axis='x', length=1.2)
        axes.append(ax)
    return axes
COX_LAYOUT = [['DLBCL'], ['FL', 'MCL', 'BL']]
COX_DATASET = {'DLBCL': 'Cohort_O', 'FL': 'Cohort_G', 'MCL': 'Cohort_Q', 'BL': 'Cohort_O'}

def cox_column(cv, coef, mod, diseases, x, y, p, n_max, wl, wf, wt, axis_lo=0.5, axis_hi=200.0, gap=1.2):
    lines = []
    for dz in diseases:
        lines.append(('head', dz, None))
        for _, r in coef[coef.disease == dz].sort_values('row_order').iterrows():
            lines.append(('row', dz, r))
        lines.append(('cidx', dz, None))
    n = len(lines)
    assert n <= n_max
    pc = p * (n_max - 1) / (n - 1)
    y_first = y + 0.5 * p
    y_bot = y_first + (n - 1) * pc + 0.5 * p + gap
    hh = y_bot - y
    ylim = ((y_bot - y_first) / pc, -(y_first - y) / pc)
    lab = cv.ax(x, y, wl, hh)
    lab.axis('off')
    lab.set_xlim(0, wl)
    lab.set_ylim(*ylim)
    ax = cv.ax(x + wl + 1.0, y, wf, hh, row='G')
    ax.patch.set_alpha(0)
    num = cv.ax(x + wl + wf + 2.0, y, wt, hh)
    num.axis('off')
    num.set_xlim(0, wt)
    num.set_ylim(*ylim)
    yline = lambda i: y_first + i * pc
    cv.text(x + wl + wf + 2.0 + wt - 13.4, y_first - p, 'HR (95% CI)', ha='right', va='center', color=TXT)
    cv.text(x + wl + wf + 2.0 + wt - 0.2, y_first - p, 'P', ha='right', va='center', color=TXT)
    blocks = {}
    for i, (kind, dz, r) in enumerate(lines):
        m = mod.loc[dz]
        if kind == 'head':
            cv.text(x + 0.2, yline(i), f'{dz} | {COX_DATASET[dz]} | n = {int(m.n)} | {int(m.events)} deaths', ha='left', va='center', color=TXT)
            blocks.setdefault(dz, [i + 1, None])
        elif kind == 'cidx':
            blocks[dz][1] = i - 1
            t = f'C-index {m.C_clinical:.2f} -> {m.C_full:.2f} with CMI-riskscore (P = {fmtP(m.dC_p)})'
            cv.text(x + 0.2, yline(i), t, ha='left', va='center', color=TXT)
        else:
            col = HIGH if r.term == 'cmi' else 'black'
            lab.text(wl - 0.2, i, r.label, ha='right', va='center', fontsize=SZ, color=col)
            a, b = (max(r.lower, axis_lo), min(r.upper, axis_hi))
            ax.plot([a, b], [i, i], color=HIGH if r.term == 'cmi' else '#55595E', lw=LW_ERR, zorder=2, solid_capstyle='butt')
            if r.upper > axis_hi:
                ax.annotate('', xy=(axis_hi, i), xytext=(axis_hi * 0.6, i), arrowprops=dict(arrowstyle='-|>', lw=LW_ERR, color=HIGH, mutation_scale=4.5), zorder=2)
            ax.scatter([r.HR], [i], s=13 if r.term == 'cmi' else 8, marker='s', c=col if r.term == 'cmi' else '#55595E', lw=0, zorder=3)
            t1 = f'{fmt_ci(r.HR)} ({fmt_ci(r.lower)}-{fmt_ci(r.upper)})'
            t2 = fmtP(r.p)
            num.text(wt - 13.4, i, t1, ha='right', va='center', fontsize=SZ, color=col)
            num.text(wt - 0.2, i, t2, ha='right', va='center', fontsize=SZ, color=col)
            cv.tables.append({'panel': f'G-{diseases[0]}column', 'texts': [r.label, t1, t2]})
    for i, (kind, dz, r) in enumerate(lines[:-1]):
        if kind == 'cidx':
            cv.hline(x, x + wl + wf + 2.0 + wt, yline(i) + 0.5 * pc, color=MID, lw=0.3)
    for dz, (a, b) in blocks.items():
        ax.plot([1, 1], [a - 0.45, b + 0.45], color=MID, lw=0.45, ls=(0, (2, 2)), zorder=0)
    ax.set_xscale('log')
    ax.set_xlim(axis_lo, axis_hi)
    ax.set_ylim(*ylim)
    ax.set_xticks([1, 10, 100])
    ax.set_xticklabels(['1', '10', '100'])
    ax.minorticks_off()
    ax.tick_params(left=False, labelleft=False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Hazard ratio')
    return y_bot

def panel_G(cv, D, x0, y0, h, col_w):
    cv.letter(x0, y0, 'G')
    coef = D['c20.1_r01_multivariable_Cox_coefficients']
    mod = D['c20.1_r02_multivariable_Cox_models'].set_index('disease')
    cv.label(x0 + 5.0, y0 + 0.5, f'{len(mod)} lymphoma entities | {len(set(COX_DATASET.values()))} cohorts | n = {int(mod.n.sum()):,} | {int(mod.events.sum()):,} deaths')
    top = y0 + 4.6 + LH
    nmax = max((sum((1 + (coef.disease == dz).sum() + 1 for dz in col)) for col in COX_LAYOUT))
    gap = 1.2
    p = (h - (top - y0) - 7.0 - gap) / nmax
    bots = [cox_column(cv, coef, mod, col, xx, top, p, nmax, wl=31.0, wf=19.0, wt=col_w - 31.0 - 19.0 - 3.0, gap=gap) for col, xx in zip(COX_LAYOUT, (x0, x0 + col_w + 4.0))]
    assert max(bots) - min(bots) < 1e-06, bots

def draw(fig, fw, fh, ox, oy, D):
    cv = Canvas(fig, fw, fh, ox, oy)
    cv.tables = []
    cv.n_points = {}
    PL = 15.0
    y = 0.0
    panel_A(cv, D, 0.0, y, 47.0, 41.0)
    panel_B(cv, D, 50.0, y, 117.0, 15.4)
    y = 43.5
    panel_C(cv, D, 0.0, y, PL, 20.0)
    y = 89.0
    rows_y, rows_h = (y + 8.2, 13 * LH * 1.05)
    xr = panel_D(cv, D, 0.0, y, rows_y, rows_h)
    panel_E(cv, D, xr + 3.0, y, rows_y, rows_h, W - xr - 3.0)
    y = rows_y + rows_h + 9.0
    panel_F(cv, D, 0.0, y, 40.0, PL, 23.4, 4.6)
    y = y + 41.5
    panel_G(cv, D, 0.0, y, H - y, (W - 4.0) / 2)
    return cv

def check_alignment(cv):
    import collections
    rows = collections.defaultdict(list)
    for r, x, yb in cv.reg:
        rows[r].append((x, yb))
    out = {'row_baseline_deviation_mm': {}, 'left_edge_deviation_mm': {}}
    for r, v in rows.items():
        ys = [b for _, b in v]
        out['row_baseline_deviation_mm'][r] = round(max(ys) - min(ys), 4)
    lefts = {r: min((x for x, _ in v)) for r, v in rows.items()}
    out['left_edge_deviation_mm']['C and F leftmost plotting regions'] = round(abs(lefts['C'] - lefts['F']), 4)
    bad = {k: v for k, v in {**out['row_baseline_deviation_mm'], **out['left_edge_deviation_mm']}.items() if v > 1e-06}
    if bad:
        sys.exit(f'Alignment validation failed: {bad}')
    return out

def next_version(d: Path, name: str) -> int:
    pat = re.compile(f'^{re.escape(name)}_v(\\d+)\\.(pdf|svg|png)$')
    vs = [int(m.group(1)) for f in d.iterdir() if (m := pat.match(f.name))]
    return max(vs) + 1 if vs else 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', default='submission', choices=['submission', 'internal'])
    ap.add_argument('--version', type=int, default=None, help='Automatically increment the version.')
    ap.add_argument('--preview', default=None, help='Save a preview PNG.')
    a = ap.parse_args()
    fam = setup_fonts()
    D = load()
    MM = 1 / 25.4
    if a.preview:
        fig = plt.figure(figsize=(W * MM, H * MM))
        draw(fig, W, H, 0.0, 0.0, D)
        fig.savefig(a.preview, dpi=300)
        plt.close(fig)
        print(f'Preview: {a.preview}')
        return
    OUTDIR.mkdir(parents=True, exist_ok=True)
    v = a.version or next_version(OUTDIR, figure_name)
    fig = plt.figure(figsize=(W * MM, H * MM))
    cv = draw(fig, W, H, 0.0, 0.0, D)
    alignment = check_alignment(cv)
    highlights = highlight_artist_count(fig, HIGH)
    if highlights == 0:
        sys.exit('The highlight color is absent from the figure.')
    if a.mode == 'internal':
        fig.text(0.995, 0.002, f'v{v} | {datetime.date.today()}', ha='right', va='bottom', fontsize=6, color=MID)
    base = OUTDIR / f'{figure_name}_v{v}'
    fig.savefig(f'{base}.pdf')
    fig.savefig(f'{base}.svg')
    fig.savefig(f'{base}.png', dpi=600)
    plt.close(fig)
    fig = plt.figure(figsize=(A4W * MM, A4H * MM))
    draw(fig, A4W, A4H, A4_L, A4_T, D)
    base4 = OUTDIR / f'{figure_name}_A4_v{v}'
    fig.savefig(f'{base4}.pdf')
    fig.savefig(f'{base4}.png', dpi=300)
    plt.close(fig)
    rec = {'version': f'v{v}', 'figure_name': figure_name, 'mode': a.mode, 'font': fam, 'canvas_mm': [W, H], 'a4_page_mm': [A4W, A4H], 'font_size_pt': {'ticks': SZ, 'axis_titles': SZ_T, 'panel_letter': SZ_TAG}, 'panel': ['A study design', 'B 101 algorithms in 13 cohorts', 'C training and test KM and ROC', 'D HR in 13 cohorts', 'E C-index comparison with 34 published models', 'F clinical and molecular features', 'G disease-specific multivariable Cox models'], 'highlight_color': HIGH, 'highlight_used': highlights > 0, 'highlight_artist_count': highlights, 'check_alignment': alignment, 'table_rows': cv.tables, 'panel_E_point_counts': cv.n_points}
    (OUTDIR / f'{figure_name}_v{v}_validation.json').write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[OK] {os.path.basename(str(base))}.pdf / .svg / .png; A4 {base4}.pdf ;font {fam}; font size {SZ} pt; canvas {W:.0f} x {H:.0f} mm')
if __name__ == '__main__':
    main()
