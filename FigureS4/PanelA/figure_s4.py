from publication_data import read_csv, read_parquet, numeric_id, numeric_ids, configured_id, validate_cohorts
import os
import argparse, datetime, hashlib, json, re, shutil, sys
from pathlib import Path
import numpy as np
import pandas as pd
HERE = Path(__file__).resolve().parent
STEP = HERE
INDIR = Path(os.environ.get('SCI_INPUT_DIR', str(HERE / 'inputs'))) / 'plot_data'
OUTDIR = Path(os.environ.get('SCI_OUTPUT_DIR', str(HERE / 'results')))
figure_name = 'FigureS4'
sys.path.append(str(HERE))
from fig6v2_style import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
LH = SZ * 0.3528 * 1.32
D_ORDER = ['Cohort_H_DLBCL', 'Cohort_I_DLBCL', 'Cohort_J_DLBCL', 'Cohort_R_DLBCL', 'Cohort_M_DLBCL', 'Cohort_A_DLBCL', 'Cohort_F_DLBCL', 'Cohort_G_FL', 'Cohort_M_BL', 'Cohort_M_MCL', 'Cohort_Q_MCL']
EF_ORDER = ['Cohort_H_DLBCL', 'Cohort_I_DLBCL', 'Cohort_J_DLBCL', 'Cohort_R_DLBCL', 'Cohort_M_DLBCL', 'Cohort_A_DLBCL', 'Cohort_G_FL', 'Cohort_M_MCL', 'Cohort_Q_MCL']
LABEL = {'Training': ('Mix', 'Train'), 'Cohort_H_DLBCL': ('DLBCL', 'Cohort_H'), 'Cohort_I_DLBCL': ('DLBCL', 'Cohort_I'), 'Cohort_M_DLBCL': ('DLBCL', 'Cohort_M'), 'Cohort_F_DLBCL': ('DLBCL', 'Cohort_F'), 'Cohort_A_DLBCL': ('DLBCL', 'Cohort_A'), 'Cohort_J_DLBCL': ('DLBCL', 'Cohort_J'), 'Cohort_R_DLBCL': ('DLBCL', 'Cohort_R'), 'Cohort_G_FL': ('FL', 'Cohort_G'), 'Cohort_M_BL': ('BL', 'Cohort_M'), 'Cohort_M_MCL': ('MCL', 'Cohort_M'), 'Cohort_Q_MCL': ('MCL', 'Cohort_Q')}

def load():
    D = {f.stem: validate_cohorts(read_csv(f)) for f in INDIR.glob('*.csv')}
    need = ['c20_r02_KM_subject_data', 'c20_r03_KM_cohort_statistics', 'c20_r04_risk_table_counts', 'c20_r05c_time_dependent_AUC_95CI', 'c20_r07_published_model_integrated_AUC', 'c20_r07c_published_model_integrated_AUC_cohort_summary', 'c20_r08_time_dependent_AUC_long', 'c20_r10_clinical_associations_group_statistics', 'c20_r11_clinical_associations_subject_data', 'c20_r12_clinical_analysis_sample_counts', 'c20.1_r01_multivariable_Cox_coefficients', 'c20.1_r02_multivariable_Cox_models', 'c20.1_r03_multivariable_Cox_sensitivity']
    miss = [k for k in need if k not in D]
    if miss:
        sys.exit(f'Missing plot data: {miss}')
    return D

def fmt_ci(v):
    if v >= 1000:
        return f"{float(f'{v:.3g}'):,.0f}"
    if v >= 100:
        return f'{v:.0f}'
    return f'{v:#.3g}'.rstrip('.')

def page1(fig, fw, fh, ox, oy, D):
    cv = Canvas(fig, fw, fh, ox, oy)
    cv.tables = []
    overlaps = []
    pts = D['c20_r02_KM_subject_data']
    st = D['c20_r03_KM_cohort_statistics'].set_index('group')
    rt = D['c20_r04_risk_table_counts']
    ntot = int(st.loc[D_ORDER, 'n'].sum())
    etot = int(st.loc[D_ORDER, 'events'].sum())
    cv.letter(0, 0, 'A')
    cv.label(5.5, 0.6, f'{len(D_ORDER)} cohorts | n = {ntot:,} | {etot:,} deaths')
    ncol, pw, ph = (4, 33.5, 26.0)
    cell_w, cell_h = (44.5, 48.5)
    left, top0 = (12.0, 10.0)
    for i, g in enumerate(D_ORDER):
        r, c = divmod(i, ncol)
        x = left + c * cell_w
        y = top0 + r * cell_h
        s = st.loc[g]
        ent, lab = LABEL[g]
        cv.text(x, y - 3.6, f'{ent} {lab} ({int(s.n)})')
        ax = cv.ax(x, y, pw, ph, row=f'S6A-{r + 1}')
        curves = []
        for grp, col in (('Low', LOW), ('High', HIGH)):
            q = pts[(pts.group == g) & (pts.risk_group == grp)]
            X, Y, cx, cy = km_steps(q.os_months, q.os_status, min_risk=MIN_RISK)
            if X.max() > 120:
                k = X <= 120
                X, Y = (np.append(X[k], 120.0), np.append(Y[k], Y[k][-1]))
            curves.append((X, Y))
            ax.plot(X, Y, color=col, lw=LW_DATA * 0.9, zorder=3)
            ck = cx <= min(120, X.max())
            ax.plot(cx[ck], cy[ck], ls='none', marker='|', ms=2.0, mew=0.4, color=col, zorder=3)
        ax.set_xlim(0, 120)
        ax.set_ylim(0, 1.04)
        ax.set_xticks([0, 24, 48, 72, 96, 120])
        ax.set_xticklabels(['0', '', '48', '', '96', ''])
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(['0', '0.5', '1'])
        if c != 0:
            ax.tick_params(labelleft=False)
        hr = f'HR {fmt_ci(s.HR_show)} ({fmt_ci(s.HR_lower_show)}-{fmt_ci(s.HR_upper_show)})'
        l2 = f'P = {fmtP(s.logrank_p)}'
        variants = [[hr, l2], [f'HR {fmt_ci(s.HR_show)}', f'({fmt_ci(s.HR_lower_show)}-{fmt_ci(s.HR_upper_show)})', l2]]
        ok = place_lines(ax, variants, curves, x_right=118, x_left=4, y_cands=[0.4, 0.3, 0.5, 0.2, 0.6, 0.1, 0.04, 0.68], pad_mm=0.45, fontsize=SZ)
        if not ok:
            overlaps.append(g)
        rax = cv.ax(x, y + ph + 4.2, pw, 5.8)
        rax.set_xlim(0, 120)
        rax.set_ylim(0, 2)
        rax.axis('off')
        rr = rt[rt.group == g]
        for row, (grp, col) in enumerate((('High', HIGH), ('Low', LOW))):
            q = rr[rr.risk_group == grp].sort_values('time')
            for t, n in zip(q.time, q.n_at_risk):
                if t in (0, 48, 96):
                    rax.text(t, 1.5 - row, f'{int(n)}', ha='center', va='center', fontsize=SZ, color=col)
        if c == 0:
            for row, (nm, col) in enumerate((('High', HIGH), ('Low', LOW))):
                cv.text(x - 3.4, y + ph + 4.2 + 1.45 + row * 2.9, nm, ha='right', va='center', color=col)
            cv.text(2.4, y + ph / 2, 'Overall survival', rotation=90, ha='center', va='center', fontsize=SZ_T)
    r, c = divmod(len(D_ORDER), ncol)
    x = left + c * cell_w
    y = top0 + r * cell_h
    for k, (nm, col) in enumerate((('CMI-high (High)', HIGH), ('CMI-low (Low)', LOW))):
        yy = y + ph * 0.35 + k * 3.8
        cv.hline(x + 1.5, x + 5.0, yy, color=col, lw=LW_DATA * 1.4, solid_capstyle='butt')
        cv.text(x + 6.0, yy, nm, ha='left', va='center', color=col)
    cv.text(left + 2 * cell_w - 5.5, top0 + 2 * cell_h + ph + 10.8, 'Months', ha='center', va='top', fontsize=SZ_T)
    ci = D['c20_r05c_time_dependent_AUC_95CI']
    y0 = top0 + 3 * cell_h + 1.0
    cv.letter(0, y0, 'B')
    cv.label(5.5, y0 + 0.6, f'{len(D_ORDER)} cohorts | n = {ntot:,} | {etot:,} deaths')
    n = len(D_ORDER)
    row_h = 5.6
    top = y0 + 12.0
    lab = cv.ax(0, top, 47.0, n * row_h)
    lab.axis('off')
    lab.set_xlim(0, 47)
    lab.set_ylim(n - 0.5, -0.5)
    ax = cv.ax(49.0, top, 52.0, n * row_h, row='S6B')
    num = cv.ax(103.0, top, 75.0, n * row_h)
    num.axis('off')
    num.set_xlim(0, 75)
    num.set_ylim(n - 0.5, -0.5)
    cols = {12: 24.0, 36: 49.5, 60: 75.0}
    off = {12: -0.27, 36: 0.0, 60: 0.27}
    prev = None
    for i, g in enumerate(D_ORDER):
        ent, lb = LABEL[g]
        lab.text(0.5, i, '' if ent == prev else ent, ha='left', va='center', fontsize=SZ, color=TXT)
        lab.text(46.5, i, f"{lb} ({int(st.loc[g, 'n'])})", ha='right', va='center', fontsize=SZ)
        if prev is not None and ent != prev:
            ax.axhline(i - 0.5, color=MID, lw=0.35, zorder=0)
        prev = ent
        q = ci[ci.group == g].set_index('eval_time')
        txts = []
        for t in (12, 36, 60):
            if t not in q.index:
                num.text(cols[t], i, '-', ha='right', va='center', fontsize=SZ, color=MID)
                txts.append('-')
                continue
            r = q.loc[t]
            if np.isfinite(r.lower):
                ax.plot([r.lower, r.upper], [i + off[t]] * 2, color=ROC3[t], lw=LW_ERR, zorder=2, solid_capstyle='butt')
                txt = f'{r.auc:.2f} ({r.lower:.2f}-{r.upper:.2f})'
            else:
                txt = f'{r.auc:.2f} (-)'
            ax.scatter([r.auc], [i + off[t]], s=9, c=ROC3[t], lw=0, zorder=3, clip_on=False)
            num.text(cols[t], i, txt, ha='right', va='center', fontsize=SZ)
            txts.append(txt)
        cv.tables.append({'panel': 'S6B', 'texts': [f"{lb} ({int(st.loc[g, 'n'])})"] + txts})
    for t, nm in ((12, '1-y'), (36, '3-y'), (60, '5-y')):
        num.text(cols[t], -0.95, nm, ha='right', va='center', fontsize=SZ)
        num.plot([cols[t] - 9.6, cols[t] - 6.6], [-0.95, -0.95], color=ROC3[t], lw=LW_DATA * 1.4, solid_capstyle='butt', clip_on=False)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0.65, 1.02)
    ax.set_xticks([0.7, 0.8, 0.9, 1.0])
    ax.set_xticklabels(['0.7', '0.8', '0.9', '1'])
    ax.tick_params(left=False, labelleft=False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('AUC')
    return (cv, overlaps)

def line_key(cv, x0, y, color, lw, label, marker=None, w=5.0):
    from matplotlib.lines import Line2D
    X = [(cv.ox + x0 + w * f) / cv.fw for f in (0, 0.5, 1)]
    Y = [1 - (cv.oy + y) / cv.fh] * 3
    cv.fig.add_artist(Line2D(X, Y, transform=cv.fig.transFigure, color=color, lw=lw, marker=marker or '', ms=2.0, markevery=[1], solid_capstyle='butt'))
    return cv.text(x0 + w + 1.2, y, label, va='center', color=TXT)

def page2(fig, fw, fh, ox, oy, D):
    cv = Canvas(fig, fw, fh, ox, oy)
    cv.tables = []
    a = D['c20_r07_published_model_integrated_AUC']
    a = a[a.group.isin(EF_ORDER)]
    fs = D['c20_r07c_published_model_integrated_AUC_cohort_summary'].set_index('group')
    if set(EF_ORDER) != set(fs.index[fs.show_in_figure]):
        sys.exit(f'Displayed cohorts do not match the input display flags: {sorted(set(EF_ORDER) ^ set(fs.index[fs.show_in_figure]))}')
    n_pub = int(a[~a.is_RSF].model.nunique())
    nn = D['c20_r08_time_dependent_AUC_long'].groupby('group').n_common.first().astype(int)
    cv.letter(0, 0, 'C')
    ev_e = int(D['c20_r03_KM_cohort_statistics'].set_index('group').loc[EF_ORDER, 'events'].sum())
    cv.label(5.5, 0.6, f'{len(EF_ORDER)} cohorts | n = {int(nn.loc[EF_ORDER].sum()):,} | {ev_e:,} deaths')
    n = len(EF_ORDER)
    row_h = 5.2
    top = 11.0
    lab = cv.ax(0, top, 44.0, n * row_h)
    lab.axis('off')
    lab.set_xlim(0, 44)
    lab.set_ylim(n - 0.5, -0.5)
    ax = cv.ax(46.0, top, 100.0, n * row_h, row='S6C')
    qa = cv.ax(148.0, top, 30.0, n * row_h)
    qa.axis('off')
    qa.set_xlim(0, 30)
    qa.set_ylim(n - 0.5, -0.5)
    yi = {g: i for i, g in enumerate(EF_ORDER)}
    lo, hi = (a.integrated_mean_auc.min(), a.integrated_mean_auc.max())
    if lo < 0.3 or hi > 1.0:
        sys.exit(f'Integrated AUC {lo:.3f}-{hi:.3f} exceeds the axis range')
    pub = a[~a.is_RSF]
    rsf = a[a.is_RSF]
    ax.scatter(pub.integrated_mean_auc, pub.group.map(yi), s=7, c=DOT, lw=0, zorder=2)
    ax.scatter(rsf.integrated_mean_auc, rsf.group.map(yi), s=20, marker='D', c=HIGH, edgecolors='black', linewidths=0.4, zorder=4)
    prev = None
    for i, g in enumerate(EF_ORDER):
        ent, lb = LABEL[g]
        lab.text(0.5, i, '' if ent == prev else ent, ha='left', va='center', fontsize=SZ, color=TXT)
        lab.text(43.5, i, f'{lb} ({nn[g]:,})', ha='right', va='center', fontsize=SZ)
        if prev is not None and ent != prev:
            ax.axhline(i - 0.5, color=MID, lw=0.35, zorder=1)
        prev = ent
        ax.axhline(i, color='#E6E7E9', lw=0.4, zorder=0)
        t = fmtP(fs.loc[g, 'q_best'])
        qa.text(29.5, i, t, ha='right', va='center', fontsize=SZ)
        cv.tables.append({'panel': 'S6C', 'texts': [f'{lb} ({nn[g]:,})', t]})
    qa.text(29.5, -1.05, 'q vs best', ha='right', va='center', fontsize=SZ, color=TXT)
    key(ax, 0.335, -1.05, 'D', HIGH, 'CMI-riskscore', size=16, edge='black', dx=0.012)
    key(ax, 0.62, -1.05, 'o', DOT, 'published', size=10, dx=0.012)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0.3, 1.0)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(['0.4', '0.6', '0.8', '1.0'])
    ax.tick_params(left=False, labelleft=False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Integrated AUC (12-60 months)')
    td = D['c20_r08_time_dependent_AUC_long']
    y0 = top + n * row_h + 14.0
    cv.letter(0, y0, 'D')
    td = td[td.group.isin(EF_ORDER)]
    pubs = td[~td.is_RSF].model.nunique()
    cv.label(5.5, y0 + 0.6, f'{len(EF_ORDER)} cohorts | n = {int(nn.loc[EF_ORDER].sum()):,} | {ev_e:,} deaths')
    ncol, pw, ph = (5, 29.5, 21.0)
    cw, ch = (33.6, 30.5)
    left, top0 = (12.0, y0 + 10.0)
    for i, g in enumerate(EF_ORDER):
        r, c = divmod(i, ncol)
        x = left + c * cw
        y = top0 + r * ch
        ent, lb = LABEL[g]
        q = td[td.group == g]
        cv.text(x, y - 3.6, (f'{lb}' if ent == 'Mix' else f'{ent} {lb}') + f' ({int(q.n_common.iloc[0]):,})')
        ax = cv.ax(x, y, pw, ph, row=f'S6D-{r + 1}')
        for m, z in q[~q.is_RSF].groupby('model'):
            z = z.sort_values('evaluation_time')
            ax.plot(z.evaluation_time, z.auc, color=DOT, lw=0.5, alpha=0.9, zorder=2)
        z = q[q.is_RSF].sort_values('evaluation_time')
        ax.plot(z.evaluation_time, z.auc, color=HIGH, lw=LW_DATA * 1.1, zorder=4, marker='o', ms=2.0)
        ax.axhline(0.5, color=MID, lw=0.4, ls=(0, (2, 2)), zorder=1)
        ax.set_xlim(8, 64)
        ax.set_ylim(0.2, 1.04)
        ax.set_xticks([12, 24, 36, 48, 60])
        ax.set_xticklabels(['12', '', '36', '', '60'])
        ax.set_yticks([0.2, 0.6, 1.0])
        ax.set_yticklabels(['0.2', '0.6', '1'])
        if c != 0:
            ax.tick_params(labelleft=False)
        if i + ncol < n:
            ax.tick_params(labelbottom=False)
        if c == 0:
            cv.text(2.4, y + ph / 2, 'AUC', rotation=90, ha='center', va='center', fontsize=SZ_T)
    cv.text(left + (ncol * cw - (cw - pw)) / 2, top0 + ch + ph + 4.4, 'Months', ha='center', va='top', fontsize=SZ_T)
    line_key(cv, 118.0, y0 + 2.0, HIGH, LW_DATA * 1.1, 'CMI-riskscore', marker='o')
    line_key(cv, 144.0, y0 + 2.0, DOT, 0.8, 'published')
    y0 = top0 + 2 * ch + 12.0
    cv.letter(0, y0, 'E')
    mod = D['c20.1_r02_multivariable_Cox_models'].set_index('disease')
    coef = D['c20.1_r01_multivariable_Cox_coefficients']
    sen = D['c20.1_r03_multivariable_Cox_sensitivity']
    sen = sen[sen.term.isin(['cmi', 'z'])]
    cv.label(5.5, y0 + 0.6, f'{len(mod)} lymphoma entities | n = {int(mod.n.sum()):,} | {int(mod.events.sum()):,} deaths')
    lines = []
    for dz in ('DLBCL', 'FL', 'MCL', 'BL'):
        m = coef[(coef.disease == dz) & (coef.term == 'cmi')].iloc[0]
        lines.append((dz, 'Main model', dict(n=int(mod.loc[dz, 'n']), ev=int(mod.loc[dz, 'events']), HR=m.HR, lo=m.lower, hi=m.upper, p=m.p)))
        for a_ in SENS_ORDER:
            z = sen[(sen.disease == dz) & (sen.analysis == a_)]
            if len(z):
                z = z.iloc[0]
                lines.append((dz, SENS[a_], dict(n=int(z.n), ev=int(z.events), HR=z.HR, lo=z.lower, hi=z.upper, p=z.p)))
    used = set(sen.analysis)
    assert used <= set(SENS_ORDER), f'Unrecognized sensitivity analysis: {used - set(SENS_ORDER)}'
    n = len(lines)
    top = y0 + 9.5
    hh = n * 3.55
    if top + hh + 11.0 > H:
        sys.exit(f'Panel requires {top + hh + 11.0:.1f} mm; canvas height is {H} mm')
    lab = cv.ax(0, top, 56.0, hh)
    lab.axis('off')
    lab.set_xlim(0, 56)
    lab.set_ylim(n - 0.5, -0.5)
    ne = cv.ax(57.0, top, 20.0, hh)
    ne.axis('off')
    ne.set_xlim(0, 20)
    ne.set_ylim(n - 0.5, -0.5)
    ax = cv.ax(79.0, top, 42.0, hh, row='S6E')
    ax.patch.set_alpha(0)
    num = cv.ax(123.0, top, 55.0, hh)
    num.axis('off')
    num.set_xlim(0, 55)
    num.set_ylim(n - 0.5, -0.5)
    ne.text(19.5, -1.05, 'n / events', ha='right', va='center', fontsize=SZ, color=TXT)
    num.text(36.0, -1.05, 'HR (95% CI)', ha='right', va='center', fontsize=SZ, color=TXT)
    num.text(54.5, -1.05, 'P', ha='right', va='center', fontsize=SZ, color=TXT)
    lo_x, hi_x = (1.0, 500.0)
    prev = None
    for i, (dz, nm, v) in enumerate(lines):
        if dz != prev:
            lab.text(0.3, i, dz, ha='left', va='center', fontsize=SZ, color=TXT)
            if prev is not None:
                ax.axhline(i - 0.5, color=MID, lw=0.3, zorder=0)
            prev = dz
        col = HIGH if nm == 'Main model' else 'black'
        lab.text(55.5, i, nm, ha='right', va='center', fontsize=SZ, color=col)
        ne.text(19.5, i, f"{v['n']:,} / {v['ev']:,}", ha='right', va='center', fontsize=SZ)
        a_, b_ = (max(v['lo'], lo_x), min(v['hi'], hi_x))
        ax.plot([a_, b_], [i, i], color=HIGH if nm == 'Main model' else '#55595E', lw=LW_ERR, zorder=2, solid_capstyle='butt')
        if v['hi'] > hi_x:
            ax.annotate('', xy=(hi_x, i), xytext=(hi_x * 0.6, i), arrowprops=dict(arrowstyle='-|>', lw=LW_ERR, color=col, mutation_scale=5), zorder=2)
        ax.scatter([v['HR']], [i], s=16 if nm == 'Main model' else 11, marker='s', c=HIGH if nm == 'Main model' else '#55595E', lw=0, zorder=3)
        t1 = f"{fmt_ci(v['HR'])} ({fmt_ci(v['lo'])}-{fmt_ci(v['hi'])})"
        t2 = fmtP(v['p'])
        num.text(36.0, i, t1, ha='right', va='center', fontsize=SZ, color=col)
        num.text(54.5, i, t2, ha='right', va='center', fontsize=SZ, color=col)
        cv.tables.append({'panel': 'S6E', 'texts': [nm, t1, t2]})
    ax.axvline(1, color=MID, lw=0.45, ls=(0, (2, 2)), zorder=0)
    ax.set_xscale('log')
    ax.set_xlim(0.8, hi_x)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks([1, 10, 100])
    ax.set_xticklabels(['1', '10', '100'])
    ax.minorticks_off()
    ax.tick_params(left=False, labelleft=False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('Hazard ratio')
    return cv
SENS = {'Continuous CMI-riskscore (per 1 SD)': 'Continuous score (per 1 SD)', 'DLBCL all treatment intents': 'All treatment intents', 'DLBCL pooled cohorts (stratified by cohort)': 'Pooled 4 cohorts (stratified)', 'FL failure-free survival endpoint': 'Failure-free survival', 'MCL adjusted for MIPI only': 'Adjusted for MIPI only', 'MCL diagnostic biopsies only': 'Diagnostic biopsies only', 'BL additionally adjusted for LDH > ULN (Firth)': 'Plus LDH > ULN', 'Follow-up interval: 0-12 months': 'Follow-up 0-12 months', 'Follow-up interval: 12-36 months': 'Follow-up 12-36 months', 'Follow-up interval: >36 months': 'Follow-up > 36 months', 'FL follow-up interval: 0-36 months (Firth)': 'Follow-up 0-36 months', 'FL follow-up interval: >36 months (landmark)': 'Follow-up > 36 months'}
SENS_ORDER = ['Continuous CMI-riskscore (per 1 SD)', 'DLBCL all treatment intents', 'DLBCL pooled cohorts (stratified by cohort)', 'FL failure-free survival endpoint', 'MCL adjusted for MIPI only', 'MCL diagnostic biopsies only', 'BL additionally adjusted for LDH > ULN (Firth)', 'Follow-up interval: 0-12 months', 'Follow-up interval: 12-36 months', 'Follow-up interval: >36 months', 'FL follow-up interval: 0-36 months (Firth)', 'FL follow-up interval: >36 months (landmark)']

def page3(fig, fw, fh, ox, oy, D):
    cv = Canvas(fig, fw, fh, ox, oy)
    cv.tables = []
    rows = D['c20_r10_clinical_associations_group_statistics']
    pts = D['c20_r11_clinical_associations_subject_data']
    gm = D['c20_r12_clinical_analysis_sample_counts'].set_index('disease')
    cv.letter(0, 0, 'F')
    ents = ['DLBCL', 'BL', 'FL', 'MCL', 'MZL']
    cv.label(5.5, 0.6, f'{rows.row_order.nunique()} clinical / molecular features | ' + ' | '.join((f"{k} {int(gm.loc[k, 'n_scored'] - gm.loc[k, 'n_QC_excluded']):,}" for k in ents)) + ' scored tumours')
    keys = list(dict.fromkeys(rows.sort_values('row_order').row_order))
    ncol, bw, bh = (5, 28.0, 31.0)
    cw, ch = (33.2, 55.0)
    slot = 0
    for i, k in enumerate(keys):
        rr = rows[rows.row_order == k].sort_values('order')
        span = 2 if len(rr) > 3 else 1
        if slot % ncol + span > ncol:
            slot += ncol - slot % ncol
        r, c = divmod(slot, ncol)
        slot += span
        x = 12.0 + c * cw
        y = 14.0 + r * ch
        cv.text(x, y - 7.0, rr.disease.iloc[0], color=TXT)
        cv.text(x, y - 3.6, rr.row_label.iloc[0])
        ax = cv.ax(x, y, bw + (span - 1) * cw, bh, row=f'S6F-{r + 1}')
        lv = list(rr.level_label)
        cols = [LOW, HIGH] if len(lv) == 2 else [LOW] + [DOT] * (len(lv) - 2) + [HIGH]
        data = [pts[(pts.row_order == k) & (pts.level_label == L)].score.values for L in lv]
        pos = violin_box(ax, data, cols, seed=7 + i, pt_size=1.1)
        top = max((np.max(v) for v in data))
        p_txt = f'P = {fmtP(float(rr.p_show.iloc[0]))}'
        p_bracket(ax, pos[0], pos[-1], top + 0.45, p_txt)
        ax.set_xlim(0.4, len(lv) + 0.6)
        ax.set_ylim(-3.3, 8.2)
        ax.set_yticks([-3, 0, 3, 6])
        ax.set_xticks(pos)
        ax.set_xticklabels([f'{L}\n{int(n):,}' for L, n in zip(lv, rr.n)], linespacing=1.15)
        ax.tick_params(axis='x', length=1.2)
        if c != 0:
            ax.tick_params(labelleft=False)
        if c == 0:
            cv.text(2.4, y + bh / 2, 'CMI-riskscore (z)', rotation=90, ha='center', va='center', fontsize=SZ_T)
    return cv

def next_version(d: Path, name: str) -> int:
    pat = re.compile(f'^{re.escape(name)}_v(\\d+)\\.(pdf|svg|png)$')
    vs = [int(m.group(1)) for f in d.iterdir() if (m := pat.match(f.name))]
    return max(vs) + 1 if vs else 1

def check_alignment(cvs):
    import collections
    rows = collections.defaultdict(list)
    for cv in cvs:
        for r, x, yb in cv.reg:
            rows[r].append((x, yb))
    out = {'row_baseline_deviation_mm': {r: round(max((b for _, b in v)) - min((b for _, b in v)), 4) for r, v in rows.items()}}
    lefts = {r: min((x for x, _ in v)) for r, v in rows.items() if r.startswith(('S6A', 'S6D', 'S6F'))}
    out['left_edge_deviation_mm'] = round(max(lefts.values()) - min(lefts.values()), 4)
    bad = {k: v for k, v in out['row_baseline_deviation_mm'].items() if v > 1e-06}
    if bad or out['left_edge_deviation_mm'] > 1e-06:
        sys.exit(f"Alignment validation failed: {bad} / left edge {out['left_edge_deviation_mm']}")
    return out
PAGES = ((1, page1), (2, page2), (3, page3))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', type=int, default=None)
    ap.add_argument('--draft_dir', default=None, help='Save draft renders in this directory.')
    a = ap.parse_args()
    fam = setup_fonts()
    D = load()
    MM = 1 / 25.4
    if a.draft_dir:
        out = Path(a.draft_dir)
        out.mkdir(parents=True, exist_ok=True)
        for k, fn in PAGES:
            fig = plt.figure(figsize=(W * MM, H * MM))
            res = fn(fig, W, H, 0, 0, D)
            if k == 1:
                print('Overlapping panels: ', res[1])
            fig.savefig(out / f'S6_p{k}.png', dpi=300)
            fig.savefig(out / f'S6_p{k}.pdf')
            plt.close(fig)
        return
    OUTDIR.mkdir(parents=True, exist_ok=True)
    v = a.version or next_version(OUTDIR, figure_name)
    base = OUTDIR / f'{figure_name}_v{v}'
    cvs = []
    overlaps = []
    highlights = []
    table_rows = []
    with PdfPages(f'{base}.pdf') as pdf:
        for k, fn in PAGES:
            fig = plt.figure(figsize=(W * MM, H * MM))
            res = fn(fig, W, H, 0, 0, D)
            if k == 1:
                cv, overlaps = res
            else:
                cv = res
            cvs.append(cv)
            highlights.append(highlight_artist_count(fig, HIGH))
            table_rows += [dict(t, page=k) for t in cv.tables]
            pdf.savefig(fig)
            fig.savefig(OUTDIR / f'{figure_name}_page_{k}_v{v}.svg')
            fig.savefig(OUTDIR / f'{figure_name}_page_{k}_v{v}.png', dpi=600)
            plt.close(fig)
    alignment = check_alignment(cvs)
    if overlaps:
        sys.exit(f'No non-overlapping HR and P label position in panels: {overlaps}')
    shutil.copy2(OUTDIR / f'{figure_name}_page_1_v{v}.png', f'{base}.png')
    with PdfPages(OUTDIR / f'{figure_name}_A4_v{v}.pdf') as pdf:
        for k, fn in PAGES:
            fig = plt.figure(figsize=(A4W * MM, A4H * MM))
            fn(fig, A4W, A4H, A4_L, A4_T, D)
            pdf.savefig(fig)
            fig.savefig(OUTDIR / f'{figure_name}_A4_page_{k}_v{v}.png', dpi=300)
            plt.close(fig)
    rec = {'version': f'v{v}', 'figure_name': figure_name, 'font': fam, 'canvas_mm': [W, H], 'n_pages': len(PAGES), 'font_size_pt': {'ticks': SZ, 'axis_titles': SZ_T, 'panel_letter': SZ_TAG}, 'highlight_color': HIGH, 'highlight_used': all((k > 0 for k in highlights)), 'highlight_artist_counts_by_page': highlights, 'check_alignment': alignment, 'HR_P_overlap_panels': overlaps, 'table_rows': table_rows}
    (OUTDIR / f'{figure_name}_v{v}_validation.json').write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'[OK] {os.path.basename(str(base))}.pdf({len(PAGES)} pages); per-page SVG and PNG; A4 {figure_name}_A4_v{v}.pdf ;font {fam}; alignment {alignment}')
if __name__ == '__main__':
    main()
