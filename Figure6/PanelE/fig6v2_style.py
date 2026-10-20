import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
W = 180.0
H = float(os.environ.get('FIG6_H', '256.0'))
A4W, A4H, A4_L, A4_T = (210.0, 297.0, 15.0, 20.0)
SZ = float(os.environ.get('FIG6_SZ', '8.0'))
SZ_T = SZ * 1.1
SZ_TAG = SZ * 1.4
LW_DATA, LW_ERR, LW_AXIS = (1.0, 0.65, 0.55)
HIGH, LOW = ('#D8503C', '#3C7DB1')
DOT = '#B9BDC3'
MID = '#8A9197'
TXT = '#5F656B'
BAND = '#F1F2F4'
ROC3 = {12: '#EFA394', 36: '#D8503C', 60: '#8C2B1F'}

def fmtP(p):
    return '' if p is None or (isinstance(p, float) and np.isnan(p)) else f'{p:.2e}'

def setup_fonts():
    fams = [f.name for f in font_manager.fontManager.ttflist]
    fam = 'Arial' if 'Arial' in fams else 'Liberation Sans' if 'Liberation Sans' in fams else None
    if fam is None:
        sys.exit('Arial or Liberation Sans is required.')
    plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'pdf.use14corefonts': False, 'svg.fonttype': 'none', 'font.family': fam, 'font.size': SZ, 'axes.labelsize': SZ_T, 'xtick.labelsize': SZ, 'ytick.labelsize': SZ, 'axes.linewidth': LW_AXIS, 'xtick.major.width': LW_AXIS, 'ytick.major.width': LW_AXIS, 'xtick.major.size': SZ * 0.25, 'ytick.major.size': SZ * 0.25, 'xtick.major.pad': SZ * 0.3, 'ytick.major.pad': SZ * 0.3, 'axes.labelpad': SZ * 0.35, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False, 'legend.frameon': False, 'axes.unicode_minus': True, 'lines.solid_capstyle': 'butt'})
    return fam

class Canvas:

    def __init__(self, fig, fw, fh, ox, oy):
        self.fig, self.fw, self.fh, self.ox, self.oy = (fig, fw, fh, ox, oy)
        self.reg = []

    def ax(self, x, y, w, h, row=None, **kw):
        if row is not None:
            self.reg.append((row, x, y + h))
        return self.fig.add_axes([(self.ox + x) / self.fw, 1 - (self.oy + y + h) / self.fh, w / self.fw, h / self.fh], **kw)

    def text(self, x, y, s, **kw):
        kw.setdefault('fontsize', SZ)
        kw.setdefault('va', 'top')
        kw.setdefault('ha', 'left')
        return self.fig.text((self.ox + x) / self.fw, 1 - (self.oy + y) / self.fh, s, **kw)

    def letter(self, x, y, s):
        dy = max(0.6, SZ_TAG * 0.079)
        return self.text(x + 0.6, y + dy, s, fontsize=SZ_TAG, fontweight='bold', va='top', ha='left')

    def label(self, x, y, s, **kw):
        kw.setdefault('color', TXT)
        return self.text(x, y, s, **kw)

    def hline(self, x0, x1, y, **kw):
        from matplotlib.lines import Line2D
        X = [(self.ox + x0) / self.fw, (self.ox + x1) / self.fw]
        Y = [1 - (self.oy + y) / self.fh] * 2
        return self.fig.add_artist(Line2D(X, Y, transform=self.fig.transFigure, **kw))

def clean(ax, left=True, bottom=True):
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)
    if not left:
        ax.tick_params(left=False, labelleft=False)
    if not bottom:
        ax.tick_params(bottom=False, labelbottom=False)
MIN_RISK = 3

def km_steps(t, e, min_risk=1):
    t = np.asarray(t, float)
    e = np.asarray(e, int)
    order = np.lexsort((-e, t))
    t, e = (t[order], e[order])
    uniq = np.unique(t)
    n = len(t)
    s = 1.0
    xs = [0.0]
    ys = [1.0]
    cx, cy = ([], [])
    at = n
    last = 0.0
    for u in uniq:
        if at < min_risk:
            break
        m = t == u
        d = int(e[m].sum())
        c = int((1 - e[m]).sum())
        if d > 0:
            xs += [u, u]
            ys += [s, s * (1 - d / at)]
            s = s * (1 - d / at)
        if c > 0:
            cx.append(u)
            cy.append(s)
        at -= d + c
        last = u
    xs.append(last)
    ys.append(s)
    return (np.array(xs), np.array(ys), np.array(cx), np.array(cy))

def text_band(curves, x_from, x_to, h, cands, pad=0.03):
    best, best_hit = (cands[0], 1000000000.0)
    grid = np.linspace(x_from, x_to, 160)
    for yb in cands:
        hit = 0
        for X, Y in curves:
            idx = np.searchsorted(X, grid, side='right') - 1
            ok = (idx >= 0) & (grid <= X.max())
            yy = Y[np.clip(idx, 0, len(Y) - 1)][ok]
            hit += int(((yy >= yb - pad) & (yy <= yb + h + pad)).sum())
        if hit == 0:
            return yb
        if hit < best_hit:
            best, best_hit = (yb, hit)
    return best

def highlight_artist_count(fig, color):
    from matplotlib.colors import to_rgba
    target = np.array(to_rgba(color))[:3]

    def same(c):
        try:
            arr = np.atleast_2d(np.asarray(to_rgba(c) if isinstance(c, str) else c, float))
        except (ValueError, TypeError):
            return False
        return arr.shape[-1] >= 3 and bool(np.any(np.all(np.abs(arr[:, :3] - target) < 0.001, axis=1)))
    n = 0
    for a in fig.findobj():
        for getter in ('get_color', 'get_facecolor', 'get_edgecolor', 'get_markerfacecolor'):
            if hasattr(a, getter):
                try:
                    if same(getattr(a, getter)()):
                        n += 1
                        break
                except Exception:
                    pass
    return n

def _segments_hit(curves, x0, x1, y0, y1):
    for X, Y in curves:
        X = np.asarray(X, float)
        Y = np.asarray(Y, float)
        xa, xb, ya, yb_ = (X[:-1], X[1:], Y[:-1], Y[1:])
        lo_x, hi_x = (np.minimum(xa, xb), np.maximum(xa, xb))
        lo_y, hi_y = (np.minimum(ya, yb_), np.maximum(ya, yb_))
        if np.any((hi_x >= x0) & (lo_x <= x1) & (hi_y >= y0) & (lo_y <= y1)):
            return True
    return False

def place_lines(ax, variants, curves, x_right, x_left, y_cands, line_mm=3.25, pad_mm=0.7, **kw):
    fig = ax.figure
    r = fig.canvas.get_renderer()
    bb = ax.get_window_extent(r)
    (xa, xb), (ya, yb) = (ax.get_xlim(), ax.get_ylim())
    px_per_mm = fig.dpi / 25.4
    sx = (xb - xa) / bb.width
    sy = (yb - ya) / bb.height
    lh = line_mm * px_per_mm * sy
    padx = pad_mm * px_per_mm * sx
    pady = pad_mm * px_per_mm * sy

    def width(s):
        t = ax.text(0, 0, s, fontsize=kw.get('fontsize', SZ))
        w = t.get_window_extent(r).width * sx
        t.remove()
        return w

    def clear(ws, ha, xx, yb0):
        n = len(ws)
        for k, w in enumerate(ws):
            y0 = yb0 + (n - 1 - k) * lh
            x0, x1 = (xx - w, xx) if ha == 'right' else (xx, xx + w)
            if x0 < xa or x1 > xb or y0 < ya or (y0 + lh > yb):
                return False
            if _segments_hit(curves, x0 - padx, x1 + padx, y0 - pady, y0 + lh + pady):
                return False
        return True

    def draw(lines, ha, xx, yb0, **extra):
        for k, s in enumerate(lines):
            ax.text(xx, yb0 + (len(lines) - 1 - k) * lh, s, ha=ha, va='bottom', zorder=6, **extra, **kw)
    best = None
    for lines in variants:
        ws = [width(s) for s in lines]
        for yb0 in y_cands:
            for ha, xx in (('right', x_right), ('left', x_left)):
                if clear(ws, ha, xx, yb0):
                    draw(lines, ha, xx, yb0)
                    return True
                if best is None:
                    best = (lines, ha, xx, yb0)
        free = []
        hh = lh * len(lines)
        for yb0 in np.arange(ya + 0.02 * (yb - ya), yb - hh, 0.01 * (yb - ya)):
            for x1 in np.arange(x_right, xa + max(ws), -0.01 * (xb - xa)):
                if clear(ws, 'right', x1, yb0):
                    free.append((x1, yb0))
                    break
        if free:
            yc_pref = ya + 0.5 * (yb - ya)
            x1, yb0 = max(free, key=lambda t: (round(t[0], 6), -abs(t[1] + hh / 2 - yc_pref)))
            draw(lines, 'right', x1, yb0)
            return True
    lines, ha, xx, yb0 = best
    draw(lines, ha, xx, yb0, bbox=dict(boxstyle='square,pad=0.05', fc='white', ec='none'))
    return False

def text_mm(s, size=SZ):
    w = 0.0
    for ch in s:
        if ch in 'MWmw':
            w += 0.83
        elif ch.isupper():
            w += 0.69
        elif ch in " .,|:;|()[]---'":
            w += 0.3
        else:
            w += 0.556
    return w * size * 0.3528

def key(ax, x, y, marker, color, label, size=12, edge=None, dx=1.6, trans='data'):
    ax.scatter([x], [y], s=size, marker=marker, c=color, edgecolors=edge or 'none', linewidths=0.4 if edge else 0, zorder=6, clip_on=False)
    return ax.text(x + dx, y, label, ha='left', va='center', fontsize=SZ, color=TXT, clip_on=False)

def violin_box(ax, data, cols, pos=None, jitter_max_n=300, seed=20261011, pt_size=1.4):
    pos = np.arange(1, len(data) + 1) if pos is None else np.asarray(pos, float)
    vp = ax.violinplot(data, positions=pos, widths=0.78, showextrema=False, bw_method='scott')
    for body, col in zip(vp['bodies'], cols):
        body.set_facecolor(col)
        body.set_edgecolor('none')
        body.set_alpha(0.22)
        body.set_zorder(1)
    bp = ax.boxplot(data, positions=pos, widths=0.26, showfliers=False, patch_artist=True, whis=1.5, medianprops=dict(color='white', lw=0.9), whiskerprops=dict(lw=0.55), capprops=dict(lw=0))
    for patch, col in zip(bp['boxes'], cols):
        patch.set_facecolor(col)
        patch.set_edgecolor(col)
        patch.set_linewidth(0.5)
        patch.set_zorder(3)
    for med in bp['medians']:
        med.set_zorder(4.5)
        med.set_solid_capstyle('butt')
    for j, col in enumerate(cols):
        for wsk in bp['whiskers'][2 * j:2 * j + 2]:
            wsk.set_color(col)
            wsk.set_zorder(2)
    rng = np.random.default_rng(seed)
    for j, (v, col) in enumerate(zip(data, cols)):
        if len(v) <= jitter_max_n:
            ax.scatter(pos[j] + rng.uniform(-0.17, 0.17, len(v)), v, s=pt_size, c=col, lw=0, alpha=0.75, zorder=4)
    return pos

def p_bracket(ax, x0, x1, y, text, h=0.25):
    ax.plot([x0, x0, x1, x1], [y - h, y, y, y - h], color='black', lw=0.5, zorder=5)
    return ax.text((x0 + x1) / 2, y + 0.15, text, ha='center', va='bottom', fontsize=SZ, zorder=5)
