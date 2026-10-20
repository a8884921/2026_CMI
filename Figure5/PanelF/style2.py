import os, re, json, hashlib, shutil, datetime, pathlib
import numpy as np
import logging
import matplotlib as mpl
mpl.use('Agg')
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.spatial import cKDTree
MM = 1 / 25.4
H_PLOT = 33.3
SZ_TICK = float(np.clip(round(H_PLOT / 3.7, 2), 8, 12))
SZ_SMALL = round(SZ_TICK * 0.9, 2)
SZ_TITLE = round(SZ_TICK * 1.1, 2)
SZ_TAG = round(SZ_TICK * 1.4, 2)
assert SZ_SMALL >= 8.0, 'self-arranged A4 figures must stay >= 8 pt'
PT = 25.4 / 72
LW_DATA, LW_ERR, LW_AXIS = (1.1, 0.65, 0.55)
C = {'Mal_1': '#3C5488', 'Mal_2': '#F39B7F', 'Mal_3': '#E64B35', 'Mal_3_light': '#F6C1B8', 'Mal_1_light': '#B9C3DC', 'Mal_2_light': '#FAD3C6', 'CD8T': '#00A087', 'CD8T_light': '#9FD8CB', 'Mac': '#7E6148', 'CD4T': '#8491B4', 'ink': '#1A1A1A', 'dark': '#4D4D4D', 'neutral': '#8A9197', 'light': '#C4C8CB', 'faint': '#E6E7E8', 'spot_bg': '#E3E4E5', 'zero': '#9AA0A4', 'ns': '#8491B4'}
STATES = ['Mal_1', 'Mal_2', 'Mal_3']
DISEASES = ['DLBCL', 'PCNSL', 'MZL', 'BL']
MARK = {'DLBCL': 'o', 'PCNSL': 's', 'MZL': '^', 'BL': 'D'}
MARK_SCALE = {'o': 1.0, 's': 0.85, '^': 1.05, 'D': 0.75}
CMAP_DIV = LinearSegmentedColormap.from_list('div', ['#3C5488', '#8491B4', '#F7F7F7', '#F39B7F', '#E64B35'])
CMAP_STATE = LinearSegmentedColormap.from_list('cmi', ['#3C5488', '#D9D9D9', '#E64B35'])

def apply():
    S = SZ_TICK
    mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans'], 'font.size': S, 'axes.labelsize': SZ_TITLE, 'axes.titlesize': SZ_TITLE, 'xtick.labelsize': S, 'ytick.labelsize': S, 'legend.fontsize': SZ_SMALL, 'axes.labelpad': S * 0.45, 'xtick.major.pad': S * 0.3, 'ytick.major.pad': S * 0.3, 'xtick.major.size': S * 0.25, 'ytick.major.size': S * 0.25, 'xtick.minor.size': S * 0.15, 'ytick.minor.size': S * 0.15, 'axes.linewidth': LW_AXIS, 'xtick.major.width': LW_AXIS, 'ytick.major.width': LW_AXIS, 'xtick.minor.width': LW_AXIS, 'ytick.minor.width': LW_AXIS, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': False, 'legend.frameon': False, 'axes.edgecolor': C['ink'], 'xtick.color': C['ink'], 'ytick.color': C['ink'], 'text.color': C['ink'], 'axes.labelcolor': C['ink'], 'lines.linewidth': LW_DATA, 'lines.solid_capstyle': 'butt', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'pdf.use14corefonts': False, 'svg.fonttype': 'none', 'savefig.dpi': 600, 'savefig.bbox': 'standard', 'savefig.pad_inches': 0, 'axes.unicode_minus': True, 'mathtext.default': 'regular'})

def fmtP(p):
    return '' if p is None or p != p else f'{p:.2e}'

def Ptxt(p):
    return f'P = {fmtP(p)}'

def stars(p):
    if p is None or p != p:
        return ''
    return '****' if p < 0.0001 else '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''

def star_colour(rgb):
    L = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return 'black' if L > 0.55 else 'white'

def minus(s):
    return str(s).replace('-', '-')

def nfmt(n):
    return f'{int(n):,}'

class Page:

    def __init__(self, w, h, off=(0.0, 0.0)):
        self.w, self.h = (w, h)
        self.dx, self.dy = off
        self.fig = plt.figure(figsize=(w * MM, h * MM))
        self.ink_boxes = []

    def X(self, x):
        return (x + self.dx) / self.w

    def Y(self, y):
        return 1 - (y + self.dy) / self.h

    def ax(self, x, y, w, h, **kw):
        return self.fig.add_axes([self.X(x), self.Y(y + h), w / self.w, h / self.h], **kw)

    def text(self, x, y, s, size=None, **kw):
        kw.setdefault('va', 'top')
        kw.setdefault('ha', 'left')
        kw.setdefault('fontsize', size or SZ_TICK)
        return self.fig.text(self.X(x), self.Y(y), s, **kw)

    def letter(self, x, y, L):
        return self.text(x, y, L, size=SZ_TAG, fontweight='bold', va='top')

    def tag(self, x, y, s, **kw):
        kw.setdefault('color', C['dark'])
        return self.text(x, y, s, size=SZ_SMALL, va='top', **kw)

    def line(self, xs, ys, **kw):
        kw.setdefault('color', C['ink'])
        kw.setdefault('lw', LW_AXIS)
        self.fig.add_artist(mpl.lines.Line2D([self.X(x) for x in xs], [self.Y(y) for y in ys], transform=self.fig.transFigure, **kw))

    def rect(self, x, y, w, h, **kw):
        kw.setdefault('lw', 0)
        self.fig.add_artist(Rectangle((self.X(x), self.Y(y + h)), w / self.w, h / self.h, transform=self.fig.transFigure, **kw))

    def dot(self, x, y, r, **kw):
        kw.setdefault('lw', 0)
        self.fig.add_artist(mpl.patches.Ellipse((self.X(x), self.Y(y)), 2 * r / self.w, 2 * r / self.h, transform=self.fig.transFigure, **kw))

def text_w(s, size, bold=False):
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    fp = FontProperties(family=['Arial', 'Liberation Sans'], weight='bold' if bold else 'normal', size=size)
    return TextPath((0, 0), s, prop=fp).get_extents().width * PT

def clean(ax, left=True, bottom=True):
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)
    if not left:
        ax.set_yticks([])
    if not bottom:
        ax.set_xticks([])

def blank(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

def state_ticks(ax, ns=None, axis='x', small_n=True):
    labs = []
    for i, s in enumerate(STATES):
        labs.append(s if ns is None else f'{s}\n')
    if axis == 'x':
        ax.set_xticks(range(3))
        ax.set_xticklabels(labs)
        tl = ax.get_xticklabels()
    else:
        ax.set_yticks(range(3))
        ax.set_yticklabels(labs)
        tl = ax.get_yticklabels()
    for t, s in zip(tl, STATES):
        t.set_color(C[s])
    if ns is not None and axis == 'x':
        import matplotlib.transforms as mtrans
        tr = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
        off = mtrans.ScaledTranslation(0, -(SZ_TICK * 0.3 + SZ_TICK * 0.25 + SZ_TICK * 1.12) / 72, ax.figure.dpi_scale_trans)
        for i, n in enumerate(ns):
            ax.text(i, 0, f'n = {n}', transform=tr + off, ha='center', va='top', fontsize=SZ_SMALL, color=C['dark'], clip_on=False)

def bracket(ax, x0, x1, y, h, txt, color=None, size=None, transform='axes_y'):
    import matplotlib.transforms as mtrans
    color = color or C['ink']
    tr = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([x0, x0, x1, x1], [y - h, y, y, y - h], color=color, lw=LW_AXIS, transform=tr, clip_on=False, solid_capstyle='butt')
    ax.text((x0 + x1) / 2, y + 0.012, txt, ha='center', va='bottom', fontsize=size or SZ_SMALL, color=color, transform=tr, clip_on=False)

def omnibus(ax, x0, x1, y, txt, color=None, size=None):
    import matplotlib.transforms as mtrans
    color = color or C['ink']
    tr = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([x0, x1], [y, y], color=color, lw=LW_AXIS, transform=tr, clip_on=False, solid_capstyle='butt')
    ax.text((x0 + x1) / 2, y + 0.012, txt, ha='center', va='bottom', fontsize=size or SZ_SMALL, color=color, transform=tr, clip_on=False)

def group_p(ax, ps, xs=None, levels=None, y0=1.04, color=None, size=None, ha=None):
    import matplotlib.transforms as mtrans
    xs = list(range(len(ps))) if xs is None else xs
    levels = [0] * len(ps) if levels is None else levels
    ha = ['center'] * len(ps) if ha is None else ha
    size = size or SZ_SMALL
    h_mm = ax.get_position().height * ax.figure.get_size_inches()[1] * 25.4
    step = size * 1.22 * PT / h_mm
    tr = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
    for x, p, lv, a in zip(xs, ps, levels, ha):
        ax.text(x, y0 + lv * step, Ptxt(p), ha=a, va='bottom', fontsize=size, color=color or C['ink'], transform=tr, clip_on=False)

def zero_line(ax, y=0, **kw):
    kw.setdefault('color', C['zero'])
    kw.setdefault('lw', LW_AXIS)
    kw.setdefault('ls', (0, (2.5, 2)))
    kw.setdefault('zorder', 0)
    ax.axhline(y, **kw)

def hd_map(ax, case, kind, meta, maps, edge=None):
    q = maps[maps['map'] == f'{case}_{kind}']
    h = q.side.values / 2
    verts = np.stack([np.c_[q.u - h, q.v - h], np.c_[q.u + h, q.v - h], np.c_[q.u + h, q.v + h], np.c_[q.u - h, q.v + h]], axis=1)
    ax.add_collection(PolyCollection(verts, facecolors=q.color.values, edgecolors=q.color.values, linewidths=0.05))
    segs = np.array(meta['maps'][f'{case}_zone']['segments']).reshape(-1, 2, 2)
    ax.add_collection(LineCollection(segs, colors=edge or ('#5F6B70' if kind == 'zone' else 'white'), linewidths=LW_ERR))
    ax.set_xlim(0, 1)
    ax.set_ylim(1, 0)
    blank(ax)
TINT = {'#d6e9f5': 'FL', '#e4e4e4': 'transition', '#fbe3d6': 'tFL'}

def zone_anchor(maps, case):
    z = maps[maps['map'] == f'{case}_zone'].copy()
    z['zone'] = z.color.map(TINT)
    step = float(z.side.median())
    gi = np.round(z.u / step).astype(int)
    gj = np.round(z.v / step).astype(int)
    occ = set(zip(gi, gj))
    I, J = np.meshgrid(np.arange(gi.min() - 2, gi.max() + 3), np.arange(gj.min() - 2, gj.max() + 3))
    empty = np.array([(i, j) for i, j in zip(I.ravel(), J.ravel()) if (i, j) not in occ], float) * step
    out = {}
    for name in ['FL', 'tFL']:
        a = z[z.zone == name][['u', 'v']].values
        b = np.vstack([z[z.zone != name][['u', 'v']].values, empty])
        dd, _ = cKDTree(b).query(a)
        out[name] = a[np.argmax(dd)]
    return out

def vector_cbar(ax, cmap, vmin, vmax, ticks, ticklabels=None, n=128, orientation='horizontal'):
    cm = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    edges = np.linspace(vmin, vmax, n + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        col = cm((0.5 * (a + b) - vmin) / (vmax - vmin))
        if orientation == 'horizontal':
            ax.add_patch(Rectangle((a, 0), (b - a) * 1.04, 1, fc=col, ec='none', lw=0))
        else:
            ax.add_patch(Rectangle((0, a), 1, (b - a) * 1.04, fc=col, ec='none', lw=0))
    labs = ticklabels if ticklabels is not None else [minus(f'{t:g}') for t in ticks]
    if orientation == 'horizontal':
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks(ticks)
        ax.set_xticklabels(labs)
        for s in ['top', 'right', 'left']:
            ax.spines[s].set_visible(False)
        ax.tick_params(axis='x', labelsize=SZ_SMALL, length=SZ_TICK * 0.2, pad=SZ_TICK * 0.2)
    else:
        ax.set_ylim(vmin, vmax)
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.yaxis.tick_right()
        ax.set_yticks(ticks)
        ax.set_yticklabels(labs)
        for s in ['top', 'left', 'bottom']:
            ax.spines[s].set_visible(False)
        ax.spines['right'].set_visible(True)
        ax.tick_params(axis='y', labelsize=SZ_SMALL, length=SZ_TICK * 0.2, pad=SZ_TICK * 0.2)

def heat_cells(ax, M, cmap, vmin, vmax, P=None, fmt=None, center=None, grid='white', gw=0.9):
    norm = Normalize(vmin, vmax)
    nr, nc = M.shape
    for i in range(nr):
        for j in range(nc):
            v = float(M[i, j])
            rgba = cmap(norm(np.clip(v, vmin, vmax)))
            ax.add_patch(Rectangle((j, nr - 1 - i), 1, 1, fc=rgba, ec=grid, lw=gw))
            txt = []
            if fmt:
                t_ = fmt.format(v)
                if float(t_) == 0:
                    t_ = fmt.format(0.0)
                txt.append(minus(t_))
            if P is not None:
                stv = stars(P[i, j])
                if stv:
                    txt.append(stv)
            if txt:
                ax.text(j + 0.5, nr - 1 - i + 0.5, '\n'.join(txt), ha='center', va='center', fontsize=SZ_SMALL, color=star_colour(rgba[:3]), linespacing=0.95)
    ax.set_xlim(0, nc)
    ax.set_ylim(0, nr)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)

def scatter_marker(ax, x, y, disease, color, size, **kw):
    m = MARK[disease]
    return ax.scatter(x, y, s=size * MARK_SCALE[m] ** 2, marker=m, color=color, **kw)

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def snapshot(out_dir, scripts, v):
    arch = pathlib.Path(out_dir) / 'versions'
    arch.mkdir(parents=True, exist_ok=True)
    rec = []
    for s in scripts:
        dst = arch / f'{pathlib.Path(s).stem}_v{v}{pathlib.Path(s).suffix}'
        if dst.exists():
            if sha256(dst) != sha256(s):
                raise FileExistsError(f'{dst} exists with different content; bump the version')
        else:
            shutil.copy2(s, dst)
        rec.append({'sha256': sha256(dst)})
    return rec

def write_selfcheck(path, **kw):
    kw.setdefault('created_at', datetime.datetime.now().astimezone().isoformat(timespec='seconds'))
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(kw, f, ensure_ascii=False, indent=1)

def colours_used(fig):
    used = set()
    for a in fig.findobj():
        for getter in ('get_facecolor', 'get_edgecolor', 'get_color'):
            if hasattr(a, getter):
                try:
                    c = getattr(a, getter)()
                except Exception:
                    continue
                arr = np.atleast_2d(mpl.colors.to_rgba_array(c)) if c is not None and len(np.shape(c)) else None
                if arr is None:
                    try:
                        arr = np.atleast_2d(mpl.colors.to_rgba_array(c))
                    except Exception:
                        continue
                for r in arr:
                    if r[3] > 0:
                        used.add(mpl.colors.to_hex(r[:3]).lower())
    return used
