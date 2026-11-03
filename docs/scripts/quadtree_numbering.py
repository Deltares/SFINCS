"""Schematic of SFINCS quadtree numbering and neighbour indices (2 refinement levels).

Makes docs/figures/SFINCS_quadtree_numbering.png, used in input_quadtree.rst.
Neighbour indices follow the conventions of the quadtree file (as written by HydroMT-SFINCS).

Usage: python quadtree_numbering.py [output.png]
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

if len(sys.argv) > 1:
    out = sys.argv[1]
else:
    out = Path(__file__).resolve().parents[1] / "figures" / "SFINCS_quadtree_numbering.png"

# Coarse grid
x0, y0, dx, dy = 0.0, 0.0, 1.0, 1.0
mmax, nmax = 4, 4
# Cells (n, m) refined per level; neighbours never differ more than 1 level
refined = {1: {(2, 3), (3, 2), (3, 3), (3, 4), (4, 3)}, 2: {(5, 5), (6, 5), (5, 6), (6, 6)}}
nlev = 3

# Build cells by splitting refined cells, then sort: level first, then column-major (m outer, n inner)
todo = [(1, n, m) for m in range(1, mmax + 1) for n in range(1, nmax + 1)]
cells = []
while todo:
    l, n, m = todo.pop()
    if (n, m) in refined.get(l, set()):
        todo += [(l + 1, 2 * n - i, 2 * m - j) for i in (0, 1) for j in (0, 1)]
    else:
        cells.append((l, n, m))
cells.sort(key=lambda c: (c[0], c[2], c[1]))
index = {c: i + 1 for i, c in enumerate(cells)}


def find(l, n, m):
    return index.get((l, n, m), 0)


def neighbours(l, n, m):
    nc, mc = (n + 1) // 2, (m + 1) // 2
    nb = {}
    # (name, same-level target, coarse target if allowed, fine targets)
    dirs = {
        "mu": ((n, m + 1), (nc, mc + 1) if m % 2 == 0 else None, [(2 * n - 1, 2 * m + 1), (2 * n, 2 * m + 1)]),
        "md": ((n, m - 1), (nc, mc - 1) if m % 2 == 1 else None, [(2 * n - 1, 2 * m - 2), (2 * n, 2 * m - 2)]),
        "nu": ((n + 1, m), (nc + 1, mc) if n % 2 == 0 else None, [(2 * n + 1, 2 * m - 1), (2 * n + 1, 2 * m)]),
        "nd": ((n - 1, m), (nc - 1, mc) if n % 2 == 1 else None, [(2 * n - 2, 2 * m - 1), (2 * n - 2, 2 * m)]),
    }
    for key, (same, coarse, fine) in dirs.items():
        flag, i1, i2 = 0, find(l, *same), 0
        if i1 == 0 and l > 1 and coarse is not None and find(l - 1, *coarse):
            flag, i1 = -1, find(l - 1, *coarse)
        if i1 == 0 and l < nlev:
            f1, f2 = find(l + 1, *fine[0]), find(l + 1, *fine[1])
            if f1 or f2:
                flag, i1, i2 = 1, f1, f2
        nb[key] = (flag, i1, i2)
    return nb


def geom(c):
    l, n, m = c
    d = dx / 2 ** (l - 1), dy / 2 ** (l - 1)
    return x0 + (m - 1) * d[0], y0 + (n - 1) * d[1], d[0], d[1]


def centre(c):
    x, y, wx, wy = geom(c)
    return x + 0.5 * wx, y + 0.5 * wy


levcol = {1: "#dbe7f3", 2: "#fbe3c9", 3: "#d8eed2"}
edge = "#3d4a5c"
dircol = {"mu": "#c0392b", "md": "#2471a3", "nu": "#1e8449", "nd": "#7d3c98"}


def draw_grid(ax, labels="full", highlight=None, nbrs=()):
    for c in cells:
        x, y, wx, wy = geom(c)
        fc, lw = levcol[c[0]], 1.2
        if c == highlight:
            fc, lw = "#ffd84d", 2.5
        ax.add_patch(Rectangle((x, y), wx, wy, fc=fc, ec=edge, lw=lw, zorder=1))
        cx, cy = centre(c)
        fs = {1: 15, 2: 11, 3: 9}[c[0]]
        if labels == "full":
            ax.text(cx, cy + 0.06 * wy, str(index[c]), ha="center", va="center", fontsize=fs, weight="bold", zorder=3, clip_on=True)
            ax.text(cx, cy - 0.27 * wy, f"n={c[1]}, m={c[2]}", ha="center", va="center",
                    fontsize={1: 8.5, 2: 6.5, 3: 4.0}[c[0]], color="#333", zorder=3, clip_on=True)
        else:
            ax.text(cx, cy, str(index[c]), ha="center", va="center", fontsize=fs - 2, weight="bold", zorder=3,
                    color="#000" if (highlight == c or index[c] in nbrs) else "#888", clip_on=True)


fig = plt.figure(figsize=(18, 13))
gs = fig.add_gridspec(3, 3, width_ratios=[1.45, 0.55, 0.36], hspace=0.14, wspace=0.04,
                      left=0.01, right=0.995, top=0.95, bottom=0.02)

# Panel (a): overview
ax = fig.add_subplot(gs[:, 0])
draw_grid(ax)
for m in range(1, mmax + 1):
    ax.text(x0 + (m - 0.5) * dx, y0 - 0.12, f"m = {m}", ha="center", va="top", fontsize=11)
for n in range(1, nmax + 1):
    ax.text(x0 - 0.1, y0 + (n - 0.5) * dy, f"n = {n}", ha="right", va="center", fontsize=11, rotation=90)
ax.plot(x0, y0, "ko", ms=8, zorder=5)
ax.annotate("(x0, y0)", (x0, y0), (x0 - 0.42, y0 - 0.42), fontsize=12, ha="center",
            arrowprops=dict(arrowstyle="-", color="k", lw=1))
# m and n direction arrows
ax.annotate("", (x0 + 1.6, y0 - 0.55), (x0, y0 - 0.55), arrowprops=dict(arrowstyle="-|>", lw=2, color="k"))
ax.text(x0 + 1.65, y0 - 0.55, "m direction (u)", va="center", fontsize=12)
ax.annotate("", (x0 - 0.5, y0 + 1.6), (x0 - 0.5, y0), arrowprops=dict(arrowstyle="-|>", lw=2, color="k"))
ax.text(x0 - 0.5, y0 + 1.68, "n direction (v)", ha="center", va="bottom", fontsize=12, rotation=90)
# dx / dy dimensions of the coarse grid
ax.annotate("", (x0 + 4 * dx, y0 + nmax * dy + 0.12), (x0 + 3 * dx, y0 + nmax * dy + 0.12),
            arrowprops=dict(arrowstyle="<->", lw=1))
ax.text(x0 + 3.5 * dx, y0 + nmax * dy + 0.16, "dx", ha="center", va="bottom", fontsize=11)
ax.annotate("", (x0 + mmax * dx + 0.12, y0 + dy), (x0 + mmax * dx + 0.12, y0),
            arrowprops=dict(arrowstyle="<->", lw=1))
ax.text(x0 + mmax * dx + 0.18, y0 + 0.5 * dy, "dy", va="center", fontsize=11, rotation=90)
# Legend
for l, txt in [(1, "level 1 (dx, dy)"), (2, "level 2 (dx/2, dy/2)"), (3, "level 3 (dx/4, dy/4)")]:
    yl = y0 - 0.38 - 0.2 * (l - 1)
    ax.add_patch(Rectangle((x0 + 2.9, yl - 0.07), 0.18, 0.14, fc=levcol[l], ec=edge))
    ax.text(x0 + 3.13, yl, txt, va="center", fontsize=10.5)
ax.text(x0 - 0.9, y0 - 1.05,
        "Cell numbers: level 1 first, then level 2, then level 3; within a level, n runs fastest (column by column).\n"
        "n, m are the indices within the cell's own level (level l has 2^(l-1) * nmax rows and 2^(l-1) * mmax columns).",
        fontsize=10.5, va="top", color="#333", linespacing=1.5)
ax.text(x0 - 0.9, y0 - 1.4,
        "mu / md / nu / nd :  0 = neighbour at same level,  -1 = neighbour is coarser,  1 = neighbour is finer\n"
        "mu1, mu2 ... nd1, nd2 :  cell number of the neighbour(s);  0 = no neighbour\n"
        "With finer neighbours:  mu1/md1 = lower cell, mu2/md2 = upper cell;  nu1/nd1 = left cell, nu2/nd2 = right cell",
        fontsize=10.5, va="top", ha="left", color="#222", linespacing=1.5)
ax.set_xlim(x0 - 0.95, x0 + mmax * dx + 0.4)
ax.set_ylim(y0 - 2.15, y0 + nmax * dy + 0.4)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("(a) Quadtree grid with 2 refinement levels: origin, directions and cell numbering", fontsize=13, loc="left")


# Panels (b), (c): neighbour stencils
def stencil(ax, axt, c, title, half):
    ic = index[c]
    nb = neighbours(*c)
    targets = {k: [i for i in v[1:] if i > 0] for k, v in nb.items()}
    draw_grid(ax, labels="short", highlight=c, nbrs=[i for v in targets.values() for i in v])
    cx, cy = centre(c)
    for key, ids in targets.items():
        for j, i in enumerate(ids):
            tx, ty = centre(cells[i - 1])
            ax.add_patch(FancyArrowPatch((cx, cy), (tx, ty), arrowstyle="-|>", mutation_scale=16,
                                         color=dircol[key], lw=2.2, shrinkA=14, shrinkB=16, zorder=4))
            sub = f"{key}{j + 1}"
            # label beside the arrow (perpendicular offset), at its midpoint
            vx, vy = tx - cx, ty - cy
            ln = (vx ** 2 + vy ** 2) ** 0.5
            px, py = -vy / ln, vx / ln
            if px < 0 or (px == 0 and py < 0):
                px, py = -px, -py
            off, frac = 0.09 * half, 0.55
            if len(ids) == 2:
                # two finer neighbours: *1 label below/left, *2 label above/right
                sgn = -1 if j == 0 else 1
                px, py = (0, sgn) if key in ("mu", "md") else (sgn, 0)
                off, frac = 0.14 * half, (0.72 if key in ("mu", "md") else 0.36)
            ax.text(cx + frac * vx + off * px, cy + frac * vy + off * py, sub, color=dircol[key],
                    fontsize=11, weight="bold", ha="center", va="center", zorder=6,
                    bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9))
    # keep the view inside the grid
    cxv = min(max(cx, x0 + half - 0.1 * dx), x0 + mmax * dx - half + 0.1 * dx)
    ax.set_xlim(cxv - half, cxv + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("#999")
    ax.set_title(title, fontsize=13, loc="left")
    flagtxt = {0: "same level", -1: "coarser", 1: "finer"}
    lines = [f"cell {ic}", f"level {c[0]}, n={c[1]}, m={c[2]}", ""]
    for key in ["mu", "md", "nu", "nd"]:
        f, i1, i2 = nb[key]
        lines.append(f"{key}  = {f:2d}  ({flagtxt[f]})")
        lines.append(f"{key}1 = {i1:2d}")
        lines.append(f"{key}2 = {i2:2d}")
        lines.append("")
    axt.axis("off")
    axt.text(0.04, 0.5, "\n".join(lines[:-1]), transform=axt.transAxes, fontsize=11, family="monospace",
             va="center", ha="left", bbox=dict(fc="#f7f7f7", ec="#bbb", boxstyle="round,pad=0.6"))


for row, (c, txt, half) in enumerate([
        ((1, 2, 2), "next to finer cells", 1.4 * dx),
        ((2, 6, 4), "next to coarser and finer cells", 0.95 * dx),
        ((3, 9, 9), "next to coarser cells", 0.55 * dx)]):
    stencil(fig.add_subplot(gs[row, 1]), fig.add_subplot(gs[row, 2]), c,
            f"({'bcd'[row]}) Level-{c[0]} cell {index[c]} {txt}", half)

fig.savefig(out, dpi=150)
