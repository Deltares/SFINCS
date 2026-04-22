"""
Visualisation of SFINCS vegetation drag lookup table concept.
Run with:  python vegetation_lookup_figure.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# --------------------------------------------------------------------------
# Example vegetation: 3 vertical sections stacked from bed
# --------------------------------------------------------------------------
sections = [
    {"ah": 0.5, "cd_wd": 0.8,  "label": "section 1\n(dense canopy base)", "color": "#4CAF50"},
    {"ah": 0.8, "cd_wd": 0.4,  "label": "section 2\n(mid canopy)",         "color": "#8BC34A"},
    {"ah": 0.4, "cd_wd": 0.15, "label": "section 3\n(sparse canopy top)",  "color": "#CDDC39"},
]

# Cumulative section boundaries
bottoms = [0.0]
for s in sections:
    bottoms.append(bottoms[-1] + s["ah"])
tops = bottoms[1:]
h_veg_total = bottoms[-1]   # = 1.7 m

# --------------------------------------------------------------------------
# Build lookup table (vegetation_nlookup = 20)
# --------------------------------------------------------------------------
nlookup = 20
dh = h_veg_total / nlookup
h_nodes = np.arange(0, nlookup + 1) * dh   # shape (21,)

table = np.zeros(nlookup + 1)
for k in range(1, nlookup + 1):
    h_k = k * dh
    sec_bot = 0.0
    for s in sections:
        sec_top = sec_bot + s["ah"]
        table[k] += s["cd_wd"] * max(0.0, min(sec_top, h_k) - sec_bot)
        sec_bot = sec_top

slope_table = np.diff(table)   # length nlookup

# --------------------------------------------------------------------------
# Runtime lookup example: water depth hu
# --------------------------------------------------------------------------
hu = 1.1   # example water depth at a timestep
hu_eff = min(hu, h_veg_total)
frac = hu_eff / dh
ik = min(int(frac), nlookup - 1)
frac_r = frac - ik
veg_cd_eff = table[ik] + frac_r * slope_table[ik]

# --------------------------------------------------------------------------
# Plot
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(13, 6))
fig.suptitle("SFINCS vegetation drag — lookup table concept", fontsize=13, fontweight="bold")
gs = GridSpec(1, 3, figure=fig, wspace=0.40, left=0.06, right=0.97, top=0.88, bottom=0.10)

# ── Panel A: vegetation sections ─────────────────────────────────────────
ax1 = fig.add_subplot(gs[0])
ax1.set_title("A  Vegetation sections\n(single uv-point)", fontsize=10)

for i, s in enumerate(sections):
    ax1.barh(
        bottoms[i] + s["ah"] / 2, 1.0,
        height=s["ah"], left=0,
        color=s["color"], edgecolor="k", linewidth=0.8, alpha=0.85
    )
    ax1.text(
        0.5, bottoms[i] + s["ah"] / 2,
        f'cd·b·N = {s["cd_wd"]:.2f}\nah = {s["ah"]} m',
        ha="center", va="center", fontsize=7.5
    )
    # section boundary lines
    ax1.axhline(bottoms[i], color="gray", lw=0.7, ls="--")

ax1.axhline(h_veg_total, color="gray", lw=0.7, ls="--")
ax1.set_xlim(0, 1); ax1.set_xticks([])
ax1.set_ylim(-0.15, 2.1)
ax1.set_ylabel("Height above bed (m)")
ax1.set_xlabel("← stem density schematic →")

# water surface at hu
ax1.axhline(hu, color="#1565C0", lw=2, ls="-", label=f"water depth hu = {hu} m")
ax1.fill_betweenx([0, hu], 0, 1, color="#90CAF9", alpha=0.25)
ax1.legend(fontsize=7.5, loc="upper right")

# ── Panel B: lookup table ─────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[1])
ax2.set_title("B  Pre-computed lookup table\n(built once in initialize_vegetation)", fontsize=10)

ax2.step(table, h_nodes, where="post", color="k", lw=1.2, label="table values")
ax2.plot(table, h_nodes, "o", color="k", ms=4, zorder=5)

# shade each section contribution band
sec_bot = 0.0
for s in sections:
    sec_top = sec_bot + s["ah"]
    ax2.axhspan(sec_bot, sec_top, color=s["color"], alpha=0.20)
    sec_bot = sec_top

# interpolation at hu
ax2.axhline(hu_eff, color="#1565C0", lw=1.5, ls="--", label=f"hu = {hu} m")
ax2.plot(veg_cd_eff, hu_eff, "*", color="red", ms=12, zorder=10,
         label=f"veg_cd_eff = {veg_cd_eff:.3f} m²/s²  (interpolated)")

# show interpolation bracket
ax2.plot([table[ik], table[ik+1]], [h_nodes[ik], h_nodes[ik+1]],
         "r-", lw=1.5, label="linear interpolation")
ax2.plot([table[ik], table[ik+1]], [h_nodes[ik], h_nodes[ik+1]],
         "rs", ms=6)

ax2.set_xlabel("Cumulative drag integral\n∑ cd·b·N · submerged thickness  (m²/s²)")
ax2.set_ylabel("Depth level h_k (m)")
ax2.set_ylim(-0.15, 2.1)
ax2.legend(fontsize=7.5, loc="lower right")

# ── Panel C: flux update (implicit) ──────────────────────────────────────
ax3 = fig.add_subplot(gs[2])
ax3.set_title("C  Explicit flux update\n(every timestep in compute_fluxes)", fontsize=10)
ax3.axis("off")

textblock = (
    r"$\bf{Runtime\ lookup\ (O(1),\ no\ inner\ loop):}$" + "\n\n"
    r"$\mathrm{frac} = \min(h_u,\; h_{veg}) \;/\; \Delta h$" + "\n"
    r"$ik = \lfloor \mathrm{frac} \rfloor$" + "\n"
    r"$cd_{eff} = \mathrm{table}[ik] + (\mathrm{frac}-ik)\times\mathrm{slope}[ik]$" + "\n\n"
    r"$\bf{Explicit\ momentum\ update:}$" + "\n\n"
    r"$F_{veg} = -\phi\;cd_{eff}\;u_0\;|u_0|$" + "\n\n"
    r"$q^{n+1} = \dfrac{q^n + (F_{ext} + F_{veg})\,\Delta t}"
    r"{1 + \dfrac{g\,n^2\,|q|}{h_u^{7/3}}\,\Delta t}$"
    + "\n\n"
    r"$\bf{Key\ properties:}$" + "\n"
    "• No inner loop over sections at runtime\n"
    "• Linear interpolation between table bins\n"
    "• Sections stacked from bed upward\n"
    "  (consistent with SnapWave swvegatt)\n"
    r"• $\phi$ = wet fraction (subgrid correction)"
)

ax3.text(0.03, 0.97, textblock, transform=ax3.transAxes,
         fontsize=9, va="top", ha="left",
         bbox=dict(boxstyle="round,pad=0.5", fc="#F5F5F5", ec="#BDBDBD"),
         linespacing=1.6)

plt.savefig("vegetation_lookup_figure.png", dpi=150)
print("Saved: vegetation_lookup_figure.png")
plt.show()
