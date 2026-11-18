"""Figure 1 for the ES&T Water manuscript: grab-and-read agreement, three panels.

A  Colilert vs candidate      B  MF (S-Endo) vs candidate      C  Colilert vs MF

Pairs come from the PANELS block published in TLF-Method/validation.html, the
same source recompute_ia_table.py reads, so figure and table cannot diverge.
Following Fig. 2's convention: FILLED points are pairs with both assays at or
above the 10 CFU/100 mL LOQ (the acceptance rows of Table 2); OPEN muted
points are sub-LOQ pairs, including substituted non-detects (0.1 recorded in
place of a zero read), shown for completeness rather than hidden.

House style: log-log, dashed 1:1, dotted 126 CFU/100 mL crosshairs, SQUARE axes.
Writes fig_altcalc_ia.pdf/.png into the Overleaf clone.
Run: python3 scripts/make_fig_altcalc.py
"""
import io, json, math, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "output"
GREEN, INK, MUTED = "#1b9e77", "#222222", "#8a8a8a"
LOQ, LIMIT = 10.0, 126.0

P = json.load(open("data/validation_panels.json"))

def ia(pairs):
    """Willmott d on log10 pairs, EPA AltCalc App. E convention (both denominator
    terms use the reference's mean)."""
    lx = [math.log10(x) for x, y in pairs]; ly = [math.log10(y) for x, y in pairs]
    mx = sum(lx) / len(lx)
    return 1 - sum((a - b) ** 2 for a, b in zip(lx, ly)) / \
               sum((abs(a - mx) + abs(b - mx)) ** 2 for a, b in zip(lx, ly))

PANEL = [("A   Colilert vs. candidate", "col2", "Colilert (CFU/100 mL)", "candidate (CFU/100 mL)"),
         ("B   MF (S-Endo) vs. candidate", "col4", "MF (CFU/100 mL)", "candidate (CFU/100 mL)"),
         ("C   Colilert vs. MF", "col1", "Colilert (CFU/100 mL)", "MF (CFU/100 mL)")]

fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.4))
for ax, (title, key, xl, yl) in zip(axes, PANEL):
    pts = [(p["x"], p["y"]) for p in P[key] if p["x"] > 0 and p["y"] > 0]
    acc = [(x, y) for x, y in pts if x >= LOQ and y >= LOQ]
    sub = [(x, y) for x, y in pts if not (x >= LOQ and y >= LOQ)]
    allv = [v for p in pts for v in p]
    lo = 10 ** (math.log10(min(allv)) - 0.2); hi = 10 ** (math.log10(max(allv)) + 0.2)
    ax.plot([lo, hi], [lo, hi], "--", color=INK, lw=0.9, zorder=1)
    ax.axvline(LIMIT, color="#c0392b", ls=":", lw=0.9); ax.axhline(LIMIT, color="#c0392b", ls=":", lw=0.9)
    ax.scatter(*zip(*sub), s=24, facecolors="none", edgecolors=MUTED, lw=0.9,
               zorder=2, label=f"below LOQ / non-detect (n={len(sub)})") if sub else None
    ax.scatter(*zip(*acc), s=28, color=GREEN, alpha=0.8, edgecolors="white", lw=0.4,
               zorder=3, label=f"both assays $\\geq$ LOQ (n={len(acc)})")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.text(0.03, 0.97, f"IA = {ia(acc):.3f} ($\\geq$ LOQ)\nIA = {ia(pts):.3f} (all)",
            transform=ax.transAxes, va="top", fontsize=9)
    ax.legend(loc="lower right", fontsize=7, frameon=False)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/fig_altcalc_ia.{ext}", dpi=300, bbox_inches="tight")
print("written", OUT + "/fig_altcalc_ia.pdf")
