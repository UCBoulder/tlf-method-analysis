"""Figure 2 for the ES&T Water manuscript: the two continuous in-situ arms,
scored on the nowcast's prequential decision log.

Panel A  Boulder Creek, action limit 126 CFU/100 mL.
Panel B  Seine and Marne, action limit 900 CFU/100 mL.
Filled green points are the stated evaluation scope (both assays within
quantification limits and uncensored; site term warmed to the serving gate's
minimum of eight scored grabs). Open muted points are the remaining scored
grabs, shown so the scope hides nothing. Every prediction was published at the
grab's sample hour, before the result was known.

House style follows fig_altcalc_ia: log-log, dashed 1:1, dotted action-limit
crosshairs, SQUARE axes. Reads TLF-Method/nowcast_insitu_pairs.csv (written by
scripts/nowcast_insitu_paper.py). Writes fig_insitu_ia.pdf/.png into the
Overleaf clone. Run: python3 scripts/make_fig_insitu.py
"""
import csv, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "output"
SRC = "data/nowcast_insitu_pairs.csv"
GREEN, INK, MUTED = "#1b9e77", "#222222", "#8a8a8a"

def ia(prs):
    lx = [a for a, b in prs]; ly = [b for a, b in prs]
    mx = sum(lx) / len(lx)
    return 1 - sum((p - q) ** 2 for p, q in zip(lx, ly)) / \
               sum((abs(p - mx) + abs(q - mx)) ** 2 for p, q in zip(lx, ly))

rows = list(csv.DictReader(open(SRC)))
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.7))
for ax, reg, lim, title in [(axes[0], "boulder", 126.0, "A   Boulder Creek, six stations"),
                            (axes[1], "paris", 900.0, "B   Seine and Marne, seven sites")]:
    rs = [r for r in rows if r["region"] == reg]
    sc = [(float(r["obs_log10"]), float(r["stack_yhat_log10"])) for r in rs if r["scoped"] == "1"]
    un = [(float(r["obs_log10"]), float(r["stack_yhat_log10"])) for r in rs if r["scoped"] == "0"]
    allv = [v for p in sc + un for v in p]
    lo = 10 ** (min(allv) - 0.15); hi = 10 ** (max(allv) + 0.15)
    ax.plot([lo, hi], [lo, hi], "--", color=INK, lw=0.9, zorder=1)
    ax.axvline(lim, color="#c0392b", ls=":", lw=0.9, zorder=1)
    ax.axhline(lim, color="#c0392b", ls=":", lw=0.9, zorder=1)
    ax.scatter([10 ** a for a, b in un], [10 ** b for a, b in un], s=26,
               facecolors="none", edgecolors=MUTED, lw=0.9, zorder=2,
               label=f"outside scope (n={len(un)})")
    ax.scatter([10 ** a for a, b in sc], [10 ** b for a, b in sc], s=30,
               color=GREEN, alpha=0.85, edgecolors="white", lw=0.4, zorder=3,
               label=f"evaluation scope (n={len(sc)})")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("Colilert (CFU/100 mL)")
    ax.set_ylabel("nowcast estimate (CFU/100 mL)")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.text(0.03, 0.97, f"IA = {ia(sc):.3f} (evaluation scope)",
            transform=ax.transAxes, va="top", fontsize=9)
    ax.legend(loc="lower right", fontsize=7.5, frameon=False)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/fig_insitu_ia.{ext}", dpi=300, bbox_inches="tight")
print("written", OUT + "/fig_insitu_ia.pdf")
