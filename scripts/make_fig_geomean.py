"""Figure 3: the method scored at the quantity a criterion is written on,
now on the nowcast's prequential record, both in-situ arms.

A  Non-overlapping 30-day windowed geometric means, record (sensor stack)
   against grabs, normalized by each arm's action limit so Boulder (126) and
   the Seine/Marne (900) share one square log-log panel; crosshairs at 1.0.
B  The frequency half: with n grabs in a window the attainable excursion
   frequencies are multiples of 1/n (the ladder); the 10% cap passes through
   the gaps for every n below 10. Filled points are the observed windows'
   grab-based frequencies; the paired open point below each shows the
   continuous record's own fraction over the identical window.

Reads TLF-Method/nowcast_gm_windows.csv (written by nowcast_insitu_paper.py).
Writes fig_geomean.pdf/.png into the Overleaf clone.
Run: python3 scripts/make_fig_geomean.py
"""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "output"
GREEN, ORANGE, INK, MUTED = "#1b9e77", "#e08a3c", "#222222", "#8a8a8a"
LIM = {"boulder": 126.0, "paris": 900.0}

rows = list(csv.DictReader(open("data/nowcast_gm_windows.csv")))
fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.2, 4.6))

# ── A: windowed GMs, normalized by the action limit ─────────────────────────
for reg, color, marker, label in [("boulder", GREEN, "o", "Boulder Creek (GM/126)"),
                                  ("paris", ORANGE, "^", "Seine and Marne (GM/900)")]:
    rs = [r for r in rows if r["region"] == reg]
    ax.scatter([float(r["gm_grabs"]) / LIM[reg] for r in rs],
               [float(r["gm_record"]) / LIM[reg] for r in rs],
               s=42, color=color, marker=marker, alpha=0.85,
               edgecolors="white", lw=0.4, zorder=3, label=f"{label}, n={len(rs)}")
lo, hi = 0.02, 4
ax.plot([lo, hi], [lo, hi], "--", color=INK, lw=0.9, zorder=1)
ax.axvline(1, color="#c0392b", ls=":", lw=0.9); ax.axhline(1, color="#c0392b", ls=":", lw=0.9)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel("30-day geometric mean of grabs / action limit")
ax.set_ylabel("30-day geometric mean of record / action limit")
ax.set_title("A   windowed geometric means", loc="left", fontsize=11, fontweight="bold")
ax.legend(loc="upper left", fontsize=8, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ── B: the frequency ladder, with paired grab/record observations ──────────
for n in range(1, 32):
    ks = range(0, n + 1)
    bx.scatter([n] * len(ks), [k / n for k in ks], s=3, color=MUTED, alpha=0.35, zorder=1)
bx.axhline(0.10, color="#c0392b", ls=":", lw=1.1, zorder=2)
bx.text(31.6, 0.104, "10% cap", color="#c0392b", fontsize=8, va="bottom", ha="right")
seen = set()
for r in rows:
    n, k = int(r["n_grabs"]), int(r["grab_exc_stv"])
    fr = float(r["rec_frac_stv"])
    color = GREEN if r["region"] == "boulder" else ORANGE
    lg = None
    if r["region"] not in seen:
        lg = "Boulder windows" if r["region"] == "boulder" else "Seine/Marne windows"
        seen.add(r["region"])
    bx.scatter([n], [k / n], s=46, color=color, zorder=4, edgecolors="white", lw=0.4, label=lg)
    bx.scatter([n], [fr], s=30, facecolors="none", edgecolors=color, lw=1.1, zorder=3)
    bx.plot([n, n], [min(k / n, fr), max(k / n, fr)], color=color, lw=0.7, alpha=0.5, zorder=2)
bx.set_xlim(0.5, 32); bx.set_ylim(-0.02, 0.45)
bx.set_xlabel("grabs in the 30-day window (n)")
bx.set_ylabel("excursion frequency")
bx.set_title("B   the frequency half of the criterion", loc="left", fontsize=11, fontweight="bold")
bx.scatter([], [], s=30, facecolors="none", edgecolors=INK, lw=1.1, label="continuous record, same window")
bx.legend(loc="upper right", fontsize=7.5, frameon=False)
for s in ("top", "right"): bx.spines[s].set_visible(False)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/fig_geomean.{ext}", dpi=300, bbox_inches="tight")
print("written", OUT + "/fig_geomean.pdf")
