"""ALL manuscript numbers for the nowcast in-situ arms (ES&T Water paper).

The paired record is the public immutable decisions log
(validation.thelume.ai/api/nowcast/decisions): one row per scored grab carrying
every arm's prediction as published at the grab's sample hour, prequentially.
The scored arm is the sensor stack. The stated evaluation scope, anchored in
rules that predate this analysis:
  - pairs within the quantification limits of both assays, uncensored
    (manuscript Sections 2.4/2.5; the AMC eligibility rule), and
  - predictions made after the site's paired record reached the serving gate's
    minimum of eight scored grabs (lume_nowcast engine/gate.py SITE_MIN; the
    method specification requires a site term fitted to that water's prior FIB
    record).
Unrestricted figures are computed alongside; the warm-up sweep is reported.

Outputs (written into output/; the committed data/ copies are the versions the paper reports):
  TLF-Method/nowcast_insitu_stats.json   every number the manuscript cites
  TLF-Method/nowcast_insitu_pairs.csv    per-grab pairs (region, scoped flag)
  TLF-Method/nowcast_gm_windows.csv      30-day windowed GM comparison rows
Run: python3 scripts/nowcast_insitu_paper.py
"""
import csv, json, math, subprocess, datetime as dt
from statistics import mean

BASE = "https://validation.thelume.ai/api/nowcast"
OUT = "output"
LIMITS = {"boulder": 126.0, "paris": 900.0}
STV = {"boulder": 410.0, "paris": 900.0}   # frequency half: RWQC STV / EU single-sample
REGION_OF = lambda sid: "boulder" if sid.startswith("BC_") else (
    "paris" if sid.startswith(("PAR_", "MAR_")) else None)
SITES = {"boulder": ["BC_13", "BC_30", "BC_55", "BC_CU", "BC_Can", "BC_Eben"],
         "paris":   ["PAR_BrasMarie", "PAR_Bercy", "PAR_Grenelle", "PAR_Villette",
                     "MAR_Neuilly", "MAR_StMaur", "MAR_Champigny"]}
WARM_MIN = 8   # engine gate SITE_MIN

def get(url):
    return json.loads(subprocess.run(["curl", "-s", url], capture_output=True, text=True).stdout)

def willmott(prs):
    if len(prs) < 3: return float("nan")
    ob = mean(o for o, p in prs)
    den = sum((abs(p - ob) + abs(o - ob)) ** 2 for o, p in prs)
    return 1 - sum((p - o) ** 2 for o, p in prs) / den if den else float("nan")

def r2(prs):
    ob = mean(o for o, p in prs)
    sst = sum((o - ob) ** 2 for o, p in prs)
    return 1 - sum((o - p) ** 2 for o, p in prs) / sst if sst else float("nan")

def rmse(prs):
    return math.sqrt(mean((o - p) ** 2 for o, p in prs))

def r3(x): return None if x != x else round(x, 3)

# ── decisions: the paired record ────────────────────────────────────────────
dec = get(f"{BASE}/decisions")["decisions"]
for r in dec:
    r["arms"] = json.loads(r["arms_json"]) if isinstance(r["arms_json"], str) else r["arms_json"]
dec = [r for r in dec if REGION_OF(r["site_id"])
       and isinstance(r.get("obs"), (int, float))
       and (r.get("arms") or {}).get("stack")
       and isinstance(r["arms"]["stack"].get("yhat"), (int, float))]
dec.sort(key=lambda r: r["grab_ts"])
seen = {}
for r in dec:
    r["rank"] = seen.get(r["site_id"], 0)
    seen[r["site_id"]] = r["rank"] + 1
    r["region"] = REGION_OF(r["site_id"])
    r["loq_ok"] = (not r["censored"]) and (r["ecoli_cfu"] or 0) >= 10
    r["scoped"] = r["loq_ok"] and r["rank"] >= WARM_MIN

stats = {"generated_from": f"{BASE}/decisions", "scope":
         f"LOQ>=10 & uncensored & site rank>={WARM_MIN} (gate SITE_MIN)", "regions": {}}

for reg in ("boulder", "paris"):
    rs = [r for r in dec if r["region"] == reg]
    S = {}
    for label, sel in [("unrestricted", rs), ("scoped", [r for r in rs if r["scoped"]])]:
        prs = [(r["obs"], r["arms"]["stack"]["yhat"]) for r in sel]
        exc = sum(r["exceeded"] for r in sel)
        miss = sum(1 for r in sel if r["exceeded"] and not r["arms"]["stack"]["alarm"])
        fa = sum(1 for r in sel if not r["exceeded"] and r["arms"]["stack"]["alarm"])
        clnN = len(sel) - exc
        sens = (exc - miss) / exc if exc else float("nan")
        spec = (clnN - fa) / clnN if clnN else float("nan")
        S[label] = {"n": len(prs), "exc": exc, "IA": r3(willmott(prs)), "R2": r3(r2(prs)),
                    "RMSE": r3(rmse(prs)), "sens": r3(sens), "spec": r3(spec),
                    "BA": r3((sens + spec) / 2)}
    # calibration arm on the same scoped grabs: does the sensor earn its place?
    sel = [r for r in rs if r["scoped"] and (r["arms"].get("calibration") or {}).get("yhat") is not None]
    prs = [(r["obs"], r["arms"]["calibration"]["yhat"]) for r in sel]
    costs = {}
    for a in ("calibration", "stack", "regression"):
        have = [r for r in rs if r["arms"].get(a)]
        costs[a] = 3 * sum(1 for r in have if r["exceeded"] and not r["arms"][a]["alarm"]) \
                   + sum(1 for r in have if not r["exceeded"] and r["arms"][a]["alarm"])
    S["calibration_scoped"] = {"n": len(prs), "IA": r3(willmott(prs)), "RMSE": r3(rmse(prs))}
    S["decision_cost_all_grabs"] = costs
    S["sweep"] = {k: {"n": len(p), "IA": r3(willmott([(r['obs'], r['arms']['stack']['yhat']) for r in p]))}
                  for k in (4, 6, 8, 10, 12)
                  for p in [[r for r in rs if r["loq_ok"] and r["rank"] >= k]]}
    S["dates"] = [min(r["grab_ts"][:10] for r in rs), max(r["grab_ts"][:10] for r in rs)]
    stats["regions"][reg] = S

# ── windowed 30-day geometric means, stack arm, per region ─────────────────
gm_rows = []
freq = {}
for reg, sites in SITES.items():
    lim = LIMITS[reg]
    grec, ggrab = [], []
    hrs_above = 0; hrs_tot = 0; g_above = 0; g_tot = 0
    for sid in sites:
        ser = get(f"{BASE}/series?site={sid}&days=400")["series"]
        rows = []
        for x in ser:
            try:
                a = json.loads(x["arms_json"]) if isinstance(x["arms_json"], str) else x["arms_json"]
                y = (a or {}).get("stack", {}).get("yhat")
            except Exception:
                y = None
            if isinstance(y, (int, float)):
                rows.append((ts(x["ts"]), y))
        grabs = [(ts(g["ts"]), g["ecoli_cfu"]) for g in get(f"{BASE}/grabs?site={sid}")["grabs"]
                 if isinstance(g.get("ecoli_cfu"), (int, float)) and g["ecoli_cfu"] > 0]
        if not rows: continue
        # non-overlapping 30-day windows from the first series hour
        t0, t1 = rows[0][0], rows[-1][0]
        w = t0
        while w + dt.timedelta(days=30) <= t1 + dt.timedelta(hours=1):
            hi = w + dt.timedelta(days=30)
            rv = [y for t, y in rows if w <= t < hi]
            gv = [v for t, v in grabs if w <= t < hi]
            if len(rv) >= 0.5 * 30 * 24 and len(gv) >= 4:   # feed's own coverage + >=4 grab rule
                rgm = 10 ** mean(rv)
                ggm = 10 ** mean(math.log10(v) for v in gv)
                stv = STV[reg]
                gm_rows.append({"region": reg, "site": sid, "start": w.date().isoformat(),
                                "n_hours": len(rv), "n_grabs": len(gv),
                                "gm_record": round(rgm, 1), "gm_grabs": round(ggm, 1),
                                "grab_exc_stv": sum(1 for v in gv if v >= stv),
                                "rec_frac_stv": round(sum(1 for y in rv if 10 ** y >= stv) / len(rv), 4)})
                grec.append(math.log10(rgm)); ggrab.append(math.log10(ggm))
            w = hi
        hrs_above += sum(1 for t, y in rows if 10 ** y >= lim); hrs_tot += len(rows)
        g_above += sum(1 for t, v in grabs if v >= lim); g_tot += len(grabs)
    prs = list(zip(ggrab, grec))
    stats["regions"][reg]["gm30"] = {
        "windows": len(prs), "IA": r3(willmott(prs)), "RMSE": r3(rmse(prs)) if prs else None,
        "agree_above_limit": sum(1 for o, p in prs
                                 if (10 ** o >= lim) == (10 ** p >= lim)) if prs else 0,
        "exceed_frac_record": r3(hrs_above / hrs_tot) if hrs_tot else None,
        "exceed_frac_grabs": r3(g_above / g_tot) if g_tot else None}
    freq[reg] = (hrs_above, hrs_tot, g_above, g_tot)

# ── dump ────────────────────────────────────────────────────────────────────
with open(f"{OUT}/nowcast_insitu_pairs.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["region", "site", "grab_ts", "ecoli_cfu", "obs_log10",
                "stack_yhat_log10", "exceeded", "censored", "rank", "scoped"])
    for r in dec:
        w.writerow([r["region"], r["site_id"], r["grab_ts"], r["ecoli_cfu"],
                    round(r["obs"], 4), round(r["arms"]["stack"]["yhat"], 4),
                    r["exceeded"], r["censored"], r["rank"], int(r["scoped"])])
with open(f"{OUT}/nowcast_gm_windows.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(gm_rows[0].keys())); w.writeheader()
    for row in gm_rows: w.writerow(row)
with open(f"{OUT}/nowcast_insitu_stats.json", "w") as f:
    json.dump(stats, f, indent=1)
print(json.dumps(stats, indent=1))
