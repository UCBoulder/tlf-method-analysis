"""Day-blocked cross-validation of the Boulder Creek continuous in-situ arm.

The manuscript reports that arm in-sample (IA 0.79, R2 0.47). This holds it out.

Serving expression, from shared/ecoli-model.js FIELD_MODEL:
    log10(cfu) = fe[bc] + slope_tlf[bc]*tlf + slope_tof[bc]*tofr
                 + b_temp*dT + b_tlf_temp*tlf*dT
with tlf = mon2 - base20[bc], tofr = (tof - foulingOffset) - btof[bc], dT = T - 20.
/api/field-refit serves tlf, tofr, dT and y = log10(ec) per grab, so the features need
no reconstruction.

PROTOCOL: leave-one-DAY-out. Every grab sharing a calendar date is held out together,
so a fold can never be scored on a day it was fitted on. Temperature coefficients are
held at their shipped values in every fold, matching production; only the per-sensor
levels and slopes are refitted. Two parameterisations are scored, the deployed one
(per-sensor slopes) and a reduced one (shared slopes), because 27 free parameters on
~145 points is where per-sensor slopes stop being estimable inside a fold.

RESULT 2026-09-01, n=155 across 6 sensors and 34 days, 33 exceedances at 126 CFU/100 mL:

                                    IA      R2     RMSE  | sens   spec   BA
  in-sample,  per-sensor slopes    0.799   0.425   0.432 | 0.758  0.631  0.694
  DAY-BLOCKED, per-sensor slopes   0.711   0.261   0.527 | 0.727  0.598  0.663
  in-sample,  shared slopes        0.770   0.377   0.449 | 0.758  0.623  0.690
  DAY-BLOCKED, shared slopes       0.741   0.320   0.479 | 0.727  0.607  0.667

THE ARM SURVIVES. Held out it is IA 0.711 in the deployed parameterisation and 0.741
with shared slopes; both clear the 0.70 acceptance threshold. The in-sample 0.799 here
matches the 0.79 published on n=145, so the record has grown, not the method changed.

PER-SENSOR SLOPES OVERFIT. They win in sample (0.799 vs 0.770) and LOSE held out
(0.711 vs 0.741). The simpler shared-slope form generalises better on this record.

R2 falls to 0.26-0.32 held out, well below the pathway's 0.60 secondary route, which the
manuscript already states for the in-sample figure. Classification barely moves
(BA 0.694 -> 0.663), so the class decision is far more robust than the magnitude fit.

Run: python3 scripts/boulder_insitu_cv.py   (expects fr.json from GET /api/field-refit)
"""
import json, math
from collections import defaultdict

S = "data"
BOULDER = {"50046", "50048", "50052", "50059", "50062", "50066"}
B_TEMP, B_TLF_TEMP = -0.0204932, -0.0001321
LIMIT, CUT = 126.0, 79.0   # FIELD_MODEL.exc_t

rows = [r for r in json.load(open(S + "/boulder_field_refit_20260923.json"))["grabs"]
        if str(r.get("bc")) in BOULDER and r.get("y") is not None
        and r.get("tlf") is not None and r.get("tofr") is not None and r.get("dT") is not None]
for r in rows:
    r["bc"] = str(r["bc"])
    # move the fixed temperature terms to the left-hand side; they are never refitted
    r["resid_y"] = r["y"] - B_TEMP*r["dT"] - B_TLF_TEMP*r["tlf"]*r["dT"]
print("Boulder grabs served by /api/field-refit: %d across %d sensors, %d days"
      % (len(rows), len({r['bc'] for r in rows}), len({r['date'] for r in rows})))

def solve(A, b):
    n = len(b); M = [A[i][:] + [b[i]] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-12: return None
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c]/M[c][c]
                for k in range(c, n+1): M[r][k] -= f*M[c][k]
    return [M[i][n]/M[i][i] for i in range(n)]

def fit(train, per_sensor_slopes):
    """Design: per-sensor intercepts, then either per-sensor or shared tlf/tof slopes."""
    bcs = sorted({r["bc"] for r in train}); bi = {b: i for i, b in enumerate(bcs)}
    k = len(bcs) + (2*len(bcs) if per_sensor_slopes else 2)
    X, y = [], []
    for r in train:
        row = [0.0]*k
        row[bi[r["bc"]]] = 1.0
        if per_sensor_slopes:
            row[len(bcs) + 2*bi[r["bc"]]]     = r["tlf"]
            row[len(bcs) + 2*bi[r["bc"]] + 1] = r["tofr"]
        else:
            row[len(bcs)]   = r["tlf"]
            row[len(bcs)+1] = r["tofr"]
        X.append(row); y.append(r["resid_y"])
    XtX = [[sum(X[i][a]*X[i][b] for i in range(len(X))) for b in range(k)] for a in range(k)]
    for i in range(k): XtX[i][i] += 1e-9          # tiny ridge: keeps a fold with a rank-deficient sensor solvable
    beta = solve(XtX, [sum(X[i][a]*y[i] for i in range(len(X))) for a in range(k)])
    if beta is None: return None
    return bcs, bi, beta, per_sensor_slopes

def predict(model, r):
    bcs, bi, beta, ps = model
    if r["bc"] not in bi: return None
    v = beta[bi[r["bc"]]]
    if ps:
        v += beta[len(bcs) + 2*bi[r["bc"]]]*r["tlf"] + beta[len(bcs) + 2*bi[r["bc"]] + 1]*r["tofr"]
    else:
        v += beta[len(bcs)]*r["tlf"] + beta[len(bcs)+1]*r["tofr"]
    return v + B_TEMP*r["dT"] + B_TLF_TEMP*r["tlf"]*r["dT"]

def score(pairs, label):
    """pairs = [(y_true_log10, y_pred_log10, ec_true)]"""
    lx = [a for a, _, _ in pairs]; ly = [b for _, b, _ in pairs]
    mx = sum(lx)/len(lx)
    ia = 1 - sum((a-b)**2 for a, b in zip(lx, ly)) / sum((abs(a-mx)+abs(b-mx))**2 for a, b in zip(lx, ly))
    my = sum(ly)/len(ly)
    sxy = sum((a-mx)*(b-my) for a, b in zip(lx, ly))
    r2 = sxy*sxy/(sum((a-mx)**2 for a in lx)*sum((b-my)**2 for b in ly))
    rmse = math.sqrt(sum((a-b)**2 for a, b in zip(lx, ly))/len(lx))
    obs = [ec >= LIMIT for _, _, ec in pairs]; prd = [10**b >= CUT for _, b, _ in pairs]
    tp = sum(1 for o, p in zip(obs, prd) if o and p); fn = sum(1 for o, p in zip(obs, prd) if o and not p)
    tn = sum(1 for o, p in zip(obs, prd) if not o and not p); fp = sum(1 for o, p in zip(obs, prd) if not o and p)
    se = tp/(tp+fn) if tp+fn else float("nan"); sp = tn/(tn+fp) if tn+fp else float("nan")
    print("%-34s n=%3d  IA=%.3f  R2=%.3f  RMSE=%.3f | sens %.3f spec %.3f BA %.3f"
          % (label, len(pairs), ia, r2, rmse, se, sp, (se+sp)/2))
    return ia

by_day = defaultdict(list)
for r in rows: by_day[r["date"]].append(r)
print("exceedances at %d: %d of %d\n" % (LIMIT, sum(1 for r in rows if r["ec"] >= LIMIT), len(rows)))

for ps, name in ((True, "per-sensor slopes (deployed)"), (False, "shared slopes")):
    m = fit(rows, ps)
    score([(r["y"], predict(m, r), r["ec"]) for r in rows], "IN-SAMPLE, " + name)
    held = []
    for day, fold in by_day.items():
        train = [r for r in rows if r["date"] != day]
        mf = fit(train, ps)
        if mf is None: continue
        for r in fold:
            p = predict(mf, r)
            if p is not None and math.isfinite(p): held.append((r["y"], p, r["ec"]))
    score(held, "DAY-BLOCKED CV, " + name)
    print()
