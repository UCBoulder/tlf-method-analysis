"""IA for the TLF method paper's drinking-water subsection (ES&T Water ew-2026-014215 revision).
Same convention as recompute_ia_table.py and evaluate.py ia_rows: Willmott d on log10, reference
mean in both denominator terms; non-detects enter all-pairs rows as 0.1; censored CBT (>=100) at
the ceiling. Inputs are the published CBT paper's per-sample predictions (new sensor + new day)
and its replicate CBT bags (reference vs reference)."""
import numpy as np, pandas as pd
REV = "data/"
DL = 100.0
def ia(x, y):
    lx, ly = np.log10(x), np.log10(y); m = lx.mean()
    d = 1 - ((lx-ly)**2).sum() / ((abs(lx-m)+abs(ly-m))**2).sum()
    return dict(n=len(lx), ia=round(float(d),3), r2=round(float(np.corrcoef(lx,ly)[0,1]**2),3) if len(lx)>2 else None,
                rmse=round(float(np.sqrt(((lx-ly)**2).mean())),3))
def rows(obs, pred, cens):
    o = np.where(obs <= 0, 0.1, obs); p = np.maximum(pred, 0.1)
    out = {"all samples": ia(o, p)}
    for lo in (1, 10):
        m = (obs >= lo) & (obs < DL) & ~cens & (pred >= lo)
        out[f"both within limits, >={lo}"] = ia(o[m], p[m]) if m.sum() >= 3 else dict(n=int(m.sum()))
    return out
s = pd.read_csv(REV + "cbt_final_samples.csv")
print("columns:", list(s.columns))
print("candidate vs CBT:", rows(s.cbt.values.astype(float), 10**s.pred.values - 1, s.cens.values.astype(bool)))
B = pd.read_csv(REV + "cbt_replicate_bags.csv").dropna(subset=["sample_group"]); B["t"] = pd.to_datetime(B.t_utc, utc=True)
xs, ys, cs = [], [], []
for _, g in B.sort_values("t").groupby("sample_group"):
    for i in range(1, len(g)):
        xs.append(g.cbt.iloc[i-1]); ys.append(g.cbt.iloc[i]); cs.append(bool(g.censored.iloc[i-1] or g.censored.iloc[i]))
print("CBT vs replicate CBT:", rows(np.array(xs,float), np.array(ys,float), np.array(cs)))
obs = s.cbt.values; print("samples", len(obs), "zeros", int((obs<=0).sum()), "censored", int(s.cens.astype(bool).sum()), "1-99 uncensored", int(((obs>=1)&(obs<100)&~s.cens.astype(bool)).sum()))
