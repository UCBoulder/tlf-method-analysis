"""Recompute the RWQC AltCalc index-of-agreement table for the ES&T Water manuscript.

Reads the paired values published in TLF-Method/validation.html (the PANELS block behind
Figure 1), so the table is reproducible from what the site already serves.

Willmott d follows EPA-820-R-14-011 Appendix E: BOTH denominator terms use the REFERENCE's
mean of log10 values, not each series' own mean.

Full-set rows reproduce the published table exactly (209/0.96/0.86, 206/0.91/0.70,
153/0.79/0.52), which is the check that the implementation and the data source are right.

FINDING 2026-09-01: the published "within limits of quantification" rows (202/198/148) were
computed at a 1 CFU/100 mL cut, NOT at the 10 CFU/100 mL LOQ the method and the manuscript
state. At the stated 10 CFU LOQ applied to both assays the rows are 190/0.905/0.680,
184/0.881/0.618 and 134/0.705/0.417. Every one still clears the 0.70 acceptance threshold.

Run: python3 scripts/recompute_ia_table.py
"""

import re, io, json, math

P = json.load(open("data/validation_panels.json"))   # the PANELS block served on validation.thelume.ai/TLF-Method
print({k: len(v) for k, v in P.items()})

def ia_r2(pairs):
    """Willmott d and R^2 on log10 pairs. x = reference, y = candidate.
    d = 1 - sum((lx-ly)^2) / sum((|lx-mean(lx)| + |ly-mean(lx)|)^2)   [EPA AltCalc, App. E:
    BOTH denominator terms use the reference's mean]"""
    lx = [math.log10(p[0]) for p in pairs]
    ly = [math.log10(p[1]) for p in pairs]
    mx = sum(lx)/len(lx)
    num = sum((a-b)**2 for a, b in zip(lx, ly))
    den = sum((abs(a-mx) + abs(b-mx))**2 for a, b in zip(lx, ly))
    d = 1 - num/den
    my = sum(ly)/len(ly)
    sxy = sum((a-mx)*(b-my) for a, b in zip(lx, ly))
    sxx = sum((a-mx)**2 for a in lx); syy = sum((b-my)**2 for b in ly)
    r2 = (sxy*sxy)/(sxx*syy)
    return d, r2, len(pairs)

LOQ = 10.0
rows = [("Colilert vs. candidate",            "col2"),
        ("Membrane filtration vs. candidate", "col4"),
        ("Colilert vs. MF (two EPA methods)", "col1")]
print("\n%-38s %5s %6s %6s | %5s %6s %6s" % ("comparison","n","IA","R2","n>=LOQ","IA","R2"))
out = {}
for label, key in rows:
    pts = [(p["x"], p["y"]) for p in P[key]
           if p.get("x") is not None and p.get("y") is not None and p["x"] > 0 and p["y"] > 0]
    d, r2, n = ia_r2(pts)
    # LOQ restriction: BOTH values at or above the 10 CFU/100 mL quantification limit
    sub = [p for p in pts if p[0] >= LOQ and p[1] >= LOQ]
    d2, r22, n2 = ia_r2(sub)
    out[label] = (n, d, r2, n2, d2, r22)
    print("%-38s %5d %6.3f %6.3f | %6d %6.3f %6.3f" % (label, n, d, r2, n2, d2, r22))

print("\npublished table for comparison: 209/0.96/0.86 (202/0.92/0.72); 206/0.91/0.70 (198/0.87/0.60); 153/0.79/0.52 (148/0.72/0.41)")

print("\n=== which LOQ definition reproduces 202 / 198 / 148 ? ===")
defs = {
 "both >= 10":            lambda x,y: x>=10 and y>=10,
 "reference >= 10 only":  lambda x,y: x>=10,
 "candidate >= 10 only":  lambda x,y: y>=10,
 "both >= 1":             lambda x,y: x>=1 and y>=1,
 "reference >= 1 only":   lambda x,y: x>=1,
 "both > 1":              lambda x,y: x>1 and y>1,
}
for name, keep in defs.items():
    ns = []
    for label, key in rows:
        pts = [(p["x"], p["y"]) for p in P[key]
               if p.get("x") and p.get("y") and p["x"] > 0 and p["y"] > 0]
        ns.append(sum(1 for x, y in pts if keep(x, y)))
    print("%-22s -> %s   %s" % (name, ns, "*** MATCH ***" if ns == [202,198,148] else ""))
