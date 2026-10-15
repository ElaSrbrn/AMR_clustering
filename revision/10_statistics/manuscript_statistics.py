#!/usr/bin/env python3
"""Reproduce the statistics reported in the revised manuscript from the Supplementary Tables workbook.

Usage:  python3 manuscript_statistics.py Supplementary_Tables_S1-S12_clean.xlsx [outdir]

Sections
  1. Structural concordance of carbapenemase-encoding plasmids (Supplementary Table S3b):
     proportions with exact Clopper-Pearson 95% CIs, exact McNemar tests between consecutive
     depths and between assemblers, logistic GEE on log2(depth) (exchangeable working correlation,
     clustered by plasmid, robust SEs), per assembler and jointly with assembler as covariate.
  2. Read-based vs assembly-based carbapenemase typing (Supplementary Table S8a): exact McNemar.
  3. Actual vs nominal sequencing depth (Supplementary Table S6a): median, IQR, range.

Tested with Python 3.13, pandas 2.x, statsmodels 0.15.0, SciPy 1.18.1, openpyxl 3.1.
"""
import sys, os
import numpy as np
import pandas as pd
from scipy.stats import binomtest, beta
import statsmodels.api as sm
import statsmodels.formula.api as smf

CONCORDANT = {"concordant", "concordant_not_circular"}


def read_sheet(xlsx, sheet, first_header):
    """Read a sheet whose column header row starts with the cell value `first_header`
    (title and note rows above it are skipped)."""
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    hdr = raw.index[raw.iloc[:, 0].astype(str).str.strip() == first_header][0]
    return pd.read_excel(xlsx, sheet_name=sheet, header=hdr)


def clopper_pearson(k, n, alpha=0.05):
    lo = beta.ppf(alpha / 2, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - alpha / 2, k + 1, n - k) if k < n else 1.0
    return lo, hi


def mcnemar_exact(b, c):
    """b, c = discordant pair counts; two-sided exact binomial test."""
    return 1.0 if b + c == 0 else binomtest(b, b + c, 0.5).pvalue


def concordance(xlsx, out):
    d = read_sheet(xlsx, "S3b_concordance_long", "plasmid_id")
    d = d[["plasmid_id", "assembler", "depth", "class"]].dropna()
    d["conc"] = d["class"].isin(CONCORDANT).astype(int)
    d["depth_x"] = d["depth"].str.replace("x", "", regex=False).astype(int)
    depths = sorted(d["depth_x"].unique())
    n_plasmids = d["plasmid_id"].nunique()
    print(f"\n== 1. Structural concordance ({n_plasmids} plasmids) ==")

    rows = []
    for a in ["flye", "hybracter"]:
        for dp in depths:
            s = d[(d.assembler == a) & (d.depth_x == dp)]
            k, n = int(s.conc.sum()), len(s)
            lo, hi = clopper_pearson(k, n)
            rows.append(dict(assembler=a, depth=f"{dp}x", concordant=k, n=n, pct=round(100 * k / n, 1),
                             ci_low=round(100 * lo, 1), ci_high=round(100 * hi, 1)))
    prop = pd.DataFrame(rows)
    print(prop.to_string(index=False))
    prop.to_csv(os.path.join(out, "concordance_proportions.tsv"), sep="\t", index=False)

    wide = d.pivot_table(index="plasmid_id", columns=["assembler", "depth_x"], values="conc")
    tests = []
    for a in ["flye", "hybracter"]:
        for d1, d2 in zip(depths[:-1], depths[1:]):
            x, y = wide[(a, d1)], wide[(a, d2)]
            gained, lost = int(((x == 0) & (y == 1)).sum()), int(((x == 1) & (y == 0)).sum())
            tests.append(dict(comparison=f"{a} {d1}x vs {d2}x", gained=gained, lost=lost,
                              p_exact_mcnemar=mcnemar_exact(gained, lost)))
    for dp in depths:
        x, y = wide[("flye", dp)], wide[("hybracter", dp)]
        only_f, only_h = int(((x == 1) & (y == 0)).sum()), int(((x == 0) & (y == 1)).sum())
        tests.append(dict(comparison=f"flye vs hybracter {dp}x", gained=only_h, lost=only_f,
                          p_exact_mcnemar=mcnemar_exact(only_f, only_h)))
    tests = pd.DataFrame(tests)
    print("\nExact McNemar (gained/lost = concordant only at the second depth / only at the first;"
          " for assembler comparisons: concordant only with Hybracter / only with Flye)")
    print(tests.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    tests.to_csv(os.path.join(out, "concordance_mcnemar.tsv"), sep="\t", index=False)

    d["log2depth"] = np.log2(d["depth_x"])
    d["hybracter"] = (d["assembler"] == "hybracter").astype(int)
    gee_rows = []
    fam, cov = sm.families.Binomial(), sm.cov_struct.Exchangeable()
    for a in ["flye", "hybracter"]:
        s = d[d.assembler == a].sort_values(["plasmid_id", "depth_x"])
        m = smf.gee("conc ~ log2depth", "plasmid_id", s, family=fam, cov_struct=cov).fit()
        b, se, p = m.params["log2depth"], m.bse["log2depth"], m.pvalues["log2depth"]
        gee_rows.append(dict(model=a, term="log2(depth)", OR=np.exp(b), ci_low=np.exp(b - 1.96 * se),
                             ci_high=np.exp(b + 1.96 * se), p=p))
    s = d.sort_values(["plasmid_id", "assembler", "depth_x"])
    m = smf.gee("conc ~ log2depth + hybracter", "plasmid_id", s, family=fam, cov_struct=cov).fit()
    for term, lab in [("log2depth", "log2(depth)"), ("hybracter", "Hybracter vs Flye")]:
        b, se, p = m.params[term], m.bse[term], m.pvalues[term]
        gee_rows.append(dict(model="joint", term=lab, OR=np.exp(b), ci_low=np.exp(b - 1.96 * se),
                             ci_high=np.exp(b + 1.96 * se), p=p))
    gee = pd.DataFrame(gee_rows)
    print("\nLogistic GEE (exchangeable, clustered by plasmid, robust SE); OR per doubling of depth")
    print(gee.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    gee.to_csv(os.path.join(out, "concordance_gee.tsv"), sep="\t", index=False)


def read_vs_assembly(xlsx, out):
    """Exact McNemar from Supplementary Table S8a. Reads typed all isolates correctly at both
    depths, so every discordant pair is 'reads correct, assembly not'."""
    s = read_sheet(xlsx, "S8a_bart_carbapenemases", "Measure")
    row = s[s.iloc[:, 0].astype(str).str.startswith("All carbapenemase alleles")].iloc[0]
    print("\n== 2. Read-based vs assembly-based carbapenemase typing (S8a) ==")
    res = []
    for dp, rc, ac in [("10x", "Reads 10×", "Assembly 10×"), ("20x", "Reads 20×", "Assembly 20×")]:
        rk, n = map(int, str(row[rc]).split("/"))
        ak, _ = map(int, str(row[ac]).split("/"))
        # discordant counts follow from the margins only when one approach is correct in all isolates
        assert rk == n, "reads not correct in all isolates: per-isolate data needed (S8c and S2c)"
        b, c = n - ak, 0
        res.append(dict(depth=dp, reads=f"{rk}/{n}", assemblies=f"{ak}/{n}", reads_only=b, assemblies_only=c,
                        p_exact_mcnemar=mcnemar_exact(b, c)))
    res = pd.DataFrame(res)
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    res.to_csv(os.path.join(out, "reads_vs_assemblies_mcnemar.tsv"), sep="\t", index=False)


def actual_depth(xlsx, out):
    s = read_sheet(xlsx, "S6a_actual_depth", "Isolate")
    s = s[s["Cohort"] == "S1a"].copy()
    s["nominal"] = s["Nominal depth"].astype(str)
    print(f"\n== 3. Actual depth (S6a; {s['Isolate'].nunique()} carbapenemase-producing isolates) ==")
    rows = []
    for nom, g in s.groupby("nominal", sort=False):
        v = g["Actual depth (×)"].astype(float).to_numpy()
        q1, med, q3 = np.percentile(v, [25, 50, 75], method="weibull")
        rows.append(dict(nominal=nom, n=len(v), median=round(med, 1), iqr_low=round(q1, 1),
                         iqr_high=round(q3, 1), min=round(v.min(), 1), max=round(v.max(), 1)))
    tab = pd.DataFrame(rows)
    print(tab.to_string(index=False))
    r = s["Actual / nominal"].astype(float)
    print(f"Actual/nominal ratio: {r.min():.2f} to {r.max():.2f}")
    tab.to_csv(os.path.join(out, "actual_depth_summary.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    xlsx = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "statistics_output"
    os.makedirs(out, exist_ok=True)
    concordance(xlsx, out)
    read_vs_assembly(xlsx, out)
    actual_depth(xlsx, out)
