#!/usr/bin/env python3
"""Step 14d: provenance of the nearest PLSDB plasmids (BioSample location, collection date, isolation source, host),
joined to each reference plasmid, with summary counts per similarity group.
Inputs : results/plsdb_nearest.tsv, results/plsdb_nearest_metadata.tsv (step 14c)
Output : results/plsdb_nearest_provenance.tsv; summary printed (also save with | tee results/plsdb_provenance_summary.txt)
Metadata are as recorded in BioSample via PLSDB 2023_11_03_v2 and not curated: the location may be the submitting
institution or culture collection rather than the place of isolation (e.g. FDA-ARGOS strains)."""
import csv, sys, collections
sys.path.insert(0, "lib")
from common import read_tsv, write_tsv  # noqa: E402
csv.field_size_limit(10**9)
MISSING = {"", "nan", "na", "n/a", "missing", "not collected", "not applicable", "unknown", "none"}
def val(s):
    s = (s or "").strip()
    return "" if s.lower() in MISSING else s
rd = csv.DictReader(open("results/plsdb_nearest_metadata.tsv"), delimiter="\t")
acc = next(c for c in rd.fieldnames if "ACC" in c.upper())
meta = {r[acc].split(".")[0]: r for r in rd}
cols = ["BIOSAMPLE_Location", "BIOSAMPLE_CollectionDate", "BIOSAMPLE_IsolationSource", "BIOSAMPLE_Host", "BIOSAMPLE_HostDisease", "NUCCORE_Description"]
rows = []
for r in read_tsv("results/plsdb_nearest.tsv"):
    m = meta.get(r["plsdb_hit"].split(".")[0], {})
    d, s = float(r["mash_distance"]), int(r["shared_hashes"].split("/")[0])
    grp = "near_identical (<=0.001, >=950/1000)" if d <= 0.001 and s >= 950 else ("<=0.01" if d <= 0.01 else ">0.01")
    out = {"plasmid_id": r["plasmid_id"], "ref_rep_type": r["ref_rep_type"], "plsdb_hit": r["plsdb_hit"], "mash_distance": r["mash_distance"],
           "shared_hashes": r["shared_hashes"], "similarity_group": grp, "species": m.get("TAXONOMY_species", "")}
    out.update({c: val(m.get(c, "")) for c in cols})
    loc = out["BIOSAMPLE_Location"]; out["country"] = loc.split(":")[0].strip() if loc else ""
    date = out["BIOSAMPLE_CollectionDate"]; out["year"] = date[:4] if date[:4].isdigit() else ""
    rows.append(out)
fields = list(rows[0].keys())
write_tsv("results/plsdb_nearest_provenance.tsv", rows, fields)
print(f"written: results/plsdb_nearest_provenance.tsv ({len(rows)} plasmids)")

def summary(sub, label):
    hits = {}
    for x in sub: hits[x["plsdb_hit"]] = x
    h = list(hits.values())
    yrs = [int(x["year"]) for x in h if x["year"]]
    print(f"\n== {label}: {len(sub)} plasmids, {len(h)} distinct nearest hits")
    print(f"  hits with location: {sum(1 for x in h if x['country'])}/{len(h)}; with collection year: {len(yrs)}/{len(h)}"
          + (f" (range {min(yrs)}-{max(yrs)})" if yrs else ""))
    print("  countries (distinct hits):", dict(collections.Counter(x["country"] or "NA" for x in h).most_common()))
    print("  countries (per plasmid):  ", dict(collections.Counter(x["country"] or "NA" for x in sub).most_common()))
    print("  host (per plasmid):       ", dict(collections.Counter(x["BIOSAMPLE_Host"] or "NA" for x in sub).most_common()))
    print("  species (per plasmid):    ", dict(collections.Counter(x["species"] for x in sub).most_common()))

summary(rows, "all")
for g in ["near_identical (<=0.001, >=950/1000)", "<=0.01", ">0.01"]:
    summary([x for x in rows if x["similarity_group"] == g], g)
print("\n== hits without host but with isolation source (non-human candidates):")
for x in {x["plsdb_hit"]: x for x in rows}.values():
    if not x["BIOSAMPLE_Host"] or x["BIOSAMPLE_Host"].lower() != "homo sapiens":
        print(f"  {x['plsdb_hit']}\t{x['species']}\t{x['country'] or 'NA'}\thost={x['BIOSAMPLE_Host'] or 'NA'}\tsource={x['BIOSAMPLE_IsolationSource'] or 'NA'}")
