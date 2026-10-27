"""
FAERS ETL for dermatologic drug sex-difference pharmacovigilance study.

Reads quarterly zips directly from memory (no intermediate extraction),
keeps only dermatologic drug records, and writes cleaned tables for R.

Key methodological choices:
  - Dedup: for each caseid keep the record with the largest caseversion.
  - Drug role: primary suspect (PS) for main analysis.
  - Drug matching: word-boundary regex on drugname + prod_ai.
"""
import csv
import io
import os
import re
import sys
import traceback
import zipfile
from glob import glob

import numpy as np
import pandas as pd

ROOT = r"C:\Users\lenovo\WorkBuddy\2026-10-17-00-49-45"
RAW = os.path.join(ROOT, "data", "raw")
PROC = os.path.join(ROOT, "data", "processed")
os.makedirs(PROC, exist_ok=True)

CATALOG = os.path.join(ROOT, "code", "drug_catalog.csv")

READ_KW = dict(
    sep="$", dtype=str, encoding="latin-1",
    quoting=csv.QUOTE_NONE, on_bad_lines="skip", engine="c",
)

# ---------------------------------------------------------------- drug terms
cat = pd.read_csv(CATALOG)
cat["syn"] = cat["synonyms"].fillna("").str.lower().str.split("|")
syn2grp = {}
for _, r in cat.iterrows():
    for s in r["syn"]:
        s = s.strip()
        if s:
            syn2grp.setdefault(s, set()).add(r["drug_group"])

ALL_SYN = sorted(syn2grp, key=len, reverse=True)
BIG_RE = re.compile(r"\b(?:%s)\b" % "|".join(re.escape(s) for s in ALL_SYN), re.I)

GRP2CAT = dict(zip(cat["drug_group"], cat["category"]))
GRP2LBL = dict(zip(cat["drug_group"], cat["drug_label"]))

# per-group regex
GRP_RE = {}
for _, r in cat.iterrows():
    syns = [s.strip() for s in r["syn"] if s.strip()]
    GRP_RE[r["drug_group"]] = re.compile(
        r"\b(?:%s)\b" % "|".join(re.escape(s) for s in syns), re.I)


def read_tbl(zf, key):
    """Locate and read a table inside a quarterly zip."""
    cands = [n for n in zf.namelist()
             if n.lower().endswith(".txt") and key.lower() in os.path.basename(n).lower()]
    if not cands:
        return None
    with zf.open(cands[0]) as f:
        blob = f.read()
    # some quarters are UTF-8 with BOM; strip it or the first column name breaks
    if blob[:3] == b"\xef\xbb\xbf":
        blob = blob[3:]
    df = pd.read_csv(io.BytesIO(blob), **READ_KW)
    df.columns = [c.strip().lstrip("\ufeff").lower() for c in df.columns]
    return df


def norm(df, cols):
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    return df


def main():
    zips = sorted(glob(os.path.join(RAW, "faers_ascii_*.zip")))
    print(f"[info] {len(zips)} quarterly zips found", flush=True)

    drug_parts, demo_parts, reac_parts = [], [], []
    outc_parts, indi_parts = [], []
    seen = set()

    for zpath in zips:
        tag = os.path.basename(zpath).replace("faers_ascii_", "").replace(".zip", "")
        # ---- DRUG first: identify dermatologic reports
        try:
            with zipfile.ZipFile(zpath) as zf:
                d = read_tbl(zf, "DRUG")
                if d is None:
                    print(f"[warn] {tag}: no DRUG table", flush=True)
                    continue
                d = norm(d, ["primaryid", "caseid", "drug_seq", "role_cod",
                             "drugname", "prod_ai", "route"])
                # trim to suspect drugs before the (expensive) regex pass
                rc = d["role_cod"].astype(str).str.upper().str.strip()
                d = d[rc.isin(["PS", "SS"])]
                if d.empty:
                    print(f"[ok]   {tag}: 0 suspect drug rows", flush=True)
                    continue
                txt = (d["drugname"].fillna("") + " | " +
                       d["prod_ai"].fillna("").astype(str))
                hit = txt.str.contains(BIG_RE, na=False)
                d = d[hit]
                if d.empty:
                    print(f"[ok]   {tag}: 0 derm drug rows", flush=True)
                    continue

                rows = []
                sub = txt[hit]
                d = d.reset_index(drop=True)
                sub = sub.reset_index(drop=True)
                for grp, rx in GRP_RE.items():
                    m = sub.str.contains(rx, na=False)
                    if m.any():
                        t = d.loc[m, ["primaryid", "caseid", "role_cod", "route"]].copy()
                        t["drug_group"] = grp
                        rows.append(t)
                if not rows:
                    continue
                derm = pd.concat(rows, ignore_index=True)
                derm["category"] = derm["drug_group"].map(GRP2CAT)
                derm["quarter"] = tag
                drug_parts.append(derm)
                pids = set(derm["primaryid"])

                # ---- DEMO
                dem = read_tbl(zf, "DEMO")
                if dem is not None:
                    # FAERS renamed gndr_cod -> sex around 2014Q3 and uses
                    # occp_cod (not reporter_cod) for reporter occupation.
                    dem = norm(dem, ["primaryid", "caseid", "caseversion", "sex",
                                     "gndr_cod", "age", "age_cod", "event_dt",
                                     "fda_dt", "occp_cod", "occr_country"])
                    dem = dem[dem["primaryid"].isin(pids)]
                    if not dem.empty:
                        dem["quarter"] = tag
                        demo_parts.append(dem)

                        fin_pids = set(dem["primaryid"])
                        # ---- REAC
                        rc = read_tbl(zf, "REAC")
                        if rc is not None:
                            rc = norm(rc, ["primaryid", "caseid", "pt"])
                            rc = rc[rc["primaryid"].isin(fin_pids)]
                            if not rc.empty:
                                reac_parts.append(rc[["primaryid", "pt"]])
                        # ---- OUTC
                        oc = read_tbl(zf, "OUTC")
                        if oc is not None:
                            oc = norm(oc, ["primaryid", "caseid", "outc_cod",
                                           "outc_code"])
                            if "outc_code" in oc.columns:
                                oc["outc_cod"] = oc["outc_cod"].fillna(
                                    oc["outc_code"])
                            oc = oc[oc["primaryid"].isin(fin_pids)]
                            if not oc.empty:
                                outc_parts.append(oc[["primaryid", "outc_cod"]])
                        # ---- INDI
                        ind = read_tbl(zf, "INDI")
                        if ind is not None:
                            ind = norm(ind, ["primaryid", "caseid", "indi_pt"])
                            ind = ind[ind["primaryid"].isin(fin_pids)]
                            if not ind.empty:
                                indi_parts.append(ind[["primaryid", "indi_pt"]])
            print(f"[ok]   {tag}: derm drug rows={len(derm):,}", flush=True)
            seen.add(tag)
        except Exception as e:
            print(f"[ERR]  {tag}: {e}", flush=True)

    print("\n[info] concatenating ...", flush=True)
    drug = pd.concat(drug_parts, ignore_index=True)
    demo = pd.concat(demo_parts, ignore_index=True)
    reac = pd.concat(reac_parts, ignore_index=True).drop_duplicates()
    outc = pd.concat(outc_parts, ignore_index=True).drop_duplicates()
    indi = pd.concat(indi_parts, ignore_index=True).drop_duplicates()

    print(f"raw derm drug rows : {len(drug):,}")
    print(f"raw demo rows      : {len(demo):,}")
    print(f"raw reac rows      : {len(reac):,}")

    # ---------------------------------------------------------- dedup by caseid
    print("[info] deduplicating DEMO (caseid -> max caseversion) ...", flush=True)
    demo["caseversion"] = pd.to_numeric(demo["caseversion"], errors="coerce").fillna(-1)
    demo["_fda"] = demo["fda_dt"].fillna("").astype(str)
    demo = demo.sort_values(["caseid", "caseversion", "_fda"])
    demo = demo.groupby("caseid", as_index=False, sort=False).tail(1)
    demo = demo.drop(columns=["_fda"])
    print(f"after dedup        : {len(demo):,}")

    # ---------------------------------------------------------- sex / age
    sx = demo["sex"].astype(str).str.strip().str.upper()
    gx = demo["gndr_cod"].astype(str).str.strip().str.upper()
    demo["sex"] = sx.where(sx.isin(["M", "F"]), gx)
    n_sex_ok = int(demo["sex"].isin(["M", "F"]).sum())
    print(f"sex resolved       : {n_sex_ok:,} / {len(demo):,} "
          f"({100*n_sex_ok/max(len(demo),1):.1f}%)")
    demo = demo[demo["sex"].isin(["M", "F"])].copy()
    print(f"after sex filter   : {len(demo):,}")

    def to_years(a, c):
        a = pd.to_numeric(a, errors="coerce")
        c = c.astype(str).str.strip().str.upper()
        yr = a.where(c == "YR")
        yr = yr.fillna(a.where(c == "MON") / 12)
        yr = yr.fillna(a.where(c == "WK") / 52.0)
        yr = yr.fillna(a.where(c == "DY") / 365.0)
        yr = yr.fillna(a.where(c == "DEC") * 10)
        return yr

    demo["age_yr"] = to_years(demo["age"], demo["age_cod"])
    demo.loc[(demo["age_yr"] < 0) | (demo["age_yr"] > 120), "age_yr"] = np.nan
    demo["event_year"] = pd.to_numeric(
        demo["event_dt"].astype(str).str.slice(0, 4), errors="coerce")
    demo.loc[demo["event_year"].isna(), "event_year"] = pd.to_numeric(
        demo["fda_dt"].astype(str).str.slice(0, 4), errors="coerce")

    keep_pids = set(demo["primaryid"])
    drug = drug[drug["primaryid"].isin(keep_pids)].copy()
    reac = reac[reac["primaryid"].isin(keep_pids)].copy()
    outc = outc[outc["primaryid"].isin(keep_pids)].copy()
    indi = indi[indi["primaryid"].isin(keep_pids)].copy()

    print(f"drug rows (final)  : {len(drug):,}")
    print(f"reac rows (final)  : {len(reac):,}")
    print(f"unique cases       : {len(demo):,}")

    # ---------------------------------------------------------- save
    cols_demo = ["primaryid", "caseid", "sex", "age_yr", "event_year",
                 "occp_cod", "occr_country", "quarter"]
    demo[[c for c in cols_demo if c in demo.columns]].to_csv(
        os.path.join(PROC, "demo.csv.gz"), index=False, compression="gzip")
    drug[["primaryid", "caseid", "drug_group", "category", "role_cod", "route", "quarter"]].to_csv(
        os.path.join(PROC, "drug.csv.gz"), index=False, compression="gzip")
    reac.to_csv(os.path.join(PROC, "reac.csv.gz"), index=False, compression="gzip")
    outc.to_csv(os.path.join(PROC, "outc.csv.gz"), index=False, compression="gzip")
    indi.to_csv(os.path.join(PROC, "indi.csv.gz"), index=False, compression="gzip")

    print("\n[done] wrote to", PROC)
    print("\n--- drug_group counts (PS only) ---")
    ps = drug[drug["role_cod"].astype(str).str.upper() == "PS"]
    print(ps.groupby(["category", "drug_group"]).size().sort_values(
        ascending=False).head(30).to_string())
    print("\n--- sex distribution ---")
    print(demo["sex"].value_counts().to_string())


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
