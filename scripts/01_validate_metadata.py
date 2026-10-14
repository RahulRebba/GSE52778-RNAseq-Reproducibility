#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import sys

ROOT = Path(__file__).resolve().parents[1]
meta_path = ROOT / "data" / "metadata" / "sample_metadata_cleaned.csv"
out_dir = ROOT / "results" / "generated"
out_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(meta_path)
required = ["sra_run","geo_sample","donor","treatment","dex","study","include_primary_analysis"]

problems = []
missing_cols = [c for c in required if c not in df.columns]
if missing_cols:
    problems.append(f"Missing columns: {missing_cols}")

if not missing_cols:
    if df[required].isna().any().any():
        problems.append("Missing values found in required fields.")
    if df["sra_run"].duplicated().any():
        problems.append("Duplicate SRA run IDs found.")
    if df["geo_sample"].duplicated().any():
        problems.append("Duplicate GEO sample IDs found.")
    if not set(df["dex"]).issubset({"untrt","trt"}):
        problems.append("Unexpected dex values.")
    # Every donor should have one untreated and one treated sample
    pairs = df.groupby("donor")["dex"].apply(lambda x: set(x))
    bad_pairs = [donor for donor, vals in pairs.items() if vals != {"untrt","trt"}]
    if bad_pairs:
        problems.append(f"Incomplete donor pairs: {bad_pairs}")
    if len(df) != 8:
        problems.append(f"Expected 8 primary samples; found {len(df)}.")

report = []
report.append("GSE52778 METADATA VALIDATION")
report.append("=" * 34)
report.append(f"Rows: {len(df)}")
report.append(f"Unique SRA runs: {df['sra_run'].nunique() if 'sra_run' in df else 'NA'}")
report.append(f"Unique GEO samples: {df['geo_sample'].nunique() if 'geo_sample' in df else 'NA'}")
report.append(f"Donors: {df['donor'].nunique() if 'donor' in df else 'NA'}")
report.append(f"Treatment counts: {df['dex'].value_counts().to_dict() if 'dex' in df else 'NA'}")
report.append("")
if problems:
    report.append("STATUS: INVESTIGATE")
    report.extend([f"- {p}" for p in problems])
else:
    report.append("STATUS: PASS")
    report.append("- 8 unique samples")
    report.append("- 4 complete untreated/treated donor pairs")
    report.append("- No missing required metadata")
    report.append("- No duplicate sample IDs")

text = "\n".join(report) + "\n"
print(text)
(out_dir / "metadata_validation_report.txt").write_text(text, encoding="utf-8")
sys.exit(1 if problems else 0)
