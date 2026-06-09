import polars as pl
import pathlib
import re

# Base directory
base_dir = pathlib.Path("qc/qualimap")

# Regex patterns for extracting metrics
patterns = {
    "Mean mapping quality (Q)": r"mean mapping quality\s*=\s*([\d.]+)",
    "Duplication Rate (%)": r"duplication rate\s*=\s*([\d.]+)%",
    "Median insert size (bp)": r"median insert size\s*=\s*([\d.]+)",
    "Mean coverage": r"mean coverageData\s*=\s*([\d.]+)",
    "Mapped Reads (%)": r"number of mapped reads\s*=\s*[\d,]+\s*\(([\d.]+)%",
}

rows = []

for folder in base_dir.iterdir():
    # Skip hidden files/folders, those with lowercase "s"/"p", and ".done" dirs
    if not folder.is_dir() or "s" in folder.name or "p" in folder.name or folder.name.endswith(".done"):
        continue

    genome_file = folder / "genome_results.txt"
    if not genome_file.exists():
        continue

    sample_id = folder.name
    text = genome_file.read_text()

    values = {"sample_id": sample_id}
    for field, pat in patterns.items():
        m = re.search(pat, text)
        if m:
            values[field] = float(m.group(1).replace(",", ""))
        else:
            values[field] = None

    rows.append(values)

# Build Polars DataFrame
df = pl.DataFrame(rows)

# Round all numeric columns to whole numbers
numeric_cols = [c for c in df.columns if c != "sample_id"]
df = df.with_columns([pl.col(c).round(0).cast(int) for c in numeric_cols if c in df.columns])

# Add "x" suffix to coverage
df = df.with_columns(
    (pl.col("Mean coverage").cast(str) + "x").alias("Mean coverage")
)

# Sort alphabetically
df = df.sort("sample_id")

# Save and show
df.write_csv("qualimap_summary.csv")