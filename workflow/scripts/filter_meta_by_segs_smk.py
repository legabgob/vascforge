# scripts/filter_meta_by_segs_smk.py
"""
Filter meta.csv to only include images that have corresponding segmentations.
This prevents downstream issues where ROI masks are created for non-existent images.
"""

from pathlib import Path
import pandas as pd

# Snakemake inputs
segs_dir = Path(str(snakemake.input.segs_dir))
meta_in = Path(str(snakemake.input.meta_csv))

# Snakemake output
meta_out = Path(str(snakemake.output.meta_filtered))

# Extension for segmentation files
ext = str(getattr(snakemake.params, "ext", ".png"))

# Get list of actual segmentation filenames (stems only, no extension)
seg_files = {p.stem for p in segs_dir.glob(f"*{ext}") if p.is_file()}

if not seg_files:
    raise SystemExit(f"No segmentation files found in {segs_dir}")

print(f"Found {len(seg_files)} segmentation files in {segs_dir}")

# Read meta.csv
df = pd.read_csv(meta_in)
n_original = len(df)

# Determine which column contains the image names
# Try common patterns: first column, or columns named 'id', 'name', 'image', 'filename'
id_col = None
for col in ['id', 'name', 'image', 'filename']:
    if col in df.columns:
        id_col = col
        break

if id_col is None:
    # Assume first column contains image IDs
    id_col = df.columns[0]

print(f"Using column '{id_col}' as image identifier")

# Filter to only include rows where the image exists in segmentations
# Handle both with and without extension in meta.csv
df_filtered = df[df[id_col].apply(
    lambda x: str(x).replace(ext, '') in seg_files
)]

n_filtered = len(df_filtered)
n_removed = n_original - n_filtered

# Save filtered meta.csv
meta_out.parent.mkdir(parents=True, exist_ok=True)
df_filtered.to_csv(meta_out, index=False)

print(f"Filtered meta.csv: {n_original} → {n_filtered} rows ({n_removed} removed)")

if n_removed > 0:
    removed_ids = set(df[id_col]) - set(df_filtered[id_col])
    print(f"Removed entries: {sorted(removed_ids)}")
