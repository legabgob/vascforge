"""
Filter out non-square images from meta.csv
"""
from pathlib import Path
import pandas as pd
from PIL import Image

# Get inputs
meta_in = Path(snakemake.input.meta)
segs_dir = Path(snakemake.input.segs_dir)
masks_dir = Path(snakemake.input.masks_dir)
meta_out = Path(snakemake.output.meta)

# Read metadata
df = pd.read_csv(meta_in)
print(f"Original entries: {len(df)}")

# Check which images are square
square_mask = []
removed = []

for idx, row in df.iterrows():
    img_id = row['id']
    seg_path = segs_dir / f"{img_id}.png"
    mask_path = masks_dir / f"{img_id}.png"
    
    # Check segmentation
    seg_square = False
    seg_size = None
    if seg_path.exists():
        try:
            img = Image.open(seg_path)
            seg_size = img.size
            seg_square = img.size[0] == img.size[1]
        except Exception as e:
            removed.append(f"{img_id} (seg error: {e})")
            square_mask.append(False)
            continue
    else:
        removed.append(f"{img_id} (seg missing)")
        square_mask.append(False)
        continue
    
    # Check mask
    mask_square = False
    mask_size = None
    if mask_path.exists():
        try:
            img = Image.open(mask_path)
            mask_size = img.size
            mask_square = img.size[0] == img.size[1]
        except Exception as e:
            removed.append(f"{img_id} (mask error: {e})")
            square_mask.append(False)
            continue
    else:
        removed.append(f"{img_id} (mask missing)")
        square_mask.append(False)
        continue
    
    # Only keep if BOTH are square
    if seg_square and mask_square:
        square_mask.append(True)
    else:
        square_mask.append(False)
        reasons = []
        if not seg_square:
            reasons.append(f"seg {seg_size[0]}x{seg_size[1]}")
        if not mask_square:
            reasons.append(f"mask {mask_size[0]}x{mask_size[1]}")
        removed.append(f"{img_id} ({', '.join(reasons)})")

# Filter to square images only
df_filtered = df[square_mask].copy()

print(f"Square images (both seg AND mask): {len(df_filtered)}")
print(f"Removed {len(removed)} entries:")
for r in removed[:20]:  # Show first 20
    print(f"  - {r}")
if len(removed) > 20:
    print(f"  ... and {len(removed) - 20} more")

# Save filtered metadata
df_filtered.to_csv(meta_out, index=False)

# Also write detailed log
with open(snakemake.log[0], 'w') as logfile:
    logfile.write(f"Original entries: {len(df)}\n")
    logfile.write(f"Square entries (both seg AND mask): {len(df_filtered)}\n")
    logfile.write(f"Removed: {len(removed)}\n\n")
    logfile.write("=== REMOVED ENTRIES ===\n")
    for r in removed:
        logfile.write(f"{r}\n")
