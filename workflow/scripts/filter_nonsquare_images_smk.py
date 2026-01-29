"""
Filter out non-square images from meta.csv
"""
from pathlib import Path
import pandas as pd
from PIL import Image

# Get inputs
meta_in = Path(snakemake.input.meta)
images_dir = Path(snakemake.input.images_dir)
meta_out = Path(snakemake.output.meta)

# Read metadata
df = pd.read_csv(meta_in)
print(f"Original entries: {len(df)}")

# Check which images are square
square_mask = []
removed = []

for idx, row in df.iterrows():
    img_id = row['id']
    img_path = images_dir / f"{img_id}.png"
    
    if not img_path.exists():
        square_mask.append(False)
        removed.append(f"{img_id} (missing)")
        continue
    
    try:
        img = Image.open(img_path)
        is_square = img.size[0] == img.size[1]
        square_mask.append(is_square)
        
        if not is_square:
            removed.append(f"{img_id} ({img.size[0]}x{img.size[1]})")
    except Exception as e:
        square_mask.append(False)
        removed.append(f"{img_id} (error: {e})")

# Filter to square images only
df_filtered = df[square_mask].copy()

print(f"Square images: {len(df_filtered)}")
print(f"Removed {len(removed)} non-square images:")
for r in removed[:10]:  # Show first 10
    print(f"  - {r}")
if len(removed) > 10:
    print(f"  ... and {len(removed) - 10} more")

# Save filtered metadata
df_filtered.to_csv(meta_out, index=False)
