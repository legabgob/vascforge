import numpy as np
from PIL import Image
from pathlib import Path

# Snakemake inputs
rgb_dir = Path(str(snakemake.input.rgb_dir))
out_dir = Path(str(snakemake.output.out_dir))

# Parameters
threshold = getattr(snakemake.params, "threshold", 10)
ext = getattr(snakemake.params, "ext", ".png")

out_dir.mkdir(parents=True, exist_ok=True)

rgb_files = sorted(rgb_dir.glob(f"*{ext}"))
if not rgb_files:
    raise SystemExit(f"No {ext} files found in {rgb_dir}")

count = 0
for rgb_path in rgb_files:
    # Load RGB image
    img = np.array(Image.open(rgb_path))
    
    # Compute std across RGB channels
    img_std = img.std(axis=2)
    
    # Create mask
    fovea_mask = img_std.copy()
    fovea_mask[fovea_mask > threshold] = 255
    
    # Center crop to square if not already square
    h, w = fovea_mask.shape
    if h != w:
        size = min(h, w)
        start_y = (h - size) // 2
        start_x = (w - size) // 2
        fovea_mask = fovea_mask[start_y:start_y+size, start_x:start_x+size]
    
    # Save mask
    out_path = out_dir / rgb_path.name
    Image.fromarray(fovea_mask.astype(np.uint8), mode="L").save(out_path)
    count += 1

print(f"Created {count} ROI masks from {rgb_dir} to {out_dir}")
