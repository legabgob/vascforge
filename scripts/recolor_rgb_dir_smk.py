# scripts/recolor_rgb_dir_smk.py
#
# Snakemake "directory -> directory" recoloring for RGB label images.
# Recolors by exact RGB matching using a dataset-specific mapping provided via params.mapping.
#
# Expected Snakemake interface:
#   input:
#     in_dir = "path/to/input_rgb_dir"
#   output:
#     out_dir = directory("path/to/output_rgb_dir")
#   params:
#     mapping = [ [[r,g,b],[r,g,b]], ... ]   # list of [srcRGB, dstRGB]
#     ext = ".png" (optional; default ".png")

from pathlib import Path
import numpy as np
from PIL import Image

def _validate_mapping(mapping):
    pairs = []
    for pair in mapping:
        if not (isinstance(pair, (list, tuple)) and len(pair) == 2):
            raise ValueError(f"Bad mapping entry (expected [src,dst]): {pair}")
        src, dst = pair
        if len(src) != 3 or len(dst) != 3:
            raise ValueError(f"RGB tuples must have length 3: {pair}")
        pairs.append((np.array(src, dtype=np.uint8), np.array(dst, dtype=np.uint8)))
    return pairs

def recolor(arr: np.ndarray, pairs):
    if arr.ndim != 3 or arr.shape[2] != 3:
        raise ValueError("Input must be HxWx3 RGB")
    if not pairs:
        return arr  # passthrough
    out = arr.copy()
    for src, dst in pairs:
        mask = (arr == src).all(axis=2)
        out[mask] = dst
    return out

in_dir = Path(str(snakemake.input.in_dir))
out_dir = Path(str(snakemake.output.out_dir))
ext = getattr(snakemake.params, "ext", ".png")
mapping = getattr(snakemake.params, "mapping", None)
if mapping is None:
    raise ValueError("Missing params.mapping")

pairs = _validate_mapping(mapping)

out_dir.mkdir(parents=True, exist_ok=True)
files = sorted(p for p in in_dir.glob(f"*{ext}") if p.is_file())
if not files:
    raise SystemExit(f"No {ext} files found in {in_dir}")

for i, src in enumerate(files, 1):
    arr = np.array(Image.open(src).convert("RGB"), dtype=np.uint8)
    out = recolor(arr, pairs)
    Image.fromarray(out, mode="RGB").save(out_dir / src.name)
    if i % 100 == 0:
        print(f"[recolor_rgb] {i}/{len(files)} done...")

print(f"[recolor_rgb] Done. Wrote {len(files)} file(s) to {out_dir}")

