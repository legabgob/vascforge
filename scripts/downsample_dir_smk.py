# scripts/downsample_dir_smk.py
#
# Snakemake "directory -> directory" batch downsampling.
# Keeps discrete labels crisp by using NEAREST resampling.
#
# Expected Snakemake interface:
#   input:
#     in_dir = "path/to/input_dir"
#   output:
#     out_dir = directory("path/to/output_dir")
#   params:
#     kind  = "segs_converted" | "roi_masks_binarized" | "gt" | "gt_vessel" | ...
#     width = 576 | 1024 | ...
#     ext   = ".png" (optional; default ".png")
#
# This script processes files directly under in_dir (non-recursive).

from pathlib import Path
from PIL import Image

RGB_KINDS = {
    "segs_converted", "segs", "predictions",
    "gt", "gts", "gt_av", "gts_av", "gt_converted", "gts_converted",
}
GRAY_KINDS = {
    "roi_masks_binarized", "roi_masks", "masks", "mask",
    "gt_vessel", "gts_vessel", "vessel_gt", "vessel", "gt_gray",
}

def normalize_kind(kind: str) -> str:
    return kind.strip().lower().replace("-", "_")

def downsample_one(src: Path, dst: Path, kind: str, target_w: int):
    kind_n = normalize_kind(kind)

    im = Image.open(src)
    if kind_n in RGB_KINDS:
        im = im.convert("RGB")
        resample = Image.NEAREST
    elif kind_n in GRAY_KINDS:
        im = im.convert("L")
        resample = Image.NEAREST
    else:
        raise ValueError(
            f"Unknown kind='{kind}'.\n"
            f"RGB kinds: {sorted(RGB_KINDS)}\n"
            f"GRAY kinds: {sorted(GRAY_KINDS)}"
        )

    w, h = im.size
    if w == 0 or h == 0:
        raise ValueError(f"Invalid image size for {src}: {(w,h)}")

    target_h = int(round(h * (target_w / float(w))))
    im_small = im.resize((target_w, target_h), resample=resample)

    dst.parent.mkdir(parents=True, exist_ok=True)
    im_small.save(dst)

# -------- Snakemake entrypoint --------
in_dir = Path(str(snakemake.input.in_dir))
out_dir = Path(str(snakemake.output.out_dir))
kind = str(snakemake.params.kind)
width = int(snakemake.params.width)
ext = getattr(snakemake.params, "ext", ".png")

out_dir.mkdir(parents=True, exist_ok=True)

files = sorted(p for p in in_dir.glob(f"*{ext}") if p.is_file())
if not files:
    raise SystemExit(f"No {ext} files found in {in_dir}")

for i, src in enumerate(files, 1):
    dst = out_dir / src.name
    downsample_one(src, dst, kind, width)
    if i % 100 == 0:
        print(f"[downsample_dir] {i}/{len(files)} done...")

print(f"[downsample_dir] Done. Wrote {len(files)} file(s) to {out_dir}")

