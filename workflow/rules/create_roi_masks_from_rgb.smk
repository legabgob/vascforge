from snakemake.io import directory
from pathlib import Path

def get_rgb_dirs(wildcards):
    """Find ALL RGB directories for a dataset (handles multiple splits)"""
    base = Path(config['legacy_root']) / wildcards.dataset
    
    # Try with middle directory (UKBB_1000_Fundus/fundus/seg_legacy/rgb)
    candidates = list(base.glob("*/seg_legacy/rgb"))
    if candidates:
        return [str(c) for c in candidates]
    
    # Try direct path (Fundus-AVSeg/seg_legacy/rgb)
    direct = base / "seg_legacy" / "rgb"
    if direct.exists():
        return [str(direct)]
    
    raise FileNotFoundError(f"Could not find RGB directory for {wildcards.dataset}")

rule make_roi_masks_from_rgb:
    input:
        rgb_dirs = get_rgb_dirs,
        segs_dir = "data/{dataset}/segs_converted"  # Use the merged directory!
    output:
        out_dir = directory("data/{dataset}/roi_masks")
    params:
        threshold = 5,
        ext = ".png"
    script:
        "../scripts/roi_from_rgb_smk.py"
