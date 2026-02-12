from snakemake.io import directory
from pathlib import Path

def get_rgb_dir(wildcards):
    base = Path(config['legacy_root']) / wildcards.dataset
    
    # Try with middle directory (UKBB_1000_Fundus/fundus/seg_legacy/rgb)
    candidates = list(base.glob("*/seg_legacy/rgb"))
    if candidates:
        return str(candidates[0])
    
    # Try direct path (Fundus-AVSeg/seg_legacy/rgb)
    direct = base / "seg_legacy" / "rgb"
    if direct.exists():
        return str(direct)
    
    raise FileNotFoundError(f"Could not find RGB directory for {wildcards.dataset}")

rule make_roi_masks_from_rgb:
    input:
        rgb_dir = get_rgb_dir
    output:
        out_dir = directory("data/{dataset}/roi_masks")
    params:
        threshold = 10,
        ext = ".png"
    script:
        "../scripts/roi_from_rgb_smk.py"
