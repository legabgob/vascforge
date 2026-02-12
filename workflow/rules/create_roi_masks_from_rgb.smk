from snakemake.io import directory
from pathlib import Path

rule make_roi_masks_from_rgb:
    input:
        rgb_dir = lambda wc: str(list(Path(config['legacy_root']).glob(f"{wc.dataset}/*/seg_legacy/rgb"))[0])
    output:
        out_dir = directory("data/{dataset}/roi_masks")
    params:
        threshold = 5,
        ext = ".png"
    script:
        "../scripts/roi_from_rgb_smk.py"

