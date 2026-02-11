# workflow/rules/create_roi_masks_from_rgb.smk
from snakemake.io import directory

rule make_roi_masks_from_rgb:
    input:
        rgb_dir = lambda wc: glob(f"{config['legacy_root']}/{wc.dataset}/*/seg_legacy/rgb")[0]
    output:
        out_dir = directory("data/{dataset}/roi_masks")
    params:
        threshold = 10,  # std threshold for background detection
        ext = ".png"
    script:
        "../scripts/roi_from_rgb_smk.py"
