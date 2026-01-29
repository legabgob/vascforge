# workflow/rules/filter_downsampled_square.smk
"""
Filter downsampled directories to contain only square images.
This ensures that the refinement step only processes square images.
"""
from snakemake.io import directory

rule filter_downsampled_square:
    input:
        meta = "data/{dataset}/meta/meta_filtered.csv",
        segs = "data/{dataset}/downsampled/{res}px/segs_converted",
        masks = "data/{dataset}/downsampled/{res}px/roi_masks_binarized",
    output:
        segs = directory("data/{dataset}/downsampled/{res}px/segs_converted_square"),
        masks = directory("data/{dataset}/downsampled/{res}px/roi_masks_binarized_square"),
    log:
        "logs/filter_downsampled_square/{dataset}_{res}px.log"
    script:
        "../scripts/filter_nonsquare_images_smk.py"

