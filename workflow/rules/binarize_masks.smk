# workflow/rules/binarize_masks.smk
from snakemake.io import directory

rule binarize_roi_masks:
    """
    Replace pixel value 1 -> 255 for ROI masks (directory batch).
    """
    input:
        in_dir = "data/{dataset}/roi_masks"
    output:
        out_dir = directory("data/{dataset}/roi_masks_binarized")
    benchmark:
        "benchmarks/binarize_roi_masks/{dataset}.tsv"
    params:
        name = "binarize_roi_masks",
        time = "00:30:00",
        mem = 4000,
        threads = 1,
        ext = ".png"
    script:

