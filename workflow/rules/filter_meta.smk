# workflow/rules/filter_meta.smk
"""
Filter meta.csv to only include images with corresponding segmentations.
Prevents downstream issues with orphan ROI masks.
"""

rule filter_meta_by_segs:
    """
    Filter meta.csv to only include rows where segmentation images exist.
    This prevents creating ROI masks for non-existent images.
    """
    input:
        segs_dir = "data/{dataset}/segs_converted",
        meta_csv = "data/{dataset}/meta/meta.csv",
    output:
        meta_filtered = "data/{dataset}/meta/meta_filtered.csv",
    params:
        ext = ".png",
    script:
        "../scripts/filter_meta_by_segs_smk.py"


rule filter_nonsquare_images:
    input:
        meta = "data/{dataset}/meta/meta_filtered.csv",
        segs_dir = "data/{dataset}/segs_converted",
        masks_dir = "data/{dataset}/roi_masks_binarized"
    output:
        meta = "data/{dataset}/meta/meta_filtered_square.csv"
    log:
        "logs/filter_nonsquare/{dataset}.log"
    script:
        "../scripts/filter_nonsquare_images_smk.py"
