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
