# workflow/rules/create_roi_masks.smk
from snakemake.io import directory

rule make_roi_masks:
    input:
        csv = "data/{dataset}/meta/meta_filtered.csv"
    output:
        out_dir = directory("data/{dataset}/roi_masks")
    params:
        img_dir = lambda wc: f"data/{wc.dataset}/images",
        img_ext = ".png",
        skip_if_flag_false = True,
    script:
        "../scripts/roi_from_csv_smk.py"

