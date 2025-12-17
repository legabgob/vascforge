import os
from snakemake.io import directory

def find_roi_csv(wc):
    base = f"data/{wc.dataset}/meta"
    bounds = f"{base}/bounds.csv"
    meta   = f"{base}/meta.csv"

    if os.path.exists(bounds):
        return bounds
    if os.path.exists(meta):
        return meta

    raise FileNotFoundError(f"No bounds.csv or meta.csv in {base}/")


rule make_roi_masks:
    input:
        meta_done = "data/{dataset}/meta/.copied",
        csv = find_roi_csv,
    output:
        out_dir = directory("data/{dataset}/roi_masks")
    params:
        img_dir = lambda wc: f"data/{wc.dataset}/images",
        img_ext = ".png",
        skip_if_flag_false = True,
    script:
        "scripts/roi_from_csv_smk.py"
