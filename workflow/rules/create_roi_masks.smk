def find_roi_csv(wc):
    """
    Automatically selects bounds.csv or meta.csv if present.
    Priority: bounds.csv > meta.csv
    """
    base = f"data/{wc.dataset}/meta"
    bounds = f"{base}/bounds.csv"
    meta   = f"{base}/meta.csv"

    if os.path.exists(bounds):
        return bounds
    if os.path.exists(meta):
        return meta

    raise FileNotFoundError(
        f"No bounds.csv or meta.csv found in {base}/"
    )


rule make_roi_masks:
    """
    Generate ROI circular masks from either bounds.csv or meta.csv,
    depending on which one exists.
    """
    input:
        csv = find_roi_csv
    output:
        out_dir = directory("data/{dataset}/roi_masks")
    params:
        img_dir = lambda wc: f"data/{wc.dataset}/images",
        img_ext = ".png",
        skip_if_flag_false = True,
    script:
        "scripts/roi_from_csv_smk.py"
