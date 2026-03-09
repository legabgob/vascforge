# workflow/rules/copy_meta.smk
import os
from pathlib import Path

def find_meta_like_csv(wc):
    """
    Find bounds.csv (preferred) or meta.csv under:
      A) legacy_root/dataset/seg_legacy/
      B) legacy_root/dataset/*/seg_legacy/
    """
    LEGACY_ROOT = Path(config["legacy_root"]).resolve()
    droot = LEGACY_ROOT / wc.dataset

    # Prefer A
    base_a = droot / "seg_legacy"
    bounds_a = base_a / "bounds.csv"
    meta_a = base_a / "meta.csv"
    if bounds_a.exists():
        return str(bounds_a)
    if meta_a.exists():
        return str(meta_a)

    # Else search B
    seg_legacy_dirs = sorted(droot.glob("*/seg_legacy"))
    for sd in seg_legacy_dirs:
        bounds = sd / "bounds.csv"
        meta = sd / "meta.csv"
        if bounds.exists():
            return str(bounds)
        if meta.exists():
            return str(meta)

    raise FileNotFoundError(f"No meta.csv or bounds.csv found for dataset {wc.dataset} under {droot}")

rule copy_meta_csv:
    """
    Copy bounds.csv or meta.csv from legacy seg_legacy into workspace.
    Always written as data/{dataset}/meta/meta.csv (content may be either format).
    """
    input:
        src = find_meta_like_csv
    output:
        dst = "data/{dataset}/meta/meta.csv"
    benchmark:
        "benchmarks/copy_meta_csv/{dataset}.tsv"
    params:
        name = "copy_meta_csv",
        time = "00:10:00",
        mem = 1000,
        threads = 1,
    shell:
        r"""
        mkdir -p $(dirname "{output.dst}")
        cp "{input.src}" "{output.dst}"
        """

