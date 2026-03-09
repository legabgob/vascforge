# workflow/rules/vascxgray_to_rgb.smk
import os
from pathlib import Path
from snakemake.io import directory


def find_av_dir(wc):
    """
    Find a seg_legacy/av directory for this dataset, supporting:
      A) legacy_root/dataset/seg_legacy/av/
      B) legacy_root/dataset/*/seg_legacy/av/
    Deterministic choice: prefer A, else first sorted match in B.
    """
    LEGACY_ROOT = Path(config["legacy_root"]).resolve()
    droot = LEGACY_ROOT / wc.dataset

    cand_a = droot / "seg_legacy" / "av"
    if cand_a.is_dir():
        return str(cand_a)

    cands_b = sorted(droot.glob("*/seg_legacy/av"))
    for c in cands_b:
        if c.is_dir():
            return str(c)

    raise FileNotFoundError(f"No seg_legacy/av folder found under {droot}")

rule gray_to_rgb:
    """
    Convert grayscale AV label images (0,1,2,3) to RGB A/V/BV mapping.
    Runs in batch on a dataset directory.
    """
    input:
        av_dir = find_av_dir
    output:
        out_dir = directory("data/{dataset}/segs_converted")
    benchmark:
        "benchmarks/gray_to_rgb/{dataset}.tsv"
    params:
        name = "gray_to_rgb",
        time = "00:30:00",
        mem = 4000,
        threads = 1,
    script:
        "../scripts/gray_to_rgb_dir_smk.py"

