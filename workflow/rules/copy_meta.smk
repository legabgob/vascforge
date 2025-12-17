# rules/copy_meta.smk

import os
from pathlib import Path
from snakemake.io import glob_wildcards

# Expect config.yaml to define the single root
SEGROOT = Path(config["legacy_root"]).resolve()

# Discover all dataset bases that have seg_legacy/meta.csv or seg_legacy/bounds.csv
# Case A: ROOT/{dataset_base}/seg_legacy/(meta|bounds).csv
bases_a_meta, = glob_wildcards(str(SEGROOT / "{dataset_base}/seg_legacy/meta.csv"))
bases_a_bounds, = glob_wildcards(str(SEGROOT / "{dataset_base}/seg_legacy/bounds.csv"))

# Case B: ROOT/{dataset}/{other_dir}/seg_legacy/(meta|bounds).csv  -> dataset_base = "{dataset}/{other_dir}"
bases_b_meta, = glob_wildcards(str(SEGROOT / "{dataset}/{other_dir}/seg_legacy/meta.csv"))
bases_b_bounds, = glob_wildcards(str(SEGROOT / "{dataset}/{other_dir}/seg_legacy/bounds.csv"))

DATASET_BASES = sorted(set(bases_a_meta) | set(bases_a_bounds) |
                       set([f"{d}/{o}" for d, o in zip(*glob_wildcards(str(SEGROOT / "{dataset}/{other_dir}/seg_legacy/meta.csv"))[:2])]) |
                       set([f"{d}/{o}" for d, o in zip(*glob_wildcards(str(SEGROOT / "{dataset}/{other_dir}/seg_legacy/bounds.csv"))[:2])]))

# The two zip(*glob_wildcards(...)) lines above are a bit ugly; a safer approach is below:
# We'll just re-scan explicitly in a helper, and also pick bounds over meta if both exist.

def find_meta_like_csv(wc):
    """
    Return path to bounds.csv (preferred) or meta.csv for this dataset_base.
    dataset_base is relative under SEGROOT, e.g. 'Fundus-AVSeg' or 'FIVES/train'.
    """
    base = SEGROOT / wc.dataset_base / "seg_legacy"
    bounds = base / "bounds.csv"
    meta = base / "meta.csv"
    if bounds.exists():
        return str(bounds)
    if meta.exists():
        return str(meta)
    raise FileNotFoundError(f"No meta.csv or bounds.csv under {base}")

rule copy_meta_csv:
    """
    Copy bounds.csv or meta.csv from seg_legacy into the pipeline workspace.
    Output is always named meta.csv (the ROI script can handle both formats).
    """
    input:
        src=find_meta_like_csv
    output:
        dst="data/{dataset_base}/meta/meta.csv"
    shell:
        r"""
        mkdir -p $(dirname "{output.dst}")
        cp "{input.src}" "{output.dst}"
        """

