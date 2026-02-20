# workflow/rules/gt_conversion.smk
#
# Adds Ground-Truth (GT) ingestion + conversion + downsampling.
#
# Supports:
# - A/V RGB GT datasets (with optional {other_dir} wildcard in their source pattern)
# - vessel-only grayscale GT datasets
#
# Output layout:
# - A/V GT (no other_dir):
#     data/{dataset}/gt/raw
#     data/{dataset}/gt/converted
#     data/{dataset}/downsampled/{res}px/GTs
# - A/V GT (with other_dir):
#     data/{dataset}/{other_dir}/gt/raw
#     data/{dataset}/{other_dir}/gt/converted
#     data/{dataset}/{other_dir}/downsampled/{res}px/GTs
# - Vessel-only GT:
#     data/{dataset}/gt_vessel/raw
#     data/{dataset}/downsampled/{res}px/GTs_vessel
#
# Source locations are config-driven (see config.yaml: ground_truth section).
#
from pathlib import Path
from snakemake.io import directory
import glob

GT_CFG = config.get("ground_truth", {})
# Default GT root = legacy_root unless overridden
#GT_ROOT = Path(GT_CFG.get("root", config.get("legacy_root", "."))).resolve()

AV_CFG  = (GT_CFG.get("av_rgb", {}) or {}).get("datasets", {}) or {}
VES_CFG = (GT_CFG.get("vessel_gray", {}) or {}).get("datasets", {}) or {}

AV_OTHERDIR_DATASETS = {
    d for d, spec in AV_CFG.items()
    if "{other_dir}" in str(spec.get("pattern", ""))
}
AV_SIMPLE_DATASETS = set(AV_CFG.keys()) - AV_OTHERDIR_DATASETS

VES_OTHERDIR_DATASETS = {
    d for d, spec in VES_CFG.items()
    if "{other_dir}" in str(spec.get("pattern", ""))
}
VES_SIMPLE_DATASETS = set(VES_CFG.keys()) - VES_OTHERDIR_DATASETS



def _assert_no_other_dir_in_pattern(dataset: str, pattern: str):
    if "{other_dir}" in pattern:
        raise ValueError(
            f"Dataset '{dataset}' uses pattern with '{{other_dir}}', "
            "but you are calling a *simple* GT rule without other_dir."
        )

def _pattern_dir(pattern: str, **fmt) -> str:
    # Fill pattern enough to compute directory; sample placeholder replaced.
    fmt2 = dict(fmt)
    fmt2["sample"] = "__SAMPLE__"
    GT_ROOT = Path(GT_CFG.get("root", config.get("legacy_root", "."))).resolve()
    p = (GT_ROOT / pattern.format(**fmt2)).parent
    return str(p)

def _mapping(dataset: str):
    return AV_CFG[dataset].get("mapping", [])

# -------------------------
# A/V RGB GT (no other_dir)
# -------------------------

def av_src_dir_simple(wc):
    pat = AV_CFG[wc.dataset]["pattern"]
    _assert_no_other_dir_in_pattern(wc.dataset, pat)
    return _pattern_dir(pat, dataset=wc.dataset)

rule copy_gt_av_simple:
    wildcard_constraints:
        dataset="|".join(sorted(AV_SIMPLE_DATASETS)) if AV_SIMPLE_DATASETS else "NO_MATCH"
    input:
        src_dir = av_src_dir_simple
    output:
        out_dir = directory("data/{dataset}/gt/raw")
    shell:
        r"""
        mkdir -p {output.out_dir}
        cp -f {input.src_dir}/*.png {output.out_dir}/
        """

rule recolor_gt_av_simple:
    wildcard_constraints:
        dataset="|".join(sorted(AV_SIMPLE_DATASETS)) if AV_SIMPLE_DATASETS else "NO_MATCH"
    input:
        in_dir = "data/{dataset}/gt/raw"
    output:
        out_dir = directory("data/{dataset}/gt/converted")
    params:
        mapping = lambda wc: _mapping(wc.dataset),
        ext = ".png"
    script:
        "../scripts/recolor_rgb_dir_smk.py"

rule downsample_gt_av_simple:
    wildcard_constraints:
        dataset="|".join(sorted(AV_SIMPLE_DATASETS)) if AV_SIMPLE_DATASETS else "NO_MATCH"
    input:
        in_dir = "data/{dataset}/gt/converted"
    output:
        out_dir = directory("data/{dataset}/downsampled/{res}px/GTs")
    params:
        kind = "gt",
        width = lambda wc: int(wc.res),
        height = lambda wc: int(wc.res),
        ext = ".png"
    script:
        "../scripts/downsample_dir_smk.py"

# -------------------------
# A/V RGB GT (with other_dir)
# -------------------------

def av_src_dir_other(wc):
    pat = AV_CFG[wc.dataset]["pattern"]
    return _pattern_dir(pat, dataset=wc.dataset, other_dir=wc.other_dir)

rule copy_gt_av_otherdir:
    wildcard_constraints:
        dataset="|".join(sorted(AV_OTHERDIR_DATASETS)) if AV_OTHERDIR_DATASETS else "NO_MATCH"
    input:
        src_dir = av_src_dir_other
    output:
        out_dir = directory("data/{dataset}/{other_dir}/gt/raw")
    shell:
        r"""
        mkdir -p {output.out_dir}
        cp -f {input.src_dir}/*.png {output.out_dir}/
        """

rule recolor_gt_av_otherdir:
    wildcard_constraints:
        dataset="|".join(sorted(AV_OTHERDIR_DATASETS)) if AV_OTHERDIR_DATASETS else "NO_MATCH"
    input:
        in_dir = "data/{dataset}/{other_dir}/gt/raw"
    output:
        out_dir = directory("data/{dataset}/{other_dir}/gt/converted")
    params:
        mapping = lambda wc: _mapping(wc.dataset),
        ext = ".png"
    script:
        "../scripts/recolor_rgb_dir_smk.py"

rule downsample_gt_av_otherdir:
    wildcard_constraints:
        dataset="|".join(sorted(AV_OTHERDIR_DATASETS)) if AV_OTHERDIR_DATASETS else "NO_MATCH"
    input:
        in_dir = "data/{dataset}/{other_dir}/gt/converted"
    output:
        out_dir = directory("data/{dataset}/{other_dir}/downsampled/{res}px/GTs")
    params:
        kind = "gt",
        width = lambda wc: int(wc.res),
        height = lambda wc: int(wc.res),
        ext = ".png"
    script:
        "../scripts/downsample_dir_smk.py"

# -------------------------
# Vessel-only grayscale GT
# -------------------------

def vessel_src_dir(wc):
    pat = VES_CFG[wc.dataset]["pattern"]
    return _pattern_dir(pat, dataset=wc.dataset)

rule copy_gt_vessel:
    input:
        src_dir = vessel_src_dir
    output:
        out_dir = directory("data/{dataset}/gt_vessel/raw")
    shell:
        r"""
        mkdir -p {output.out_dir}
        cp -f {input.src_dir}/*.png {output.out_dir}/
        """

rule downsample_gt_vessel:
    input:
        in_dir = "data/{dataset}/gt_vessel/raw"
    output:
        out_dir = directory("data/{dataset}/downsampled/{res}px/GTs_vessel")
    params:
        kind = "gt_vessel",
        width = lambda wc: int(wc.res),
        ext = ".png"
    script:
        "../scripts/downsample_dir_smk.py"

def vessel_src_dir_other(wc):
    pat = VES_CFG[wc.dataset]["pattern"]
    return _pattern_dir(pat, dataset=wc.dataset, other_dir=wc.other_dir)

rule copy_gt_vessel_otherdir:
    input:
        src_dir = vessel_src_dir_other
    output:
        out_dir = directory("data/{dataset}/{other_dir}/gt_vessel/raw")
    shell:
        r"""
        mkdir -p {output.out_dir}
        cp -f {input.src_dir}/*.png {output.out_dir}/
        """

rule downsample_gt_vessel_otherdir:
    input:
        in_dir = "data/{dataset}/{other_dir}/gt_vessel/raw"
    output:
        out_dir = directory("data/{dataset}/{other_dir}/downsampled/{res}px/GTs_vessel")
    params:
        kind = "gt_vessel",
        width = lambda wc: int(wc.res),
        ext = ".png"
    script:
        "../scripts/downsample_dir_smk.py"

