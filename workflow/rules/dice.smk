# workflow/rules/dice.smk
#
# Runs scripts/compute_metrics_smk.py to produce DICE metric CSVs.
#
# Supports:
#   - A/V GT datasets (RGB GTs) in two layouts:
#       * simple:   data/{dataset}/downsampled/{res}px/GTs
#       * otherdir: data/{dataset}/{other_dir}/downsampled/{res}px/GTs
#   - vessel-only GT datasets (optional; requires metrics.vessel_gt_dirs entries)
#
# Why two A/V rules?
#   Datasets like "leuven-haifa" have multiple GT subfolders ({other_dir}) that must not be merged.
#   Therefore metrics are computed per (dataset, other_dir) and written under:
#       results/metrics/{dataset}/{other_dir}/metrics_{res}.csv

import re
from pathlib import Path
from snakemake.io import glob_wildcards

# --------------------------
# Configuration (safe defaults)
# --------------------------

METRICS_CFG = config.get("metrics", {}) if isinstance(config, dict) else {}

_default_av = {"Fundus-AVSeg", "leuven-haifa"}
_datasets = config.get("datasets", [])
if isinstance(_datasets, dict):
    _datasets = list(_datasets.keys())
_datasets = set(_datasets)

AV_GT_DATASETS = set(METRICS_CFG.get("av_gt_datasets", list(_default_av & _datasets)))

# Which datasets to compute metrics for
METRICS_DATASETS = METRICS_CFG.get("datasets", sorted(AV_GT_DATASETS))

# Output root
METRICS_DIR = METRICS_CFG.get("out_dir", "results/metrics")

# Optional per-dataset overrides for GT locations (A/V GT).
# Only applies to SIMPLE layout (dataset-level). For other_dir layouts we rely on pipeline outputs.
AV_GT_DIRS = METRICS_CFG.get("av_gt_dirs", {})

# Optional vessel-only native dirs mapping (only used if you explicitly include such datasets in METRICS_DATASETS)
VESSEL_GT_DIRS = METRICS_CFG.get("vessel_gt_dirs", {})

# Refinement grid (only used to declare dependency inputs)
RESOLUTIONS = [str(r) for r in config.get("resolutions", ["576", "1024"])]
k_start, k_end = config.get("k_range", [3, 9])
K_VALUES = list(range(int(k_start), int(k_end)))

# --------------------------
# Discover which A/V GT datasets use {other_dir} (from ground_truth config)
# --------------------------

GT_CFG = config.get("ground_truth", {}) if isinstance(config, dict) else {}
GT_ENABLED = bool(GT_CFG.get("enabled", False))
#GT_ROOT = Path(GT_CFG.get("root", config.get("legacy_root", "."))).resolve()

AV_CFG = ((GT_CFG.get("av_rgb", {}) or {}).get("datasets", {}) or {}) if GT_ENABLED else {}

AV_OTHERDIR_DATASETS = {
    d for d, spec in AV_CFG.items()
    if "{other_dir}" in str((spec or {}).get("pattern", ""))
}

# Restrict to datasets we actually compute metrics for
METRICS_AV_OTHERDIR = sorted(
    [d for d in METRICS_DATASETS if (d in AV_GT_DATASETS and d in AV_OTHERDIR_DATASETS)]
)
METRICS_AV_SIMPLE = sorted(
    [d for d in METRICS_DATASETS if (d in AV_GT_DATASETS and d not in AV_OTHERDIR_DATASETS)]
)
METRICS_VESSEL_ONLY = sorted([d for d in METRICS_DATASETS if d not in AV_GT_DATASETS])

# For other_dir A/V datasets, discover available other_dir values by globbing the configured GT pattern.
METRICS_AV_OTHERDIR_VALUES = {}  # dataset -> [other_dir,...]

# Get exclusion list from config
METRICS_EXCLUDE_SPLITS = METRICS_CFG.get("exclude_splits", {})

for d in METRICS_AV_OTHERDIR:
    spec = AV_CFG.get(d, {}) or {}
    pat = spec.get("pattern")
    if not pat:
        raise ValueError(
            f"metrics: dataset '{d}' requires ground_truth.av_rgb.datasets.{d}.pattern with '{{other_dir}}'"
        )
    GT_ROOT = Path(GT_CFG.get("root", config.get("legacy_root", "."))).resolve()
    abs_pat = str(GT_ROOT / pat)
    other_dirs, _samples = glob_wildcards(abs_pat)  # expects (other_dir, sample)
    
    # Apply exclusions from config
    excluded = METRICS_EXCLUDE_SPLITS.get(d, [])
    if excluded:
        print(f"Excluding splits for {d} (from config): {excluded}")
    
    # Filter to only include splits that have overlap with segmentations
    # For datasets with split-specific GTs but dataset-level segs, check for file overlap
    valid_other_dirs = []
    for od in set(other_dirs):
        # Skip if explicitly excluded in config
        if od in excluded:
            continue
            
        has_overlap = False
        for res in RESOLUTIONS:
            # Try split-specific segs first
            seg_dir = Path(f"data/{d}/{od}/downsampled/{res}px/segs_converted_square")
            if not seg_dir.exists():
                # Fall back to dataset-level segs
                seg_dir = Path(f"data/{d}/downsampled/{res}px/segs_converted_square")
            
            # Check for overlap with split-specific GTs
            gt_dir = Path(f"data/{d}/{od}/downsampled/{res}px/GTs")
            
            if seg_dir.exists() and gt_dir.exists():
                seg_names = {p.stem for p in seg_dir.glob("*.png")}
                gt_names = {p.stem for p in gt_dir.glob("*.png")}
                if seg_names & gt_names:  # If there's overlap
                    has_overlap = True
                    break
        
        if has_overlap:
            valid_other_dirs.append(od)
        else:
            print(f"Skipping {d}/{od} for metrics: no overlap between segmentations and GTs")
    
    METRICS_AV_OTHERDIR_VALUES[d] = sorted(valid_other_dirs)

# --------------------------
# Helper functions
# --------------------------

def _gt_av_dir(dataset: str, res: str) -> str:
    """A/V GT dir for dataset+res (simple layout)."""
    ds = AV_GT_DIRS.get(dataset, {}) if isinstance(AV_GT_DIRS, dict) else {}
    if isinstance(ds, dict) and ds.get(res):
        return str(ds[res])
    return f"data/{dataset}/downsampled/{res}px/GTs"

def _gt_av_dir_other(dataset: str, other_dir: str, res: str) -> str:
    """A/V GT dir for dataset+other_dir+res (otherdir layout)."""
    return f"data/{dataset}/{other_dir}/downsampled/{res}px/GTs"

def _vessel_native_dirs(dataset: str):
    """Return (gt_native_dir, seg_native_dir) for vessel-only datasets."""
    ds = VESSEL_GT_DIRS.get(dataset, {}) if isinstance(VESSEL_GT_DIRS, dict) else {}
    if not isinstance(ds, dict):
        return ("", "")
    return (ds.get("gt_native_dir", ""), ds.get("seg_native_dir", ""))

def _metrics_inputs_simple(wc):
    """Dependencies for dataset-only rules."""
    ds = wc.dataset

    refined_dirs = [
        f"results/refined/{ds}/k{k}/downsampled/{res}px"
        for k in K_VALUES
        for res in RESOLUTIONS
    ]
    # CHANGED: Use *_square directories
    unref_dirs = [
        f"data/{ds}/downsampled/1024px/segs_converted_square",
        f"data/{ds}/downsampled/576px/segs_converted_square",
    ]
    roi_dirs = [
        f"data/{ds}/downsampled/1024px/roi_masks_binarized_square",
        f"data/{ds}/downsampled/576px/roi_masks_binarized_square",
    ]

    if ds in AV_GT_DATASETS:
        gt_dirs = [_gt_av_dir(ds, "1024"), _gt_av_dir(ds, "576")]
    else:
        gt_native_dir, seg_native_dir = _vessel_native_dirs(ds)
        gt_dirs = [gt_native_dir, seg_native_dir]

    return [p for p in (refined_dirs + unref_dirs + roi_dirs + gt_dirs) if p]

def _metrics_inputs_otherdir(wc):
    """Dependencies for dataset+other_dir rules."""
    ds = wc.dataset
    od = wc.other_dir

    refined_dirs = [
        f"results/refined/{ds}/k{k}/downsampled/{res}px"
        for k in K_VALUES
        for res in RESOLUTIONS
    ]
    # CHANGED: Use *_square directories
    unref_dirs = [
        f"data/{ds}/downsampled/1024px/segs_converted_square",
        f"data/{ds}/downsampled/576px/segs_converted_square",
    ]
    roi_dirs = [
        f"data/{ds}/downsampled/1024px/roi_masks_binarized_square",
        f"data/{ds}/downsampled/576px/roi_masks_binarized_square",
    ]
    gt_dirs = [_gt_av_dir_other(ds, od, "1024"), _gt_av_dir_other(ds, od, "576")]

    return [p for p in (refined_dirs + unref_dirs + roi_dirs + gt_dirs) if p]

# --------------------------
# Rules
# --------------------------

# Build wildcard constraints for other_dir (only valid splits)
_valid_other_dirs = [
    od
    for d in METRICS_AV_OTHERDIR
    for od in METRICS_AV_OTHERDIR_VALUES.get(d, [])
]
_other_dir_constraint = "|".join(map(re.escape, _valid_other_dirs)) if _valid_other_dirs else "NO_MATCH"

rule compute_metrics_av_simple:
    """Compute metrics for A/V datasets whose GTs are in the simple layout."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, METRICS_AV_SIMPLE)) if METRICS_AV_SIMPLE else "NO_MATCH"
    input:
        _metrics_inputs_simple
    output:
        csv_1024 = f"{METRICS_DIR}" + "/{dataset}/metrics_1024.csv",
        csv_576  = f"{METRICS_DIR}" + "/{dataset}/metrics_576.csv",
    params:
        has_av_gt = 1,
        refined_root   = lambda wc: f"results/refined/{wc.dataset}",
        # CHANGED: Use *_square directories
        unref_1024_dir = lambda wc: f"data/{wc.dataset}/downsampled/1024px/segs_converted_square",
        unref_576_dir  = lambda wc: f"data/{wc.dataset}/downsampled/576px/segs_converted_square",
        gt_1024_dir = lambda wc: _gt_av_dir(wc.dataset, "1024"),
        gt_576_dir  = lambda wc: _gt_av_dir(wc.dataset, "576"),
        # CHANGED: Use *_square directories
        roi_1024_dir = lambda wc: f"data/{wc.dataset}/downsampled/1024px/roi_masks_binarized_square",
        roi_576_dir  = lambda wc: f"data/{wc.dataset}/downsampled/576px/roi_masks_binarized_square",
        gt_native_dir  = "",
        seg_native_dir = "",
    script:
        "../scripts/compute_metrics_smk.py"

rule compute_metrics_av_otherdir:
    """Compute metrics for A/V datasets with multiple GT subfolders ({other_dir})."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, METRICS_AV_OTHERDIR)) if METRICS_AV_OTHERDIR else "NO_MATCH",
        other_dir=_other_dir_constraint
    input:
        _metrics_inputs_otherdir
    output:
        csv_1024 = f"{METRICS_DIR}" + "/{dataset}/{other_dir}/metrics_1024.csv",
        csv_576  = f"{METRICS_DIR}" + "/{dataset}/{other_dir}/metrics_576.csv",
    params:
        has_av_gt = 1,
        refined_root   = lambda wc: f"results/refined/{wc.dataset}",
        # CHANGED: Use *_square directories
        unref_1024_dir = lambda wc: f"data/{wc.dataset}/downsampled/1024px/segs_converted_square",
        unref_576_dir  = lambda wc: f"data/{wc.dataset}/downsampled/576px/segs_converted_square",
        gt_1024_dir = lambda wc: _gt_av_dir_other(wc.dataset, wc.other_dir, "1024"),
        gt_576_dir  = lambda wc: _gt_av_dir_other(wc.dataset, wc.other_dir, "576"),
        # CHANGED: Use *_square directories
        roi_1024_dir = lambda wc: f"data/{wc.dataset}/downsampled/1024px/roi_masks_binarized_square",
        roi_576_dir  = lambda wc: f"data/{wc.dataset}/downsampled/576px/roi_masks_binarized_square",
        gt_native_dir  = "",
        seg_native_dir = "",
    script:
        "../scripts/compute_metrics_smk.py"

rule compute_metrics_vessel_only:
    """Optional: compute metrics for vessel-only datasets (must be configured in metrics.vessel_gt_dirs)."""
    wildcard_constraints:
        dataset="|".join(map(re.escape, METRICS_VESSEL_ONLY)) if METRICS_VESSEL_ONLY else "NO_MATCH"
    input:
        _metrics_inputs_simple
    output:
        csv_1024 = f"{METRICS_DIR}" + "/{dataset}/metrics_1024.csv",
        csv_576  = f"{METRICS_DIR}" + "/{dataset}/metrics_576.csv",
    params:
        has_av_gt = 0,
        refined_root   = lambda wc: f"results/refined/{wc.dataset}",
        # CHANGED: Use *_square directories
        unref_1024_dir = lambda wc: f"data/{wc.dataset}/downsampled/1024px/segs_converted_square",
        unref_576_dir  = lambda wc: f"data/{wc.dataset}/downsampled/576px/segs_converted_square",
        gt_1024_dir = "",
        gt_576_dir  = "",
        roi_1024_dir = "",
        roi_576_dir  = "",
        gt_native_dir  = lambda wc: _vessel_native_dirs(wc.dataset)[0],
        seg_native_dir = lambda wc: _vessel_native_dirs(wc.dataset)[1],
    script:
        "../scripts/compute_metrics_smk.py"
