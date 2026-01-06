# workflow/rules/dice.smk
#
# Runs scripts/compute_metrics_smk.py to produce per-dataset CSVs.
#
# IMPORTANT:
#   This rule declares *real inputs* (refinement outputs + downsample outputs),
#   so Snakemake will run refinement/downsample first.
#
# No hard-coded /SSD/... paths: everything is either pipeline-produced or configurable.

# --------------------------
# Configuration (with safe defaults)
# --------------------------

METRICS_CFG = config.get("metrics", {}) if isinstance(config, dict) else {}

# Default: only Fundus-AVSeg and leuven-haifa are considered A/V-GT datasets (if present in config["datasets"])
_default_av = {"Fundus-AVSeg", "leuven-haifa"}
_datasets = config.get("datasets", [])
if isinstance(_datasets, dict):
    _datasets = list(_datasets.keys())
_datasets = set(_datasets)

AV_GT_DATASETS = set(METRICS_CFG.get("av_gt_datasets", list(_default_av & _datasets)))

# Which datasets to *actually compute metrics for*.
# Default: compute only for A/V-GT datasets, to avoid failing on datasets that don't have the needed GT inputs yet.
METRICS_DATASETS = METRICS_CFG.get("datasets", sorted(AV_GT_DATASETS))

# Where CSVs are written
METRICS_DIR = METRICS_CFG.get("out_dir", "results/metrics")

# Optional per-dataset overrides for GT locations (A/V GT).
# If not provided, we assume GTs are located at:
#   data/{dataset}/downsampled/{res}px/GTs
AV_GT_DIRS = METRICS_CFG.get("av_gt_dirs", {})

# Optional per-dataset overrides for vessel-only GT locations.
VESSEL_GT_DIRS = METRICS_CFG.get("vessel_gt_dirs", {})

# Refinement grid (used only to build dependency inputs)
RESOLUTIONS = [str(r) for r in config.get("resolutions", ["576", "1024"])]
k_start, k_end = config.get("k_range", [3, 9])
K_VALUES = list(range(int(k_start), int(k_end)))

# --------------------------
# Helpers
# --------------------------

def _has_av_gt(wc):
    return int(wc.dataset in AV_GT_DATASETS)

def _gt_av_dir(dataset: str, res: str) -> str:
    # 1) config override: metrics.av_gt_dirs[dataset][res]
    ds = AV_GT_DIRS.get(dataset, {}) if isinstance(AV_GT_DIRS, dict) else {}
    if isinstance(ds, dict) and ds.get(res):
        return str(ds[res])
    # 2) default expected location inside repo workspace
    return f"data/{dataset}/downsampled/{res}px/GTs"

def _vessel_native_dirs(dataset: str):
    # Must be provided via config if you want vessel-only mode to work.
    # Example structure:
    # metrics:
    #   vessel_gt_dirs:
    #     FIVES:
    #       gt_native_dir: /path/to/GTs
    #       seg_native_dir: /path/to/vessel_segs
    ds = VESSEL_GT_DIRS.get(dataset, {}) if isinstance(VESSEL_GT_DIRS, dict) else {}
    return (str(ds.get("gt_native_dir", "")), str(ds.get("seg_native_dir", "")))

def _metrics_inputs(wc):
    """Declare real dependencies so compute_metrics runs after refinement + downsample."""
    ds = wc.dataset

    # All refined dirs for this dataset (for all K and all RESOLUTIONS)
    refined_dirs = [
        f"results/refined/{ds}/k{k}/downsampled/{res}px"
        for k in K_VALUES
        for res in RESOLUTIONS
    ]

    # Unrefined RGB segs produced by downsample.smk
    unref_dirs = [
        f"data/{ds}/downsampled/1024px/segs_converted",
        f"data/{ds}/downsampled/576px/segs_converted",
    ]

    # ROI masks produced by downsample.smk (masked evaluation)
    roi_dirs = [
        f"data/{ds}/downsampled/1024px/roi_masks_binarized",
        f"data/{ds}/downsampled/576px/roi_masks_binarized",
    ]

    # Ground truth inputs
    if ds in AV_GT_DATASETS:
        gt_dirs = [_gt_av_dir(ds, "1024"), _gt_av_dir(ds, "576")]
    else:
        gt_native_dir, seg_native_dir = _vessel_native_dirs(ds)
        gt_dirs = [gt_native_dir, seg_native_dir]

    # Filter out empty strings (so non-configured vessel-only datasets don't create a weird empty input)
    return [p for p in (refined_dirs + unref_dirs + roi_dirs + gt_dirs) if p]

# --------------------------
# Rule
# --------------------------

rule compute_metrics:
    """
    Compute DICE metrics for refinement outputs.

    Ground truth selection:
      - A/V GT mode for datasets listed in metrics.av_gt_datasets
      - vessel-only mode otherwise (requires metrics.vessel_gt_dirs per dataset)
    """
    input:
        _metrics_inputs
    output:
        csv_1024 = f"{METRICS_DIR}" + "/{dataset}/metrics_1024.csv",
        csv_576  = f"{METRICS_DIR}" + "/{dataset}/metrics_576.csv",
    params:
        has_av_gt = _has_av_gt,

        # Refined root produced by refinement.smk
        refined_root = lambda wc: f"results/refined/{wc.dataset}",

        # Unrefined dirs produced by downsample.smk
        unref_1024_dir = lambda wc: f"data/{wc.dataset}/downsampled/1024px/segs_converted",
        unref_576_dir  = lambda wc: f"data/{wc.dataset}/downsampled/576px/segs_converted",

        # A/V GT dirs (configurable)
        gt_1024_dir = lambda wc: _gt_av_dir(wc.dataset, "1024"),
        gt_576_dir  = lambda wc: _gt_av_dir(wc.dataset, "576"),

        # ROI masks produced by downsample.smk (masked metrics)
        roi_1024_dir = lambda wc: f"data/{wc.dataset}/downsampled/1024px/roi_masks_binarized",
        roi_576_dir  = lambda wc: f"data/{wc.dataset}/downsampled/576px/roi_masks_binarized",

        # Vessel-only GT dirs (optional; only used if has_av_gt == 0)
        gt_native_dir  = lambda wc: _vessel_native_dirs(wc.dataset)[0],
        seg_native_dir = lambda wc: _vessel_native_dirs(wc.dataset)[1],
    script:
        "scripts/compute_metrics_smk.py"

