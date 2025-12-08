import os

# Datasets and what kind of GT they have
HAS_AV_GT = {
    "Fundus-AVSeg": 1,  # A/V RGB GT and ROI masks
    "FIVES": 0,         # vessel-only GT
}

# Refined prediction roots (per dataset)
REFINED_ROOT = {
    "Fundus-AVSeg": "/SSD/home/gabriel/rrwnet/refined_predictions/Fundus-AVSeg",
    "FIVES": "/SSD/home/gabriel/rrwnet/refined_predictions/FIVES",
}

# Unrefined RGB segs (per dataset & resolution)
UNREF_1024_DIR = {
    "Fundus-AVSeg": "/SSD/home/gabriel/rrwnet/data/Fundus-AVSeg/downsampled/1024px/segs",
    "FIVES": "/SSD/home/gabriel/rrwnet/data/FIVES/train/downsampled/1024px/segs",
}

UNREF_576_DIR = {
    "Fundus-AVSeg": "/SSD/home/gabriel/rrwnet/data/Fundus-AVSeg/downsampled/576px/segs",
    "FIVES": "/SSD/home/gabriel/rrwnet/data/FIVES/train/downsampled/576px/segs",
}

# A/V GT + ROI (Fundus-AVSeg style)
GT_1024_DIR = {
    "Fundus-AVSeg": "/SSD/home/gabriel/rrwnet/data/Fundus-AVSeg/downsampled/1024px/GTs",
}
GT_576_DIR = {
    "Fundus-AVSeg": "/SSD/home/gabriel/rrwnet/data/Fundus-AVSeg/downsampled/576px/GTs",
}
ROI_1024_DIR = {
    "Fundus-AVSeg": "/SSD/home/gabriel/rrwnet/data/Fundus-AVSeg/downsampled/1024px/masks",
}
ROI_576_DIR = {
    "Fundus-AVSeg": "/SSD/home/gabriel/rrwnet/data/Fundus-AVSeg/downsampled/576px/masks",
}

# Vessel-only GT (FIVES style)
GT_NATIVE_DIR = {
    "FIVES": "/SSD/home/gabriel/rrwnet/data/FIVES/train/downsampled/GTs",
}
SEG_NATIVE_DIR = {
    "FIVES": "/SSD/home/gabriel/rrwnet/data/FIVES/train/vessel_segs",
}

# List of datasets you want to evaluate
DATASETS = ["Fundus-AVSeg", "FIVES"]


rule compute_metrics:
    """
    Compute DICE metrics for refinement:
      - If dataset has A/V GT, use RGB GT + ROI (Fundus-AVSeg-style).
      - Otherwise use vessel-only GT (FIVES-style).
    """
    output:
        csv_1024 = "metrics/{dataset}_1024.csv",
        csv_576  = "metrics/{dataset}_576.csv",
    params:
        has_av_gt = lambda wc: HAS_AV_GT[wc.dataset],
        refined_root = lambda wc: REFINED_ROOT[wc.dataset],
        unref_1024_dir = lambda wc: UNREF_1024_DIR[wc.dataset],
        unref_576_dir = lambda wc: UNREF_576_DIR[wc.dataset],
        # A/V GT branch (Fundus-AVSeg). For datasets without AV GT, these keys can be missing/unused.
        gt_1024_dir = lambda wc: GT_1024_DIR.get(wc.dataset, ""),
        gt_576_dir  = lambda wc: GT_576_DIR.get(wc.dataset, ""),
        roi_1024_dir = lambda wc: ROI_1024_DIR.get(wc.dataset, ""),
        roi_576_dir  = lambda wc: ROI_576_DIR.get(wc.dataset, ""),
        # Vessel-only GT branch (FIVES)
        gt_native_dir = lambda wc: GT_NATIVE_DIR.get(wc.dataset, ""),
        seg_native_dir = lambda wc: SEG_NATIVE_DIR.get(wc.dataset, ""),
    script:
        "scripts/compute_metrics_smk.py"

