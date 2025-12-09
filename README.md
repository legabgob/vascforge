# Retinal Segmentation Relabelling & Refinement Pipeline

Snakemake-based pipeline for **refining retinal vessel artery/vein (A/V) labels**,  
computing quality metrics, and preparing datasets for downstream analysis.

At a high level:

1. Colour fundus images are segmented with an external model (e.g. **VascX**) to obtain initial vessel + A/V maps.
2. These maps are **refined** with the RRWNet refinement module for different numbers of iterations `k`.
3. The pipeline computes **Dice-based metrics** comparing:
   - unrefined vs ground truth,
   - refined vs ground truth,
   - unrefined vs refined predictions,
   using either full A/V ground truth (Fundus-AVSeg) or vessel-only ground truth (FIVES).
4. Utility rules handle **ROI masks, label conversion, downsampling, and simple pixel relabelling**.

---

## Features

- Snakemake workflow for:
  - running RRWNet refinement over `{dataset, resolution, k}` combinations,
  - computing DICE metrics in **A/V mode** (Fundus-AVSeg) or **vessel-only mode** (FIVES).
- Support for 2 image resolutions (currently **576 px** and **1024 px**).
- Dataset-aware behaviour:
  - Fundus-AVSeg: RGB A/V ground truth + ROI masks.
  - FIVES: grayscale vessel ground truth + vessel segmentations.
- Helpers:
  - convert grayscale label maps `{0,1,2,3}` → RGB A/V labels,
  - downsample images in batch,
  - generate circular ROI masks from `bounds.csv` **or** `meta.csv`,
  - replace grayscale pixel value `1` → `255` in a directory of PNGs.

---

## Repository structure

```text
workflow/
  snakefile              # main entrypoint
  rules/
    gray_to_rgb.smk      # grayscale A/V label → RGB conversion
    downsample.smk       # image / mask downsampling
    roi.smk              # ROI mask generation from bounds/meta CSV
    replace_1_to_255.smk # simple relabelling of grayscale masks
    refinement.smk       # RRWNet refinement (grid over k, res, dataset)
    metrics.smk          # unified metrics for Fundus-AVSeg & FIVES
scripts/
  gray2rgb_smk.py
  downsample_smk.py
  roi_from_csv_smk.py
  replace_1_to_255_smk.py
  compute_metrics_smk.py
```
