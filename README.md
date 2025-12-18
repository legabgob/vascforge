# Segmentation relabelling & refinement (Snakemake)

This repository contains a **Snakemake pipeline** to standardize inputs coming from `seg_legacy` exports, generate ROI masks, downsample images/masks, run **RRWNet refinement** over multiple settings, and compute evaluation metrics.

The workflow is designed to be **portable across machines**: you should only need to set **one root directory** where all datasets’ `seg_legacy` outputs live, plus a few dataset-specific options in `workflow/config/config.yaml`.

---

## What the pipeline does

The workflow is organized into stages:

1. **Convert grayscale A/V labels → RGB**  
   Converts per-image grayscale label maps (values like `{0,1,2,3}`) into RGB A/V/crossings encoding.
2. **Copy dataset metadata (`meta.csv` or `bounds.csv`) into the repo workspace**  
   Copies the per-dataset metadata file into `data/{dataset}/meta/` so subsequent rules read from a consistent location.
3. **Generate ROI masks**  
   Builds circular ROI masks from either `bounds.csv` (dict-like column) or `meta.csv` (explicit columns).
4. **Binarize masks (optional helper)**  
   Replaces pixel value `1 → 255` for binary masks.
5. **Downsample predictions and masks**  
   Downsamples directories to target widths (e.g. 576px and 1024px).
6. **RRWNet refinement**  
   Runs refinement for combinations of `{dataset, resolution, k}`.

> Important: **Snakemake decides execution order from dependencies**, not from the order of `include:` statements.  
> You get the intended order by making each stage’s outputs become the next stage’s inputs.

---

## Repo layout

Typical structure:

```text
workflow/
  Snakefile
  config/
    config.yaml
  rules/
    vascxgray_to_rgb.smk
    copy_meta.smk
    create_roi_masks.smk
    binarize_masks.smk
    downsample.smk
    refinement.smk
scripts/
  gray2rgb_smk.py
  roi_from_csv_smk.py
  binarize_masks_smk.py
  downsample_smk.py
  ...
data/
  weights/          # model weights (tracked via Git LFS)
  {dataset}/
    meta/
    roi_masks/
    ...
```

---

## Configuration

Edit:

```bash
nano workflow/config/config.yaml
```

A typical config contains:

- `seg_legacy_root`: **the only required absolute path** (root folder containing all datasets)
- `datasets`: list of datasets to run (e.g. `["Fundus-AVSeg", "FIVES"]`)
- `dataset_subdirs`: optional mapping for datasets that live under an extra subfolder
- `resolutions`: e.g. `[576, 1024]`
- `k_values`: e.g. `[3,4,5,6,7,8]`
- `weights`: paths to RRWNet weights (preferably **relative to the repo**, e.g. `data/weights/...`)

Example:

```yaml
seg_legacy_root: /HDD/data/relabelling-project

datasets:
  - Fundus-AVSeg
  - FIVES

# Some datasets are stored as: {dataset}/{other_dir}/seg_legacy/...
# Leave empty for datasets that are directly under {dataset}/seg_legacy/...
dataset_subdirs:
  Fundus-AVSeg: ""          # direct layout
  FIVES: "train"            # example of dataset/otherdir layout

resolutions: [576, 1024]
k_values: [3, 4, 5, 6, 7, 8]

weights:
  1024: data/weights/rrwnet_HRF_0.pth
  576:  data/weights/rrwnet_RITE_refinement.pth
```

### What is `refine_datasets` (if present in your config)?

If your config contains something like `refine_datasets`, it typically means:

- **datasets for which the refinement step should run**
- useful if you want to run preprocessing for many datasets but only refine a subset

If you don’t need that distinction, keep just one `datasets` list.

---

## Running the workflow

From the repo root:

### Dry-run (recommended first)

```bash
snakemake -n
```

### Run

```bash
snakemake -j 4
```

### Run a single rule (by output)

For wildcarded rules, target a **concrete output**. Example (directory output):

```bash
snakemake -j1 data/Fundus-AVSeg/downsampled/1024px/segs
```

## Notes on datasets with an extra subdirectory

Some datasets are stored like:

- `.../{dataset}/seg_legacy/...`
- others like `.../{dataset}/{other_dir}/seg_legacy/...`

The workflow supports both by using a dataset-specific `other_dir` from the config:

- if `other_dir` is empty → direct layout
- if set → use the nested layout

This is used consistently for:
- finding `av/{sample}.png` inputs
- locating `meta.csv` / `bounds.csv`
- any other dataset-local inputs derived from `seg_legacy`

---

## Git LFS for model weights

RRWNet `.pth` files are >100MB and **must** be stored with Git LFS.

Typical setup:

```bash
git lfs install
git lfs track "*.pth"
git add .gitattributes
git add data/weights/*.pth
git commit -m "Add weights via Git LFS"
git push origin main
```

