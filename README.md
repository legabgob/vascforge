```
 /$$    /$$                              /$$$$$$$$                                     
| $$   | $$                             | $$_____/                                     
| $$   | $$ /$$$$$$   /$$$$$$$  /$$$$$$$| $$     /$$$$$$   /$$$$$$   /$$$$$$   /$$$$$$ 
|  $$ / $$/|____  $$ /$$_____/ /$$_____/| $$$$$ /$$__  $$ /$$__  $$ /$$__  $$ /$$__  $$
 \  $$ $$/  /$$$$$$$|  $$$$$$ | $$      | $$__/| $$  \ $$| $$  \__/| $$  \ $$| $$$$$$$$
  \  $$$/  /$$__  $$ \____  $$| $$      | $$   | $$  | $$| $$      | $$  | $$| $$_____/
   \  $/  |  $$$$$$$ /$$$$$$$/|  $$$$$$$| $$   |  $$$$$$/| $$      |  $$$$$$$|  $$$$$$$
    \_/    \_______/|_______/  \_______/|__/    \______/ |__/       \____  $$ \_______/
                                                                    /$$  \ $$          
                                                                   |  $$$$$$/          
                                                                    \______/           

```
A Snakemake pipeline to **convert**, **mask/ROI**, **downsample**, **refine** retinal vessel segmentations (A/V/crossings) and optionally **evaluate** with DICE metrics.

This repo is designed to be **portable across machines**: you configure a single root directory containing your dataset exports (e.g. `seg_legacy/`) and the pipeline builds a consistent `data/` workspace and `results/` outputs.

---

## What this pipeline does

Given dataset outputs produced by your legacy tooling (typically `seg_legacy/`), the workflow can:

1. **Convert** grayscale A/V label maps → RGB label maps (A/V/crossings mapping).
2. **Copy metadata** (`meta.csv` or `bounds.csv`) into the workspace.
3. **Create ROI masks** from `meta.csv` or `bounds.csv` (circular ROI).
4. **Binarize ROI masks** (e.g. convert `1 → 255` if needed).
5. **Downsample** predicted labels, ROI masks (and optionally GTs) to target resolutions (e.g. 576px, 1024px).
6. **Refine** predictions using RRWNet refinement (`k` iterations) and write refined outputs.
7. **Compute DICE metrics** (optional) for:
   - datasets with A/V ground truth (RGB GT), and/or
   - datasets with vessel-only ground truth.

The workflow supports datasets that are organized either as:

- `ROOT/{dataset}/seg_legacy/...`
- `ROOT/{dataset}/{other_dir}/seg_legacy/...`  (multiple `other_dir` per dataset, kept separate)

---

## Repository layout

- `workflow/`
  - `Snakefile` – main entrypoint
  - `config/config.yaml` – all configuration (paths, datasets, refinement grid, metrics, GT handling)
  - `rules/` – modular rules (`vascxgray_to_rgb.smk`, `copy_meta.smk`, `create_roi_masks.smk`, `binarize_masks.smk`, `downsample.smk`, `refinement.smk`, `dice.smk`, etc.)
- `scripts/`
  - Snakemake-friendly scripts (directory-to-directory batch processing, ROI creation, metrics computation)
- `data/` *(generated)* – standardized workspace outputs
- `results/` *(generated)* – refined predictions + evaluation outputs

---

## Installation

### 1) Clone

```bash
git clone <your-repo-url>
cd segmentation-relabelling
```

### 2) Environment

You need Python 3.10+ and the core packages used by the scripts:

- `snakemake`
- `numpy`
- `pillow`
- `pandas` (for metrics)

If you use conda/mamba:

```bash
mamba create -n relabelling -c conda-forge python=3.10 snakemake pandas numpy pillow
mamba activate relabelling
```

---

## Configuration (`workflow/config/config.yaml`)

The pipeline is driven by `workflow/config/config.yaml`.

At minimum, set:

- `legacy_root`: the **root directory** where dataset exports live (contains the datasets folders).
- `datasets`: which datasets to process.
- `resolutions`: list of downsample widths.
- `k_range`: refinement iteration range.

Example skeleton:

```yaml
legacy_root: "/HDD/data/relabelling-project"

datasets:
  - Fundus-AVSeg
  - leuven-haifa
  - FIVES

resolutions: ["576", "1024"]
k_range: [3, 9]   # produces k=3..8
```

### Datasets with `{other_dir}`

Some datasets are structured as:

```
legacy_root/dataset/other_dir/seg_legacy/...
```

The workflow auto-detects these and preserves separation under:

```
data/{dataset}/{other_dir}/...
```

### Model weights

Refinement uses RRWNet weights. If you keep weights in-repo, use Git LFS (recommended). If you keep them outside the repo, point to them via config.

---

## Running the workflow

From the repo root:

### Dry run (recommended)

```bash
snakemake -n
```

### Print commands + reasons (still dry)

```bash
snakemake -n -p -r
```

### Run

```bash
snakemake -j 8
```

---

## Outputs

### Workspace outputs (`data/`)

Common outputs include:

- `data/{dataset}/segs_converted/`  
  RGB converted predictions (from grayscale labels)
- `data/{dataset}/meta/meta.csv`  
  Copied metadata (`meta.csv` or `bounds.csv`)
- `data/{dataset}/roi_masks/`  
  ROI masks generated from metadata
- `data/{dataset}/roi_masks_binarized/`  
  Clean/binarized ROI masks
- `data/{dataset}/downsampled/{res}px/...`  
  Downsampled predictions / masks / (optionally) GTs

If a dataset has multiple `other_dir`, the outputs become:

- `data/{dataset}/{other_dir}/...`

### Refined predictions (`results/`)

Refinement outputs:

- `results/refined/{dataset}/k{k}/downsampled/{res}px/`

Each `k` and resolution gets its own directory.

### Metrics (`results/metrics/`)

Metrics are written to:

- Simple layout:
  - `results/metrics/{dataset}/metrics_1024.csv`
  - `results/metrics/{dataset}/metrics_576.csv`
- `{other_dir}` layout (GTs separated by `other_dir`):
  - `results/metrics/{dataset}/{other_dir}/metrics_1024.csv`
  - `results/metrics/{dataset}/{other_dir}/metrics_576.csv`

---

## Ground truth (optional)

Some datasets have A/V ground truth (RGB GT), others only vessel ground truth (grayscale).

The workflow supports GT ingestion + conversion + downsampling when enabled in config:

```yaml
ground_truth:
  enabled: true
  # root: "/path/to/gt_root"  # defaults to legacy_root if omitted

  av_rgb:
    datasets:
      Fundus-AVSeg:
        pattern: "Fundus-AVSeg/GTs/{sample}.png"
        mapping:
          - [[255, 0, 0], [255, 0, 255]]
          - [[0, 0, 255], [0, 255, 255]]
          - [[0, 255, 0], [255, 255, 255]]
          - [[255, 255, 255], [0, 0, 255]]

      leuven-haifa:
        pattern: "leuven-haifa/{other_dir}/label/{sample}.png"
        mapping: [...]

  vessel_gray:
    datasets:
      FIVES:
        pattern: "FIVES/{other_dir}/GTs/{sample}.png"
```

- A/V GTs end up as `.../downsampled/{res}px/GTs`
- vessel-only GTs end up as `.../downsampled/{res}px/GTs_vessel`

---

## Common troubleshooting

### “Missing input files … GTs”
If a dataset’s GTs use `{other_dir}`, metrics must be computed per `(dataset, other_dir)`.
Make sure:
- your `ground_truth` patterns include `{other_dir}`
- `rule all` includes the correct targets
- your metrics rules are split into *simple* vs *other_dir* variants

### “The flag 'directory' is only valid for outputs”
In Snakemake, `directory()` must only be used for **outputs**, not inputs.
Inputs should be plain strings (or functions returning strings).

### Visualize the DAG
`snakemake --dag` prints a Graphviz DOT graph to stdout:

```bash
snakemake --dag | dot -Tpng > dag.png
```

