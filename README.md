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

# VascForge

A **Snakemake pipeline** for processing, refining, and evaluating retinal vessel segmentations. VascForge automates the workflow from raw segmentation outputs through refinement using RRWNet, and produces quantitative metrics and vascular features.

## Features

- **Format Conversion**: Convert grayscale A/V label maps (0/1/2/3) to RGB label maps (arteries/veins/crossings)
- **ROI Masking**: Create and apply circular ROI masks from metadata CSV files
- **Downsampling**: Resize predictions and masks to target resolutions (e.g., 576px, 1024px)
- **Refinement**: Apply RRWNet iterative refinement with configurable k iterations
- **DICE Metrics**: Compute per-class and macro-average DICE scores against ground truth
- **Feature Extraction**: Extract vascular biomarkers using the VascX library (tortuosity, caliber, bifurcations, etc.)

## Repository Structure

```
vascforge/
├── workflow/
│   ├── Snakefile              # Main workflow entrypoint
│   ├── config/
│   │   └── config.yaml        # Pipeline configuration
│   ├── envs/
│   │   └── environment.yaml   # Conda environment specification
│   ├── rules/                 # Modular Snakemake rules
│   │   ├── vascxgray_to_rgb.smk
│   │   ├── copy_meta.smk
│   │   ├── create_roi_masks.smk
│   │   ├── binarize_masks.smk
│   │   ├── downsample.smk
│   │   ├── refinement.smk
│   │   ├── dice.smk
│   │   ├── gt_conversion.smk
│   │   └── ...
│   └── scripts/               # Python scripts called by rules
│       ├── gray_to_rgb_dir_smk.py
│       ├── downsample_dir_smk.py
│       ├── compute_metrics_smk.py
│       ├── run_full_pipeline.py
│       └── ...
├── data/
│   └── weights/               # Pre-trained RRWNet model weights
└── README.md
```

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/legabgob/vascforge.git
cd vascforge
```

### 2. Create the Environment

Using conda or mamba:

```bash
mamba create -n vascforge -c conda-forge python=3.10 snakemake pandas numpy pillow
mamba activate vascforge
```

Or use the provided environment file:

```bash
mamba env create -f workflow/envs/environment.yaml
mamba activate vascforge
```

### 3. Configure External Dependencies

The pipeline requires:
- **RRWNet**: Set the path to `get_predictions.py` in `workflow/config/config.yaml`
- **VascX** (optional): For vascular feature extraction

## Configuration

Edit `workflow/config/config.yaml` to configure the pipeline:

```yaml
# Root directory containing your datasets
legacy_root: "/path/to/your/data"

# RRWNet configuration
rrwnet:
  script: "/path/to/rrwnet/get_predictions.py"
  python: "python3"

# Model weights (shipped with repo)
weights:
  "1024": "data/weights/rrwnet_HRF_0.pth"
  "576": "data/weights/rrwnet_RITE_refinement.pth"

# Datasets to process
datasets:
  - "Fundus-AVSeg"
  - "FIVES"
  - "leuven-haifa"

# Target resolutions
resolutions: ["576", "1024"]

# Refinement iterations (k=3 through k=8)
k_range: [3, 9]
```

### Dataset Structure

VascForge expects datasets organized as:

```
legacy_root/
├── Dataset-Name/
│   └── seg_legacy/
│       ├── original/     # Original fundus images
│       ├── av/           # A/V segmentation masks
│       ├── ce/           # Contrast-enhanced images
│       ├── rgb/          # RGB preprocessed images
│       ├── discs/        # Optic disc masks
│       ├── fovea.csv     # Fovea coordinates
│       └── meta.csv      # Image metadata (bounds, centers, radii)
```

For datasets with train/test/val splits:

```
legacy_root/
├── Dataset-Name/
│   ├── train/
│   │   └── seg_legacy/...
│   ├── test/
│   │   └── seg_legacy/...
│   └── val/
│       └── seg_legacy/...
```

## Usage

### Dry Run (Recommended First)

Preview what the pipeline will do:

```bash
snakemake -n
```

Show commands and reasons:

```bash
snakemake -n -p -r
```

### Run the Pipeline

Execute with 8 parallel jobs:

```bash
snakemake -j 8
```

### Visualize the DAG

Generate a workflow graph:

```bash
snakemake --dag | dot -Tpng > workflow_dag.png
```

## Outputs

### Workspace (`data/`)

- `data/{dataset}/segs_converted/` — RGB converted predictions
- `data/{dataset}/meta/meta.csv` — Copied metadata
- `data/{dataset}/roi_masks/` — Generated ROI masks
- `data/{dataset}/roi_masks_binarized/` — Binarized ROI masks (0/255)
- `data/{dataset}/downsampled/{res}px/` — Downsampled predictions, masks, and GTs

### Refined Predictions (`results/refined/`)

- `results/refined/{dataset}/k{k}/downsampled/{res}px/` — Refined segmentations for each k value

### Metrics (`results/metrics/`)

- `results/metrics/{dataset}/metrics_{res}.csv` — DICE scores per image

### Vascular Features (`results/vascx_features_*/`)

- `results/vascx_features_refined/{dataset}/...` — Features from refined segmentations
- `results/vascx_features_unrefined/{dataset}/...` — Features from original segmentations

## Ground Truth Configuration

Enable ground truth processing in config:

```yaml
ground_truth:
  enabled: true
  
  av_rgb:
    datasets:
      Fundus-AVSeg:
        pattern: "Fundus-AVSeg/annotation/{sample}.png"
        mapping:
          - [[255, 0, 0], [255, 0, 255]]     # Red → Magenta (arteries)
          - [[0, 0, 255], [0, 255, 255]]     # Blue → Cyan (veins)
          - [[0, 255, 0], [255, 255, 255]]   # Green → White (crossings)
          - [[255, 255, 255], [0, 0, 255]]   # White → Blue (vessels)

  vessel_gray:
    datasets:
      FIVES:
        pattern: "FIVES/{other_dir}/ground_truth/{sample}.png"
```

## Troubleshooting

### "Missing input files … GTs"

Ensure your ground truth patterns include `{other_dir}` if the dataset uses train/test/val splits.

### "The flag 'directory' is only valid for outputs"

In Snakemake, `directory()` must only be used for **outputs**, not inputs.

### Permission Issues with Model Weights

The weights in `data/weights/` should be readable. If using Git LFS, run:

```bash
git lfs pull
```

## Citation

If you use VascForge in your research, please cite:

```bibtex
@software{vascforge,
  title = {VascForge: A Pipeline for Retinal Vessel Segmentation Processing},
  year = {2025},
  url = {https://github.com/legabgob/vascforge}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
